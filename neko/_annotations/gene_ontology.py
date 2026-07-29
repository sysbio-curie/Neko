import re
import logging
import os
import uuid
from dataclasses import dataclass

import omnipath as op
import requests
from requests.adapters import HTTPAdapter
from urllib3.util import Retry

import pandas as pd

from neko._cache import cache_dir


_HPA_TISSUE_RESOURCE = "HPA_tissue"
_REQUIRED_ANNOTATION_COLUMNS = frozenset({
    "genesymbol",
    "label",
    "record_id",
    "value",
})
_NOT_DETECTED_LEVELS = frozenset({"not detected"})
_HPA_CANCER_DATA_URL = (
    "https://www.proteinatlas.org/download/tsv/cancer_data.tsv.zip"
)
_HPA_CANCER_CACHE_SUBDIR = "annotations"
_HPA_CANCER_CACHE_FILENAME = "hpa_cancer_data.tsv.zip"
_MIN_HPA_CANCER_ROWS = 100_000
_HPA_CANCER_COLUMNS = frozenset({
    "Gene name",
    "Cancer",
    "High",
    "Medium",
    "Low",
    "Not detected",
})
_HPA_CANCER_TISSUES = frozenset({
    "breast cancer",
    "carcinoid",
    "cervical cancer",
    "colorectal cancer",
    "endometrial cancer",
    "glioma",
    "head and neck cancer",
    "liver cancer",
    "lung cancer",
    "lymphoma",
    "melanoma",
    "ovarian cancer",
    "pancreatic cancer",
    "prostate cancer",
    "renal cancer",
    "skin cancer",
    "stomach cancer",
    "testis cancer",
    "thyroid cancer",
    "urothelial cancer",
})

logger = logging.getLogger(__name__)

_GO_API_BASE = "https://api.geneontology.org/api"
_GO_ID_PATTERN = re.compile(r"^GO:\d{7}$")
_IEA_EVIDENCE = "ECO:0000501"
_GO_PAGE_SIZE = 100
_GO_MAX_PAGES = 1000


class AnnotationServiceError(RuntimeError):
    """Raised when tissue annotations cannot be retrieved or interpreted."""


class GeneOntologyError(RuntimeError):
    """Raised when GO data cannot be retrieved or interpreted."""


class GeneOntologyNotFoundError(GeneOntologyError):
    """Raised when a GO accession is not recognized by the GO API."""


@dataclass(frozen=True)
class GOTerm:
    """Canonical metadata for a Gene Ontology term."""

    go_id: str
    label: str


@dataclass(frozen=True)
class GOGene:
    """A gene returned by the GO association API."""

    gene_id: str | None
    symbol: str | None
    taxon_id: str
    taxon_label: str | None = None


def _normalize_annotation_value(value):
    """Normalize an OmniPath annotation value for comparison."""
    if pd.isna(value):
        return ""
    return " ".join(str(value).casefold().split())


def _hpa_cancer_cache_path():
    """Return the managed cache path for the HPA cancer table."""
    return (
        cache_dir(_HPA_CANCER_CACHE_SUBDIR)
        / _HPA_CANCER_CACHE_FILENAME
    )


def _new_hpa_session():
    """Create a bounded-retry session for Human Protein Atlas downloads."""
    session = requests.Session()
    retries = Retry(
        total=3,
        backoff_factor=0.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=frozenset({'GET'}),
    )
    session.mount('https://', HTTPAdapter(max_retries=retries))
    return session


def _read_hpa_cancer_table(path):
    """Read and validate the cached HPA cancer expression table."""
    table = pd.read_csv(path, sep='\t', compression='zip', low_memory=False)
    missing_columns = _HPA_CANCER_COLUMNS.difference(table.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"HPA cancer data is missing columns: {missing}.")
    if len(table) < _MIN_HPA_CANCER_ROWS:
        raise ValueError(
            "HPA cancer data is unexpectedly small "
            f"({len(table)} rows)."
        )
    return table


def _download_hpa_cancer_table(timeout=120):
    """Download, validate, and atomically cache the HPA cancer table."""
    path = _hpa_cancer_cache_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    token = f"{os.getpid()}.{uuid.uuid4().hex}"
    temporary_path = path.with_name(f".{path.name}.{token}.tmp")

    try:
        response = _new_hpa_session().get(
            _HPA_CANCER_DATA_URL,
            timeout=(10, timeout),
        )
        response.raise_for_status()
        temporary_path.write_bytes(response.content)
        table = _read_hpa_cancer_table(temporary_path)
        temporary_path.replace(path)
        return table
    finally:
        temporary_path.unlink(missing_ok=True)


def _load_hpa_cancer_table():
    """Load a valid cached HPA table, downloading it only when necessary."""
    path = _hpa_cancer_cache_path()
    if path.is_file():
        try:
            return _read_hpa_cancer_table(path)
        except (OSError, ValueError, pd.errors.ParserError) as exc:
            logger.warning("Ignoring invalid cached HPA cancer data: %s", exc)
    return _download_hpa_cancer_table()


def _hpa_cancer_expressed_genes(gene_symbols, tissue):
    """Return genes detected by HPA IHC in at least one cancer sample."""
    target_tissue = _normalize_annotation_value(tissue)
    table = _load_hpa_cancer_table()
    cancers = table['Cancer'].map(_normalize_annotation_value)
    symbols = table['Gene name'].map(_normalize_annotation_value)
    requested = {
        _normalize_annotation_value(symbol)
        for symbol in gene_symbols
    }
    matches = table[cancers.eq(target_tissue) & symbols.isin(requested)].copy()
    if matches.empty:
        return set()

    detected_counts = matches[['High', 'Medium', 'Low']].apply(
        pd.to_numeric,
        errors='coerce',
    ).fillna(0).sum(axis=1)
    return set(symbols.loc[matches.index[detected_counts.gt(0)]])


def _expressed_genes(annotations_df, tissue):
    """Return normalized symbols expressed in ``tissue`` in HPA records."""
    if not isinstance(annotations_df, pd.DataFrame):
        raise AnnotationServiceError(
            "OmniPath returned an invalid tissue annotation response: "
            f"expected a pandas DataFrame, got {type(annotations_df).__name__}."
        )

    missing_columns = _REQUIRED_ANNOTATION_COLUMNS.difference(
        annotations_df.columns
    )
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise AnnotationServiceError(
            "OmniPath returned an invalid tissue annotation response; "
            f"missing columns: {missing}."
        )

    if annotations_df.empty:
        return set()

    records = annotations_df.copy()
    records["_label"] = records["label"].map(_normalize_annotation_value)
    records = records[records["_label"].isin(("tissue", "level"))]
    if records.empty:
        return set()

    records = (
        records.pivot_table(
            index=["genesymbol", "record_id"],
            columns="_label",
            values="value",
            aggfunc="first",
            observed=True,
        )
        .reset_index()
    )
    if "tissue" not in records or "level" not in records:
        return set()

    target_tissue = _normalize_annotation_value(tissue)
    tissues = records["tissue"].map(_normalize_annotation_value)
    levels = records["level"].map(_normalize_annotation_value)
    detected = levels.ne("") & ~levels.isin(_NOT_DETECTED_LEVELS)
    matching_records = records[tissues.eq(target_tissue) & detected]

    return {
        _normalize_annotation_value(symbol)
        for symbol in matching_records["genesymbol"]
        if not pd.isna(symbol)
    }


def _new_go_session(user_agent):
    """Create a retrying HTTP session for the official GO API."""
    session = requests.Session()
    retries = Retry(
        total=5,
        connect=3,
        read=3,
        status=5,
        backoff_factor=1.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({'GET'}),
        respect_retry_after_header=True,
        raise_on_status=False,
    )
    session.headers.update({
        'Accept': 'application/json',
        'User-Agent': user_agent,
    })
    session.mount('https://', HTTPAdapter(max_retries=retries))
    return session


class Ontology:
    """
    class that stores some functionalities to connect phenotypes to nodes
    and to associate information for each node at tissue level

    """
    def __init__(
            self,
            taxon_id=9606,
            timeout=30.0,
            user_agent=(
                "NeKo (https://github.com/sysbio-curie/Neko)"
            ),
            session=None,
        ):
        self.taxon_id = self._normalize_taxon_id(taxon_id)
        self.timeout = timeout
        self.session = session or _new_go_session(user_agent)
        self._term_cache = {}
        self.accession_to_phenotype_dict = {
            "GO:0010718": (
                "positive regulation of epithelial to mesenchymal transition"
            ),
        }

    @staticmethod
    def _normalize_go_id(go_id):
        """Return a canonical GO CURIE or raise for invalid input."""
        if go_id is None:
            raise ValueError("A GO accession is required.")
        go_id = str(go_id).strip().upper()
        if re.fullmatch(r"\d{7}", go_id):
            go_id = f"GO:{go_id}"
        if not _GO_ID_PATTERN.fullmatch(go_id):
            raise ValueError(
                f"Invalid GO accession {go_id!r}; expected GO: followed "
                "by seven digits."
            )
        return go_id

    @staticmethod
    def _normalize_taxon_id(taxon_id):
        """Return an NCBITaxon CURIE accepted by the GO API."""
        if isinstance(taxon_id, int):
            taxon_id = str(taxon_id)
        if not isinstance(taxon_id, str):
            raise ValueError("taxon_id must be an NCBI taxonomy identifier.")
        taxon_id = taxon_id.strip()
        if taxon_id.casefold().startswith("ncbitaxon:"):
            taxon_id = taxon_id.split(":", 1)[1]
        if not taxon_id.isdigit() or int(taxon_id) <= 0:
            raise ValueError("taxon_id must be a positive NCBI taxonomy ID.")
        return f"NCBITaxon:{taxon_id}"

    def _request_json(self, url, go_id, params=None):
        """Request and decode one GO API JSON response."""
        try:
            response = self.session.get(
                url,
                params=params,
                timeout=self.timeout,
            )
            if response.status_code == 404:
                raise GeneOntologyNotFoundError(
                    f"The GO API does not recognize accession {go_id}."
                )
            response.raise_for_status()
            return response.json()
        except GeneOntologyNotFoundError:
            raise
        except requests.RequestException as exc:
            raise GeneOntologyError(
                f"GO API request failed for {go_id}: {exc}"
            ) from exc
        except ValueError as exc:
            raise GeneOntologyError(
                f"GO API returned invalid JSON for {go_id}."
            ) from exc

    def get_term(self, go_id):
        """Return canonical metadata for a GO accession."""
        go_id = self._normalize_go_id(go_id)
        if go_id in self._term_cache:
            return self._term_cache[go_id]

        payload = self._request_json(
            f"{_GO_API_BASE}/ontology/term/{go_id}",
            go_id,
        )
        if not isinstance(payload, dict):
            raise GeneOntologyError(
                "Unexpected GO term response: expected a JSON object."
            )

        returned_id = payload.get("goid")
        label = payload.get("label")
        try:
            returned_id = self._normalize_go_id(returned_id)
        except ValueError as exc:
            raise GeneOntologyError(
                f"GO term response for {go_id} has no valid accession."
            ) from exc
        if not isinstance(label, str) or not label.strip():
            raise GeneOntologyError(
                f"GO term response for {go_id} has no valid label."
            )

        term = GOTerm(returned_id, label.strip())
        self._term_cache[go_id] = term
        self._term_cache[returned_id] = term
        self.accession_to_phenotype_dict[returned_id] = term.label
        return term

    @staticmethod
    def _association_evidence_ids(association):
        """Collect evidence CURIEs from a GO association."""
        evidence_ids = set()
        evidence = association.get("evidence")
        if isinstance(evidence, str):
            evidence_ids.add(evidence)
        for item in association.get("evidence_types") or []:
            if isinstance(item, dict) and isinstance(item.get("id"), str):
                evidence_ids.add(item["id"])
        return evidence_ids

    def fetch_go_genes(
            self,
            go_id,
            *,
            taxon_id=None,
            include_descendants=False,
            exclude_automatic_assertions=False,
            page_size=_GO_PAGE_SIZE,
            max_pages=_GO_MAX_PAGES,
        ):
        """Fetch unique genes associated with a GO term."""
        term = self.get_term(go_id)
        taxon_id = self.taxon_id if taxon_id is None else (
            self._normalize_taxon_id(taxon_id)
        )
        if not isinstance(page_size, int) or page_size <= 0:
            raise ValueError("page_size must be a positive integer.")
        if not isinstance(max_pages, int) or max_pages <= 0:
            raise ValueError("max_pages must be a positive integer.")
        page_size = min(page_size, _GO_PAGE_SIZE)

        url = (
            f"{_GO_API_BASE}/bioentity/function/"
            f"{term.go_id}/genes"
        )
        start = 0
        page_count = 0
        previous_signature = None
        genes_by_id = {}

        while page_count < max_pages:
            payload = self._request_json(
                url,
                term.go_id,
                params={
                    "start": start,
                    "rows": page_size,
                    "taxon": taxon_id,
                    "relationship_type": "involved_in",
                },
            )
            if not isinstance(payload, dict):
                raise GeneOntologyError(
                    "Unexpected GO association response: expected a JSON "
                    "object."
                )
            associations = payload.get("associations")
            if not isinstance(associations, list):
                raise GeneOntologyError(
                    "Unexpected GO association response: 'associations' "
                    "must be a list."
                )
            if not associations:
                break

            signature = tuple(
                (
                    association.get("id")
                    or (
                        (association.get("subject") or {}).get("id"),
                        (association.get("object") or {}).get("id"),
                        association.get("evidence"),
                        association.get("negated"),
                    )
                )
                if isinstance(association, dict) else None
                for association in associations
            )
            if signature == previous_signature:
                raise GeneOntologyError(
                    "GO API pagination did not advance."
                )
            previous_signature = signature

            for association in associations:
                if not isinstance(association, dict):
                    raise GeneOntologyError(
                        "GO API returned a malformed association."
                    )
                if association.get("negated", False):
                    continue
                if (
                    exclude_automatic_assertions
                    and _IEA_EVIDENCE in self._association_evidence_ids(
                        association
                    )
                ):
                    continue

                subject = association.get("subject")
                obj = association.get("object")
                if not isinstance(subject, dict) or not isinstance(obj, dict):
                    continue
                if not include_descendants and obj.get("id") != term.go_id:
                    continue

                taxon = subject.get("taxon")
                if not isinstance(taxon, dict) or taxon.get("id") != taxon_id:
                    continue
                gene_id = subject.get("id")
                symbol = subject.get("label")
                gene_id = gene_id if isinstance(gene_id, str) else None
                symbol = symbol if isinstance(symbol, str) else None
                if not gene_id and not symbol:
                    continue

                key = gene_id or f"{taxon_id}:{symbol}"
                gene = GOGene(
                    gene_id=gene_id,
                    symbol=symbol,
                    taxon_id=taxon_id,
                    taxon_label=(
                        taxon.get("label")
                        if isinstance(taxon.get("label"), str)
                        else None
                    ),
                )
                existing = genes_by_id.get(key)
                if existing is None or (
                    existing.symbol is None and gene.symbol is not None
                ):
                    genes_by_id[key] = gene

            page_count += 1
            start += len(associations)
            if len(associations) < page_size:
                break
        else:
            raise GeneOntologyError(
                f"GO API pagination exceeded {max_pages} pages for "
                f"{term.go_id}."
            )

        return sorted(
            genes_by_id.values(),
            key=lambda gene: (
                gene.symbol is None,
                (gene.symbol or gene.gene_id or "").casefold(),
            ),
        )

    def get_markers(
            self,
            phenotype=None,
            id_accession=None,
            *,
            taxon_id=None,
            include_descendants=False,
            exclude_automatic_assertions=False,
        ):
        """Return gene symbols associated with a GO term."""
        id_accession = self.resolve_accession(
            phenotype=phenotype,
            id_accession=id_accession,
        )
        genes = self.fetch_go_genes(
            id_accession,
            taxon_id=taxon_id,
            include_descendants=include_descendants,
            exclude_automatic_assertions=exclude_automatic_assertions,
        )
        return sorted({gene.symbol for gene in genes if gene.symbol})

    def resolve_accession(self, phenotype=None, id_accession=None):
        """Resolve an explicit accession or a registered phenotype alias."""
        if id_accession is not None:
            return self._normalize_go_id(id_accession)
        if phenotype is None:
            raise ValueError(
                "Provide at least one of id_accession or phenotype."
            )
        normalized_phenotype = _normalize_annotation_value(phenotype)
        matches = [
            accession
            for accession, description
            in self.accession_to_phenotype_dict.items()
            if _normalize_annotation_value(description)
            == normalized_phenotype
        ]
        if not matches:
            raise ValueError(
                "No locally registered GO accession was found for "
                f"phenotype {phenotype!r}."
            )
        id_accession = matches[0]
        return self._normalize_go_id(id_accession)

    def check_tissue_annotations(self, genes_df, tissue):
        """
        Check whether genes have detected HPA expression in a tissue.

        Args:
        genes_df (DataFrame): DataFrame containing gene symbols.
        tissue (str): Tissue to match exactly after case/whitespace normalization.

        Returns:
        DataFrame: Gene symbols and their detected-expression status.

        Raises:
        AnnotationServiceError: If OmniPath is unavailable or returns an
            unexpected annotation schema.
        """

        if not isinstance(genes_df, pd.DataFrame):
            raise TypeError(
                "genes_df must be a pandas DataFrame containing a "
                "'Genesymbol' column."
            )
        if 'Genesymbol' not in genes_df:
            raise ValueError("genes_df must contain a 'Genesymbol' column.")
        if not isinstance(tissue, str) or not tissue.strip():
            raise ValueError("tissue must be a non-empty string.")

        raw_symbols = genes_df['Genesymbol']
        if raw_symbols.isna().any():
            raise ValueError("genes_df contains a missing gene symbol.")

        gene_symbols = raw_symbols.astype(str).tolist()
        if any(not symbol.strip() for symbol in gene_symbols):
            raise ValueError("genes_df contains an empty gene symbol.")
        if not gene_symbols:
            return pd.DataFrame({
                'Genesymbol': pd.Series(dtype='object'),
                'in_tissue': pd.Series(dtype='bool'),
            })

        unique_symbols = list(dict.fromkeys(gene_symbols))
        normalized_tissue = _normalize_annotation_value(tissue)

        if normalized_tissue in _HPA_CANCER_TISSUES:
            try:
                expressed = _hpa_cancer_expressed_genes(
                    unique_symbols,
                    tissue,
                )
            except Exception as exc:
                raise AnnotationServiceError(
                    "Unable to retrieve cancer expression data from the "
                    "Human Protein Atlas. No genes were classified as absent."
                ) from exc
        else:
            # A single resource-restricted request is faster, cacheable as one
            # unit, and less likely to fail than one request per gene.
            try:
                annotations_df = op.requests.Annotations.get(
                    proteins=unique_symbols,
                    resources=_HPA_TISSUE_RESOURCE,
                )
            except Exception as exc:
                raise AnnotationServiceError(
                    "Unable to retrieve HPA tissue annotations from OmniPath. "
                    "The service may be temporarily unavailable; no genes "
                    "were classified as absent."
                ) from exc
            expressed = _expressed_genes(annotations_df, tissue)

        return pd.DataFrame({
            'Genesymbol': gene_symbols,
            'in_tissue': [
                _normalize_annotation_value(symbol) in expressed
                for symbol in gene_symbols
            ],
        })
