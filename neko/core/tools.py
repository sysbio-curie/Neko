from __future__ import annotations

import re

import networkx as nx

import pandas as pd

from neko.inputs import identifier_mapping as mapping
from neko.inputs import chebi_mapping

NON_PROTEIN_ENTITY_PREFIXES = (
    "CID:",
    "COMPLEX_NAME:",
    "PROTEIN_FAMILY:",
    "PHENOTYPE:",
    "SIGNOR-",
    "STIMULUS:",
    "URS",
)
PHOSPHOSITE_PATTERN = re.compile(
    r'^(?P<protein>\S+)_(?P<residue>[STY])(?P<position>\d+)$',
    re.IGNORECASE,
)


def normalize_phosphosite_identifier(node: str | None) -> str | None:
    """Return the canonical ``GENE_RESIDUE`` phosphosite identifier."""

    if not isinstance(node, str):
        return None

    match = PHOSPHOSITE_PATTERN.fullmatch(node.strip())

    if match is None:
        return None

    return (
        f'{match.group("protein")}_'
        f'{match.group("residue").upper()}'
        f'{match.group("position")}'
    )


def is_connected(network) -> bool:
    """
    This function checks if a network is connected. It takes a Network object as input and returns True if the network
    is connected, otherwise it returns False.

    Args:
        - network: A Network object representing the network to be checked.

    Returns:
        - bool
    """
    valid_edges = network.edges.dropna(subset=['source', 'target'])
    g = nx.from_pandas_edgelist(valid_edges, 'source', 'target')
    node_identifiers = network.nodes['Uniprot'].where(
        network.nodes['Uniprot'].notna(),
        network.nodes['Genesymbol'],
    )
    g.add_nodes_from(node_identifiers.dropna().unique())

    # NetworkX deliberately rejects the null graph. For NeKo's purposes an
    # empty or one-node network has no disconnected pair to resolve.
    return len(g) <= 1 or nx.is_connected(g)

def check_sign(interaction: pd.DataFrame, consensus: bool = False) -> str:
    """
    This function checks the sign of an interaction in the Omnipath format (Pandas DataFrame or Series).
    The attribute "consensus" checks for the consistency of the sign of the interaction among the references.

    Args:
        - interaction: A pandas DataFrame or Series representing the interaction.
        - consensus: A boolean indicating whether to check for consensus among references.

    Returns:
        - A string indicating the sign of the interaction: "stimulation", "inhibition", "form complex", or "undefined".
    """
    # Handle both DataFrame and Series input
    if isinstance(interaction, pd.DataFrame):
        interaction = interaction.iloc[0]

    if consensus:
        if interaction.get("consensus_inhibition") and interaction.get("consensus_stimulation"):
            return "bimodal"
        if interaction.get("consensus_stimulation"):
            return "stimulation"
        elif interaction.get("consensus_inhibition"):
            return "inhibition"
        else:
            return "undefined"
    else:
        # Check if it is both stimulation and inhibition
        if interaction.get("is_stimulation", True) and interaction.get("is_inhibition", True):
            return "bimodal"
        if interaction.get("is_stimulation", False):
            return "stimulation"
        elif interaction.get("is_inhibition", False):
            return "inhibition"
        # Check for "form_complex" column existence
        elif interaction.get("form_complex", False):
            return "form complex"
        else:
            return "undefined"


def check_gene_list_format(gene_list: list[str]) -> bool:
    """
    This function checks the format of the gene list and returns True if the gene list is in Uniprot format,
    False if the gene list is in genesymbol format.

    Args:
        - gene_list: A list of gene identifiers. The gene identifiers can be either Uniprot identifiers or genesymbols.

    Returns:
        - A boolean indicating whether the gene list is in Uniprot format (True) or genesymbol format (False).
    """
    def canonical_resource_identifier(gene):
        return (
            normalize_phosphosite_identifier(gene)
            or chebi_mapping.normalize_identifier(gene)
            or (
                isinstance(gene, str)
                and gene.startswith(
                    NON_PROTEIN_ENTITY_PREFIXES + ("COMPLEX:",)
                )
            )
        )

    # Canonical resource identifiers behave like UniProt accessions in graph
    # operations: they are already ready for edge lookup and must not be sent
    # to a protein identifier service.
    if all(
        canonical_resource_identifier(gene) or mapping.to_uniprot(gene)
        for gene in gene_list
    ):
        return True
    # Check if the gene list contains genesymbols
    elif all(mapping.to_genesymbol(gene) for gene in gene_list):
        return False


def mapping_node_identifier(node: str) -> list[str]:
    """
    This function takes a node identifier and returns a list containing the possible identifiers for the node.
    The identifiers include a complex string, a genesymbol, and a uniprot identifier. The function uses the
    to_uniprot and to_genesymbol functions from the neko.inputs.identifier_mapping module to translate the node
    identifier into these different formats.

    Args:
        - node: A string representing the node identifier. The node identifier can be a genesymbol, a uniprot identifier,
                or a complex string.

    Returns:
        - A list containing the complex string, genesymbol, and uniprot identifier for the node. If the node identifier
          cannot be translated into one of these formats, the corresponding value in the list is None.
    """
    complex_string = None
    genesymbol = None
    uniprot = None

    canonical_chebi = chebi_mapping.normalize_identifier(node)

    if canonical_chebi:
        # ChEBI accessions are canonical non-protein identifiers. A cached
        # ASCII name is display-only; failure to obtain it must never change
        # the identifier used by resource edges or trigger a UniProt lookup.
        genesymbol = chebi_mapping.to_name(canonical_chebi) or canonical_chebi

        return [complex_string, genesymbol, canonical_chebi]

    if isinstance(node, str) and node.startswith(NON_PROTEIN_ENTITY_PREFIXES):
        # Typed group/context nodes are already normalized identifiers. Keep
        # the identifier intact for resource-edge matching and display it as
        # the node label without attempting a meaningless UniProt lookup.
        return [complex_string, node, node]

    if isinstance(node, str) and node.startswith("COMPLEX:"):
        # Check the complex prefix before attempting any translation:
        # "COMPLEX:X_Y" is never itself a valid gene symbol or UniProt
        # accession, so translating it directly would always miss (wasting
        # a live lookup); translate its individual members instead.
        node_content = node[8:]
        node_list = node_content.split("_")

        # Translate each element in node_list using mapping.to_genesymbol
        translated_node_list = [
            mapping.to_genesymbol(mapping.to_uniprot(item)) or item
            for item in node_list
        ]

        # Join the elements in node_list with "_"
        joined_node_string = "_".join(translated_node_list)

        # Add back the "COMPLEX:" prefix to the string
        complex_string = "COMPLEX:" + joined_node_string
        uniprot = node

        return [complex_string, genesymbol, uniprot]

    phosphosite = normalize_phosphosite_identifier(node)

    if phosphosite:
        # PhosphoSitePlus uses identifiers such as MAP3K4_T1494 directly in
        # resource edges. They are canonical site nodes, not gene symbols to
        # submit to UniProt.
        return [complex_string, phosphosite, phosphosite]

    node_id = mapping.to_uniprot(node)

    if node_id:
        # Convert UniProt ID to gene symbol
        uniprot = node_id
        if uniprot.startswith("MI"):
            genesymbol = uniprot
        else:
            genesymbol = mapping.to_genesymbol(uniprot)
    else:
        label = mapping.to_genesymbol(node)
        if label:
            genesymbol = label
            uniprot = mapping.to_uniprot(genesymbol)
        else:
            print("Error during translation, check syntax for ", node)

    return [complex_string, genesymbol, uniprot]


def translate_paths(paths) -> list[list[str]]:
    """
    This function translates a list of paths, where each path is a sequence of node identifiers.
    It uses the helper function `handle_complex_identifier` to translate each node identifier in the paths.

    Args:
        - paths: A list of paths, where each path is a sequence of node identifiers.
                 A node identifier can be a string or a list of strings.

    Returns:
        - A list of translated paths, where each path is a sequence of translated node identifiers.
    """
    translated_list = []

    def handle_complex_identifier(item):
        """
        This helper function translates a node identifier using the `mapping_node_identifier` function.
        It checks all possible identifiers (complex, genesymbol, uniprot) and returns the first non-None value.

        Args:
        - item: A node identifier.

        Returns:
        - The translated node identifier.
        """
        identifiers = mapping_node_identifier(item)
        return identifiers[0] or identifiers[1] or identifiers[2]

    # If input_list is a list of strings
    if isinstance(paths[0], str):
        translated_list = [handle_complex_identifier(item) for item in paths]
    # If input_list is a list of lists of strings
    elif isinstance(paths[0], list):
        for sublist in paths:
            translated_sublist = [handle_complex_identifier(item) for item in sublist]
            translated_list.append(translated_sublist)

    return translated_list


def _is_missing(value) -> bool:
    """Return whether a scalar edge value should be treated as missing."""

    if value is None or value is False:
        return True

    if not pd.api.types.is_scalar(value):
        return False
    return bool(pd.isna(value))


def _join_unique_values(values, separator: str = '; ') -> str | None:
    """Join non-empty values once while preserving their input order."""

    unique_values = []

    for value in values:
        if _is_missing(value):
            continue

        # References and types are already merged with semicolons elsewhere in
        # the package. Split them here so repeated consolidation remains
        # idempotent.
        parts = str(value).split(';')
        for part in parts:
            part = part.strip()
            if part and part not in unique_values:
                unique_values.append(part)

    return separator.join(unique_values) if unique_values else None


def join_unique(series) -> str:
    """
    This function takes a pandas Series, filters out None values, and returns a string of unique values joined by a comma.

    Args:
        - series: A pandas Series object.

    Returns: - A string of unique values in the series, joined by a comma. If a value in the series is None,
                it is not included in the output string.
    """
    return _join_unique_values(series, separator=', ') or ''


def _normalize_effect(effect) -> str:
    """Normalize effect spelling before parallel edges are consolidated."""

    if _is_missing(effect):
        return 'undefined'

    normalized = str(effect).strip().lower().replace('_', ' ')
    aliases = {
        'activate': 'stimulation',
        'activation': 'stimulation',
        'stimulate': 'stimulation',
        'inhibit': 'inhibition',
        'both': 'bimodal',
        'form complex': 'form complex',
    }
    return aliases.get(normalized, normalized or 'undefined')


def consolidate_edges(edges: pd.DataFrame) -> pd.DataFrame:
    """Merge parallel working edges without discarding conflicting signs.

    Regulatory evidence is consolidated per source-target pair. Observing both
    stimulation and inhibition produces one ``bimodal`` edge, regardless of
    evidence counts. Complex formation and non-regulatory effects remain
    separate because they cannot be represented by a regulatory sign.
    """

    if not isinstance(edges, pd.DataFrame):
        raise TypeError('edges must be a pandas DataFrame.')

    required = {'source', 'target', 'Effect'}
    missing_columns = required.difference(edges.columns)
    if missing_columns:
        detail = ', '.join(sorted(missing_columns))
        raise ValueError(f'Edges are missing required column(s): {detail}.')

    columns = ['source', 'target', 'Type', 'Effect', 'References']
    if edges.empty:
        return edges.reindex(columns=columns).copy().reset_index(drop=True)

    working = edges.copy()
    for column in ('Type', 'References'):
        if column not in working.columns:
            working[column] = None

    working['_normalized_effect'] = working['Effect'].map(_normalize_effect)
    consolidated = []
    regulatory_effects = {'stimulation', 'inhibition', 'bimodal'}

    for (source, target), group in working.groupby(
            ['source', 'target'], sort=False, dropna=False):
        categories = []
        for effect in group['_normalized_effect']:
            if effect in regulatory_effects:
                category = 'regulatory'
            elif effect == 'form complex':
                category = 'form complex'
            else:
                category = effect
            if category not in categories:
                categories.append(category)

        for category in categories:
            if category == 'regulatory':
                selected = group[
                    group['_normalized_effect'].isin(regulatory_effects)
                ]
                effects = set(selected['_normalized_effect'])
                if 'bimodal' in effects or {
                        'stimulation', 'inhibition'}.issubset(effects):
                    effect = 'bimodal'
                elif 'stimulation' in effects:
                    effect = 'stimulation'
                else:
                    effect = 'inhibition'
            else:
                selected = group[group['_normalized_effect'] == category]
                effect = category

            consolidated.append({
                'source': source,
                'target': target,
                'Type': _join_unique_values(selected['Type']),
                'Effect': effect,
                'References': _join_unique_values(selected['References']),
            })

    return pd.DataFrame(consolidated, columns=columns).reset_index(drop=True)


def determine_most_frequent_effect(effects):
    effect_counts = effects.value_counts()

    # If there's only one type of effect, return it
    if len(effect_counts) == 1:
        return effect_counts.index[0]

    # Count the occurrences of each effect type
    stimulation_count = effect_counts.get('stimulation', 0)
    inhibition_count = effect_counts.get('inhibition', 0)
    form_complex_count = effect_counts.get('form_complex', 0)

    # Calculate the total of stimulation, inhibition, and form_complex
    total_known_effects = stimulation_count + inhibition_count + form_complex_count

    # If there are more unknown effects than known ones, return 'undefined'
    if len(effects) - total_known_effects > total_known_effects:
        return 'undefined'

    # If form_complex is the most frequent, return it
    if form_complex_count > stimulation_count and form_complex_count > inhibition_count:
        return 'form_complex'

    # Handle stimulation and inhibition
    if stimulation_count > inhibition_count:
        return 'stimulation'
    elif inhibition_count > stimulation_count:
        return 'inhibition'
    elif stimulation_count == inhibition_count and stimulation_count > 0:
        return 'bimodal'

    # If we've reached this point, it means there's no clear majority
    return 'undefined'
