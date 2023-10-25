from __future__ import annotations

from typing import Literal, Sequence
import re

import pypath.utils.taxonomy as _taxonomy
import pypath.utils.mapping as _mapping
import pypath.utils.orthology as _orthology
import pypath.inputs.uniprot as _uniprot
import pypath.inputs.uniprot_db as _uniprot_db

"""
Procedures on Nodes of molecular interaction networks.
"""

_REENS = re.compile(r'ENS[A-Z]*?([A-Z])[0-9]+')
_RERSQ = re.compile(
    r'^((?:(WP|AC|AP|NC|NG|NM|NP|NR|NT|NW|XM|XP|XR|YP|ZP)_)|(?:NZ\_[A-Z]{2,4}))'
    r'(\d+)\.(\d+)?$'
)
_REMI = re.compile(r'^MI:(\d+)$')
_REMM = re.compile(r'^MIMAT:(\d+)$')
_REMIR = re.compile(r'(?i)[a-z]{3}-(?:mir|let)-.+')
_RECPX = re.compile(r'^COMPLEX:.+')
_RECBL = re.compile(r'CHEMBL:\d+')
_RECBI = re.compile(r'CHEBI:\d+')
_SMOL_TYPES = {
    'small_molecule',
    'compound',
    'drug',
    'metabolite',
    'lipid',
}


class Node:

    def __init__(
            self,
            identifier: str,
            id_type: str | None = None,
            entity_type: str | None = None,
            organism: int | str = 'human',
            original_id: str | None = None,
            original_id_type: str | None = None,
        ):
        """
        Node of a molecular interaction network to be used in queries.

        Args:
            identifier:
                Identifier of the node.
            id_type:
                Type of the identifier. If `None`, a guess will be attempted.
            entity_type:
                Type of the molecular entity.
            organism:
                Name or NCBI Taxonomy ID of the organism.
        """

        self.identifier = identifier
        self._id_type = id_type
        self.entity_type = entity_type
        self.organism = _taxonomy.ensure_ncbi_taxid(organism)

        self._bootstrap()


    @property
    def id_type(self):

        return self._id_type


    @id_type.setter
    def id_type(self, id_type: str):

        self._id_type = id_type
        self.original_id_type = self.original_id_type or id_type


    def _bootstrap(self):

        self._guess_entity_type()
        self._guess_id_type()


    def _guess_entity_type(self):

        if self.entity_type is None:

            self._guess_by_id(
                (_REMI, _REMM, _REMIR, _RECPX, _RECBL, _RECBI),
                ('mir-pre', 'mirbase', 'mir-name', None, 'chembl', 'chebi'),
                ('mirna',) * 3 + ('complex',) + ('small_molecule',) * 2,
            )

            # fall back to default
            self.entity_type = self.entity_type or 'protein'


    def _guess_id_type(self):

        if self.id_type is None:

            entity_type = 'smol' if self._is_smol else
            getattr(self, f'_guess_{entity_type}_id_type')()


    @property
    def _is_smol(self):

        return self.entity_type in _SMOL_TYPES


    def _guess_protein_id_type(self):

        if (
            self.identifier.startswith('ENS') and
            (m := _REENS.match(self.identifier))
        ):
            self.id_type = f'ens{m.groups()[0].lower()}'

        elif m := _RERSQ.match(self.identifier):
            p = (m.groups()[0][1] == 'P') * 'p'
            self.id_type = f'refseq{p}'

        elif self.identifier.isdigit():
            # best guess, it could be also GI
            self.id_type = 'pubchem' if self._is_smol else 'entrez'

        elif (
            _uniprot.valid_uniprot(self.identifier) and
            _uniprot_db.is_uniprot(self.identifier, self.organism)
        ):
            self.id_type = 'uniprot'


    def _guess_mirna_id_type(self):

        self._guess_by_id(
            (_REMI, _REMM, _REMIR),
            ('mir-pre', 'mirbase', 'mir-name'),
        )


    def _guess_smol_id_type(self):

        if self.identifier.isdigit():

            self.id_type = 'pubchem'

        else:

            self._guess_by_id((_RECBL, _RECBI), ('chembl', 'chebi'))


    def _guess_by_id(
            self,
            regex: Sequence[re.compile],
            id_types: Sequence[str],
            entity_types: = Sequence[str] | None = None,
        ):

        entity_types = entity_types or (self.entity_type,) * len(regex)

        for rex, it, et in zip(regex, id_types, entity_types):

            if rex.match(self.identifier):

                self.id_type = it
                self.entity_type = et
                break
        for rex, it in zip((_RECBL, _RECBI)), ('chembl', 'chebi')):

                if rex.match(self.identifier):

                    self.id_type = it
                    break


    def as_idtype(self, id_type: str) -> Generator[Node, None, None]:
        """
        Translate to another identifier type.

        Args:
            id_type:
                The target identifier type.

        Returns:
            Due to potential ambiguous translation, a generator of nodes is
            returned.
        """

        for _id in _mapping.map_name(
            self.identifier,
            self.id_type,
            id_type,
            entity_type = self.entity_type,
            ncbi_tax_id = self.organism,
        ):

            yield Node(
                _id,
                id_type = id_type,
                organism = self.organism,
                entity_type = self.entity_type,
                original_id = self.original_id,
                original_id_type = self.original_id_type,
            )


    def as_organism(
            self,
            organism: str | int = 10090,
        """
        Translate to another organism by orthologous gene pairs.

        Args:
            organism:
                Name or NCBI Taxonomy ID of the target organism.

        Returns:
            Due to potential ambiguous translation, a generator of nodes is
            returned.
        """
        ) -> Generator[Node, None, None]:

        organsim = taxonomy.ensure_ncbi_taxid(organism)

        for _id in _orthology.translate(
            self.identifier,
            source = self.organism,
            target = organism,
            id_type = self.id_type,
        ):

            yield Node(
                _id,
                id_type = self.id_type,
                entity_type = self.entity_type,
                organsim = organism,
                original_id = self.original_id,
                original_id_type = self.original_id_type,
            )
