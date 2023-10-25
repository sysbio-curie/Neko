from __future__ import annotations

from typing import Any
import collections

from pypath_common import misc as common

"""
Handling nodes of interest for network lookup queries.
"""

class Noi:

    def __init__(
            self,
            noi: Noi | list[str] | list[list[str]] | dict[str, list[str]],
            groups: list[str] | dict[str, str] | None = None,
            id_type: str | None = 'uniprot',
            entity_type: str | None = 'protein',
            organism: int | str | None = 9606,
        ):
        """
        Nodes of interest.

        Args:
            noi:
                Nodes of interest, accepted notations are as follows:
                    - One string: one single node
                    - Iterable of strings: multiple nodes, grouping can be
                      specified by ``groups``
                    - Iterable of interables of strings: multiple groups of
                      nodes, group names can be specified by ``groups``
        """

        self._setup(noi, groups, id_type, entity_type, organism)


    def _setup(
            self,
            noi: Noi | list[str] | list[list[str]] | dict[str, list[str]],
            groups: list[str] | dict[str, str] | None = None,
            id_type: str | None = 'uniprot',
            entity_type: str | None = 'protein',
            organism: int | str | None = 9606,
        ):

        self.nodes = self._parse(noi, groups, id_type, entity_type, organism)


    @staticmethod
    def _parse(
            noi: Noi | list[str] | list[list[str]] | dict[str, list[str]],
            groups: list[str] | dict[str, str] | None = None,
            id_type: str | None = 'uniprot',
            entity_type: str | None = 'protein',
            organism: int | str | None = 9606,
        ) -> dict[str, list[str]]:

        if isinstance(noi, Noi):

            return noi.nodes

        if isinstance(noi, dict):

            groups, noi = zip(*(
                (grp, n)
                for grp, nodes in noi.items()
                for n in nodes
            ))
            groups = list(groups)

        noi = common.to_list(noi)

        def wrong_node(n: Any):

            raise ValueError(
                f'Nodes must be Node, string or tuple, got `{n}`.',
            )

        noi = [
            n
                if isinstance(n, Node) else
            Node(
                n,
                id_type = id_type,
                entity_type = entity_type,
                organism = organism,
            )
                if isinstance(n, str) else
            Node(*(
                from_node or from_arg
                for from_node, from_arg in
                zip(
                    n + (None,) * 3,
                    (None, id_type, entity_type, organism),
                )
            ))
                if isinstance(n, tuple) else
            wrong_node(n)
            for n in noi
        ]

        if isinstance(groups, dict):

            groups = [groups[n] for n in noi]

        if isinstance(groups, list):

            if len(noi) != len(groups):

                raise ValueError('Nodes and groups are not the same length')

            _noi = collections.defaultdict(list)
            [_noi[g].append(n) for n, g in zip(noi, groups))]

            groups = dict(_noi)

        return groups or {'noi': noi}




