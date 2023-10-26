from __future__ import annotations

from typing import Any, Hashable
import collections

from pypath_common import misc as _common

"""
Handling nodes of interest for network lookup queries.
"""

class Noi(collections.abc.Mapping):

    _DEFAULTS = {
        'id_type': 'uniprot',
        'entity_type': 'protein',
        'organism': 9606,
    }

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

        self._defaults = {
            'id_type': id_type,
            'entity_type': entity_type,
            'organism': organism,
        }
        self._setup(noi, groups)


    def _setup(
            self,
            noi: Noi | list[str] | list[list[str]] | dict[str, list[str]],
            groups: list[str] | dict[str, str] | None = None,
        ):

        self.nodes = self._parse(noi, groups)


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

        noi = [self._parse_node(n, locals()) for n in _common.to_list(noi)]

        if isinstance(groups, dict):

            groups = [groups[n] for n in noi]

        if isinstance(groups, list):

            if len(noi) != len(groups):

                raise ValueError('Nodes and groups are not the same length')

            _noi = collections.defaultdict(list)
            [_noi[g].append(n) for n, g in zip(noi, groups))]

            groups = dict(_noi)

        return groups or {'noi': noi}


    def _parse_node(
            self,
            node: str | Node | tuple | dict,
            attrs: dict,
        ) -> Node:

        if not isinstance(node, Node):

            _attrs = self._defaults.copy()
            _attrs.update({
                k: v
                for k, v in (attrs or {}).items()
                if k in self._defaults and v
            })
            _attrs.update(
                node
                    if isinstance(node, dict) else
                {k: v for k, v in zip(Node._attrs, node) if v}
                    if isinstance(node, tuple) else
                {'identifier': node}
                    if isinstance(node, str) else
                {}
            )

            if 'identifier' not in _attrs:

                raise ValueError(
                    'Nodes must have identifiers. '
                    f'Could not find identifier for `{node}`.',
                )

            node = Node(**_attrs)

        return node


    def __getitem__(self, key: Hashable) -> Any:

        return self.nodes[key]


    def __iter__(self):

        return self.nodes.__iter__()


    def __len__(self) -> int:

        return len(self.nodes)


    def __delitem__(self, key: Hashable):

        del self.nodes[key]


    def __setitem__(self, key: Hashable, value: Any):

        self.nodes[key] = value


    def as_idtype(self, id_type: str, entity_type: str = 'protein') -> Noi:

        return Noi(
            {
                grp:
                for grp, nodes in self.nodes.items()
            }
        )


