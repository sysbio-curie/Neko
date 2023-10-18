from __future__ import annotations

import collections

from pypath_common import misc as common

"""
Handling nodes of interest for network lookup queries.
"""

class Noi:

    def __init__(
            self,
            noi: Noi | list[str] | list[list[str]] | dict[str, list[str]],
            groups: list[str] | dict[str, str],
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

        self._setup(noi, groups)


    def _setup(
            self,
            noi: Noi | list[str] | list[list[str]] | dict[str, list[str]],
            groups: list[str] | dict[str, str],
        ):

        self.nodes = self._parse(noi, groups)


    @staticmethod
    def _parse(
            noi: Noi | list[str] | list[list[str]] | dict[str, list[str]],
            groups: list[str] | dict[str, str],
        ) -> dict[str, list[str]]:

        # later, we can convert str to omniflow.core._node.Node instances
        # this will enable to have non-protein nodes, multiple organisms, etc

        if isinstance(noi, Noi):

            return noi.nodes

        if isinstance(noi, dict):

            return noi

        noi = common.to_list(noi)

        if all(isinstance(n, str) for it in noi):

            if isinstance(dict, groups):

                groups = [groups[n] for n in noi]

            if isinstance(list, groups):

                if len(noi) != len(groups):

                    raise ValueError('Nodes and groups are not the same length')

                _noi = collections.defaultdict(list)
                [_noi[g].append(n) for n, g in zip(noi, groups))]

                return dict(_noi)

        return {groups or 'noi': noi}




