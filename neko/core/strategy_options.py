"""Validated options and legacy migration for connection strategies."""

from __future__ import annotations

from dataclasses import dataclass
import warnings

from typing_extensions import Literal


PathPolicy = Literal[
    "one_shortest",
    "all_shortest",
    "all_bounded",
]
ReusePolicy = Literal[
    "none",
    "discovered_paths",
    "induced_subgraph",
]

PATH_POLICIES = frozenset({
    "one_shortest",
    "all_shortest",
    "all_bounded",
})
REUSE_POLICIES = frozenset({
    "none",
    "discovered_paths",
    "induced_subgraph",
})

DEFAULT_PATH_POLICY: PathPolicy = "all_bounded"
DEFAULT_REUSE_POLICY: ReusePolicy = "discovered_paths"
LEGACY_BFS_NONE_CUTOFF = 10


class _UnsetType:
    def __repr__(self) -> str:
        return "UNSET"

    def __copy__(self):
        return self

    def __deepcopy__(self, memo):
        return self


UNSET = _UnsetType()


class ConnectionStrategyMigrationWarning(FutureWarning):
    """Visible warning for the transition to explicit strategy policies."""


@dataclass(frozen=True)
class ResolvedConnectionPolicies:
    maxlen: int
    path_policy: PathPolicy
    reuse_policy: ReusePolicy
    used_legacy_arguments: bool


def _validate_maxlen(maxlen) -> int:
    if isinstance(maxlen, bool) or not isinstance(maxlen, int) or maxlen <= 0:
        raise ValueError(
            "complete_connection requires maxlen to be a positive integer.",
        )
    return maxlen


def resolve_connection_policies(
        *,
        maxlen,
        path_policy: PathPolicy | None,
        reuse_policy: ReusePolicy | None,
        algorithm=UNSET,
        minimal=UNSET,
        connect_with_bias=UNSET,
        warning_stacklevel: int = 3,
    ) -> ResolvedConnectionPolicies:
    """Resolve new policies or one complete legacy argument combination."""

    legacy_arguments = (algorithm, minimal, connect_with_bias)
    used_legacy_arguments = any(value is not UNSET for value in legacy_arguments)
    used_new_arguments = path_policy is not None or reuse_policy is not None

    if used_legacy_arguments and used_new_arguments:
        raise ValueError(
            "Do not mix legacy algorithm/minimal/connect_with_bias arguments "
            "with path_policy/reuse_policy.",
        )

    if used_legacy_arguments:
        resolved_algorithm = "dfs" if algorithm is UNSET else algorithm
        resolved_minimal = True if minimal is UNSET else minimal
        resolved_bias = False if connect_with_bias is UNSET else connect_with_bias

        if resolved_algorithm not in {"bfs", "dfs"}:
            raise ValueError("algorithm must be either 'bfs' or 'dfs'.")
        if not isinstance(resolved_minimal, bool):
            raise TypeError("minimal must be a boolean.")
        if not isinstance(resolved_bias, bool):
            raise TypeError("connect_with_bias must be a boolean.")

        resolved_path = (
            "one_shortest"
            if resolved_algorithm == "bfs"
            else "all_bounded"
        )
        if resolved_bias:
            resolved_reuse = "induced_subgraph"
        elif resolved_minimal:
            resolved_reuse = "discovered_paths"
        else:
            resolved_reuse = "none"

        if maxlen is None and resolved_algorithm == "bfs":
            resolved_maxlen = LEGACY_BFS_NONE_CUTOFF
        else:
            resolved_maxlen = _validate_maxlen(maxlen)

        warnings.warn(
            "Legacy complete_connection parameters are deprecated. Use "
            f"maxlen={resolved_maxlen}, path_policy='{resolved_path}', "
            f"reuse_policy='{resolved_reuse}'.",
            ConnectionStrategyMigrationWarning,
            stacklevel=warning_stacklevel,
        )
        return ResolvedConnectionPolicies(
            maxlen=resolved_maxlen,
            path_policy=resolved_path,
            reuse_policy=resolved_reuse,
            used_legacy_arguments=True,
        )

    resolved_path = path_policy or DEFAULT_PATH_POLICY
    resolved_reuse = reuse_policy or DEFAULT_REUSE_POLICY

    if resolved_path not in PATH_POLICIES:
        choices = ", ".join(sorted(PATH_POLICIES))
        raise ValueError(f"path_policy must be one of: {choices}.")
    if resolved_reuse not in REUSE_POLICIES:
        choices = ", ".join(sorted(REUSE_POLICIES))
        raise ValueError(f"reuse_policy must be one of: {choices}.")

    return ResolvedConnectionPolicies(
        maxlen=_validate_maxlen(maxlen),
        path_policy=resolved_path,
        reuse_policy=resolved_reuse,
        used_legacy_arguments=False,
    )
