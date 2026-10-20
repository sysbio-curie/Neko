"""Offline history/growth benchmark, not a timing-sensitive CI test.

Run in separate processes with the same interpreter and --source checkouts.
Network construction is excluded from strategy timings. Provider discovery and
identifier mapping are offline; these synthetic timings are not production SLAs.
"""

import argparse
from contextlib import redirect_stdout
import io
import hashlib
import json
import os
from pathlib import Path
import statistics
import sys
from tempfile import TemporaryDirectory
import time
from unittest.mock import patch


def main():
    # Strategy frontiers include sets. Fix their iteration order across the
    # fresh processes used for paired timing and exact history comparisons.
    if os.environ.get("PYTHONHASHSEED") != "0":
        os.environ["PYTHONHASHSEED"] = "0"
        os.execv(sys.executable, [sys.executable, *sys.argv])
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--nodes", type=int, default=100)
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--skip-strategies", action="store_true")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(Path(args.source).resolve()))
    with TemporaryDirectory() as directory:
        config = Path(directory) / ".config" / "omnipathdb.ini"
        config.parent.mkdir()
        config.write_text(
            "[https://omnipathdb.org]\nautoload = false\n"
            "cache_dir = memory\nlicense = None\n"
        )
        with patch.object(Path, "home", return_value=Path(directory)):
            import neko.core.network as library
            from neko.core import history
            from neko._outputs.exports import Exports
    import pandas as pd

    library.mapping_node_identifier = lambda value: [None, value, value]
    library.check_gene_list_format = lambda values: False
    pairs = [
        (f"N{i}", f"N{(i + gap) % args.nodes}")
        for i in range(args.nodes)
        for gap in (1, 5)
    ]
    resources = pd.DataFrame(
        [
            {
                "source": a,
                "target": b,
                "is_directed": True,
                "is_stimulation": True,
                "is_inhibition": False,
                "form_complex": False,
                "references": "123456",
            }
            for a, b in pairs
        ]
    )
    seeds = [f"N{i}" for i in range(0, args.nodes, 10)]
    all_nodes = [f"N{i}" for i in range(args.nodes)]
    report = {
        "source": args.source,
        "nodes": args.nodes,
        "steps": args.steps,
        "repeats": args.repeats,
        "hash_seed": 0,
        "strategies": {},
    }

    def build(initial):
        return library.Network(initial_nodes=initial, resources=resources)

    def measure(operation):
        samples = []
        for _ in range(args.repeats):
            start = time.perf_counter()
            operation()
            samples.append((time.perf_counter() - start) * 1000)
        return {"median_ms": statistics.median(samples), "samples_ms": samples}

    operations = {
        "connect_nodes": (all_nodes, lambda n: n.connect_nodes()),
        "complete_connection": (
            seeds,
            lambda n: n.complete_connection(
                maxlen=10,
                path_policy="all_shortest",
                reuse_policy="discovered_paths",
            ),
        ),
        "radial": (
            seeds,
            lambda n: n.connect_network_radially(
                max_len=3,
                direction="OUT",
            ),
        ),
        "upstream": (
            seeds,
            lambda n: n.connect_to_upstream_nodes(
                nodes_to_connect=seeds,
                depth=3,
                rank=3,
            ),
        ),
    }
    if args.skip_strategies:
        operations = {}
    for name, (initial, operation) in operations.items():
        modes = {}
        edge_signatures = []
        for enabled in (True, False):
            samples = []
            for _ in range(args.repeats):
                network = build(initial)
                network.set_history_tracking(enabled)
                with redirect_stdout(io.StringIO()):
                    start = time.perf_counter()
                    operation(network)
                    samples.append((time.perf_counter() - start) * 1000)
            modes["history_on" if enabled else "history_off"] = {
                "median_ms": statistics.median(samples),
                "samples_ms": samples,
                "result_nodes": len(network.nodes),
                "result_edges": len(network.edges),
                "edges_sha256": hashlib.sha256(
                    history.canonical(history.records(network.edges)).encode()
                ).hexdigest(),
            }
            if enabled:
                modes["history_on"]["history_sha256"] = hashlib.sha256(
                    history.canonical(network.history_dict()).encode()
                ).hexdigest()
            edge_signatures.append(history.records(network.edges))
        assert edge_signatures[0] == edge_signatures[1]
        network = build(initial)
        with patch.object(
            history, "frame_fingerprint", wraps=history.frame_fingerprint
        ) as fingerprints:
            with patch.object(
                network, "save_state", wraps=network.save_state
            ) as saves:
                with redirect_stdout(io.StringIO()):
                    operation(network)
                modes["fingerprint_calls"] = fingerprints.call_count
                modes["snapshots_created"] = saves.call_count
        report["strategies"][name] = modes

    with TemporaryDirectory() as directory:
        sif = Path(directory) / "fixture.sif"
        sif.write_text("".join(f"{a}\t1\t{b}\n" for a, b in pairs))
        network = library.Network(sif_file=str(sif), resources=resources)
        for index in range(args.steps):
            network.modify_node_name(
                "N0" if index % 2 == 0 else "Renamed",
                "Renamed" if index % 2 == 0 else "N0",
            )
        previous = "N0" if args.steps % 2 == 0 else "Renamed"

        def rename():
            nonlocal previous
            new = "Renamed" if previous == "N0" else "N0"
            network.modify_node_name(previous, new)
            previous = new

        report["label_mutation"] = measure(rename)
        report["history_dict"] = measure(network.history_dict)
        report["full_list"] = measure(network.list_states)
        if hasattr(network, "history_status"):
            report["history_status"] = measure(network.history_status)
            report["summary_page"] = measure(
                lambda: network.list_states(offset=2, limit=5, detail="summary")
            )
        with redirect_stdout(io.StringIO()):
            report["sif_export"] = measure(
                lambda: Exports(network).export_sif(
                    Path(directory) / "export.sif"
                )
            )
        report["companion_bytes"] = (
            (Path(directory) / "export.sif.neko.json").stat().st_size
        )
        report["retained_states"] = len(network._state_log)
    Path(args.output).write_text(json.dumps(report, indent=2) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
