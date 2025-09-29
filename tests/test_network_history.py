import importlib.util
import sys
import types
from importlib.machinery import ModuleSpec
from pathlib import Path

import pandas as pd
import pytest

# Minimal pypath stub to satisfy imports in the core module
pypath_module = types.ModuleType("pypath")
pypath_module.__spec__ = ModuleSpec("pypath", loader=None)
pypath_module.__file__ = __file__

pypath_utils = types.ModuleType("pypath.utils")
pypath_utils.__spec__ = ModuleSpec("pypath.utils", loader=None)
pypath_utils.__file__ = __file__

pypath_mapping = types.ModuleType("pypath.utils.mapping")
pypath_mapping.__spec__ = ModuleSpec("pypath.utils.mapping", loader=None)
pypath_mapping.__file__ = __file__
pypath_mapping.id_from_label0 = lambda value: value
pypath_mapping.label = lambda value: value
pypath_utils.mapping = pypath_mapping
pypath_module.utils = pypath_utils

sys.modules.setdefault("pypath", pypath_module)
sys.modules.setdefault("pypath.utils", pypath_utils)
sys.modules.setdefault("pypath.utils.mapping", pypath_mapping)

# Minimal pypath_common stub used by the inputs module
pypath_common = types.ModuleType("pypath_common")
pypath_common.__spec__ = ModuleSpec("pypath_common", loader=None)
pypath_common.__file__ = __file__

pypath_common_misc = types.ModuleType("pypath_common.misc")
pypath_common_misc.__spec__ = ModuleSpec("pypath_common.misc", loader=None)
pypath_common_misc.__file__ = __file__
pypath_common_misc.to_list = lambda value: [value] if isinstance(value, str) else list(value or [])
pypath_common.log = lambda *args, **kwargs: None

pypath_common.misc = pypath_common_misc

sys.modules.setdefault("pypath_common", pypath_common)
sys.modules.setdefault("pypath_common.misc", pypath_common_misc)

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# Ensure no pre-installed `neko` package shadows the local sources
for _mod in [name for name in list(sys.modules) if name.startswith("neko")]:
    sys.modules.pop(_mod)

neko_pkg_path = REPO_ROOT / "neko"
neko_init_path = neko_pkg_path / "__init__.py"
neko_spec = importlib.util.spec_from_file_location("neko", neko_init_path)
neko_module = importlib.util.module_from_spec(neko_spec)
neko_spec.loader.exec_module(neko_module)
neko_module.__path__ = [str(neko_pkg_path)]
sys.modules["neko"] = neko_module

from neko.core.network import Network


@pytest.fixture
def sample_resources():
    return pd.DataFrame(
        {
            "source": ["P31749", "Q9Y243", "P46087", "Q9Y243", "P01112"],
            "target": ["P17302", "P31749", "Q15796", "Q9Y6K9", "P84022"],
            "is_directed": [True, True, True, True, True],
            "is_stimulation": [True, False, True, False, True],
            "is_inhibition": [False, True, False, True, False],
            "form_complex": [False, False, False, False, False],
        }
    )


@pytest.fixture
def network(sample_resources):
    net = Network(initial_nodes=["P31749"], resources=sample_resources)
    import inspect

    source_path = inspect.getsourcefile(Network)
    assert source_path and str(REPO_ROOT) in source_path, source_path
    assert hasattr(net, "current_state_id"), dir(net)
    return net


def test_save_undo_redo(network):
    root_state = network.current_state_id
    network.add_node("Q9Y243")
    state_b = network.current_state_id

    states = {entry["id"]: entry["metadata"] for entry in network.list_states()}
    assert states[state_b]["method"] == "add_node"

    network.undo()
    assert network.current_state_id == root_state
    assert "Q9Y243" not in network.nodes["Genesymbol"].tolist()

    network.redo()
    assert network.current_state_id == state_b
    assert "Q9Y243" in network.nodes["Genesymbol"].tolist()


def test_branching_history(network):
    root_state = network.current_state_id

    network.add_node("Q9Y243")
    branch_one = network.current_state_id

    network.checkout(root_state)
    network.add_node("P46087")
    branch_two = network.current_state_id

    history = network.history_graph()
    successors = set(history.successors(root_state))
    assert successors == {branch_one, branch_two}

    # moving between branches
    network.checkout(branch_one)
    network.checkout(branch_two)
    assert set(network.history_graph().successors(root_state)) == {branch_one, branch_two}


def test_compare_states_and_list(network):
    root_state = network.current_state_id
    network.add_node("Q9Y243")
    state_b = network.current_state_id

    network.checkout(root_state)
    network.add_node("P46087")
    state_c = network.current_state_id

    diff = network.compare_states(state_b, state_c)
    assert ("P46087" in diff["added_nodes"]) ^ ("Q9Y243" in diff["added_nodes"])

    listed = network.list_states()
    listed_ids = [entry["id"] for entry in listed]
    assert listed_ids == [root_state, state_b, state_c]


def test_suspend_history(network):
    root_state = network.current_state_id

    with network.suspend_history():
        network.add_node("Q9Y243")

    assert "Q9Y243" in network.nodes["Genesymbol"].tolist()
    assert network.current_state_id == root_state

    network.add_node("P46087")
    state_c = network.current_state_id
    assert state_c != root_state

    network.set_history_tracking(False)
    network.remove_node("P46087")
    assert network.current_state_id == state_c
    network.set_history_tracking(True)

    network.add_node("P46087")
    state_after = network.current_state_id
    assert state_after != state_c


def test_history_digraph(network):
    digraph = network.history_digraph()
    label = digraph.source
    assert str(network.root_state_id) in label



def test_history_html(network):
    html = network.history_html()
    assert "<svg" in html
    assert str(network.root_state_id) in html


def test_max_history_pruning(network):
    network.set_max_history(3)
    created = []
    for node in ["Q9Y243", "P46087", "P01112", "Q15796"]:
        if node not in network.nodes["Genesymbol"].tolist():
            network.add_node(node)
            created.append(network.current_state_id)
    states = network.list_states()
    assert len(states) <= 3
    assert states[0]["id"] == network.root_state_id
    remaining_ids = {entry["id"] for entry in states}
    assert network.current_state_id in remaining_ids
    pruned = [sid for sid in created if sid not in remaining_ids]
    assert pruned
