"""Test production SIF functions without importing provider/network clients.

Compile the real methods and pure consolidation helpers from their source AST;
only identifier lookup, node storage and object synchronization are mocked.
"""
import ast
from contextlib import ExitStack
from pathlib import Path
from types import SimpleNamespace
import os
import socket
import tempfile
import unittest
from unittest import mock

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def functions_from_source(path, names, namespace, class_name=None):
    tree = ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
    body = tree.body
    if class_name:
        body = next(node.body for node in body
                    if isinstance(node, ast.ClassDef) and node.name == class_name)
    nodes = [node for node in body
             if isinstance(node, ast.FunctionDef) and node.name in names]
    assert {node.name for node in nodes} == set(names)
    module = ast.Module(body=nodes, type_ignores=[])
    exec(compile(module, str(path), 'exec'), namespace)
    return namespace


FUNCTIONS = functions_from_source(ROOT / 'neko/core/tools.py',
    ['_is_missing', '_join_unique_values', '_normalize_effect', 'consolidate_edges'],
    {'pd': pd})
CONSOLIDATE = FUNCTIONS['consolidate_edges']
PARSER = functions_from_source(ROOT / 'neko/core/network.py',
    ['_load_network_from_sif'],
    {'pd': pd, 'consolidate_edges': CONSOLIDATE,
     'mapping_node_identifier': lambda node: [None, node, node],
     'check_gene_list_format': lambda nodes: True},
    class_name='Network')['_load_network_from_sif']
EXPORT_SIF = functions_from_source(ROOT / 'neko/_outputs/exports.py',
    ['export_sif'], {'os': os}, class_name='Exports')['export_sif']


class ParserHarness:
    def __init__(self):
        self.edges = pd.DataFrame(columns=[
            'source', 'target', 'Type', 'Effect', 'References',
        ])
        self.nodes = set()

    def add_node(self, node, from_sif=False):
        assert from_sif is True
        self.nodes.add(node)

    def sync_edges_from_df(self):
        pass


def read_sif(tmp_path, content, name='input.sif'):
    source = tmp_path / name
    source.write_text(content, encoding='utf-8')
    network = ParserHarness()
    PARSER(network, str(source))
    return network


class SIFReferenceTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.tmp_path = Path(temporary.name)
        stack = ExitStack()
        self.addCleanup(stack.close)

        def no_network(*args, **kwargs):
            self.fail('SIF parser test attempted network access')

        for obj, name in ((socket.socket, 'connect'),
                          (socket.socket, 'connect_ex'),
                          (socket, 'create_connection'),
                          (socket, 'getaddrinfo')):
            stack.enter_context(mock.patch.object(obj, name, no_network))

    def test_export_import_roundtrip_preserves_each_reference_string(self):
        tmp_path = self.tmp_path
        content = (
            '# Reference PMID: 101; 202; 303\nnode_A\tstimulation\tnode_B\n'
            '# Reference PMID: 404\nnode_B\tinhibition\tnode_C\n'
            '# Reference PMID: source-record:alpha\n'
            'node_C\tstimulation\tnode_A\n'
        )
        original = read_sif(tmp_path, content)
        exported = tmp_path / 'exported.sif'
        EXPORT_SIF(SimpleNamespace(interactions=CONSOLIDATE(original.edges)), str(exported))
        assert exported.read_text() == content
        imported = read_sif(tmp_path, exported.read_text(), 'roundtrip.sif')
        pd.testing.assert_frame_equal(original.edges, imported.edges)
        assert len(imported.nodes) == 3
        assert imported.edges['References'].tolist() == [
            '101; 202; 303', '404', 'source-record:alpha',
        ]

    def test_reference_applies_once_and_ordinary_sif_keeps_fallback(self):
        tmp_path = self.tmp_path
        network = read_sif(tmp_path,
            '# Reference PMID: 101\nnode_A stimulation node_B\n'
            'node_B inhibition node_C\n'
            '# Reference PMID: 303\nnode_C stimulation node_D\n'
            'node_D inhibition node_A\n')
        assert network.edges['References'].tolist() == [
            '101', 'SIF file', '303', 'SIF file',
        ]

    def test_comments_blank_lines_and_crlf_before_edge(self):
        tmp_path = self.tmp_path
        network = read_sif(tmp_path,
            '# Reference PMID: 101; 202\r\n\r\n'
            '# explanatory comment\r\n \t\r\n'
            'node_A stimulation node_B\r\nnode_B inhibition node_C\r\n')
        assert network.edges['References'].tolist() == ['101; 202', 'SIF file']

    def test_malformed_interaction_consumes_its_annotation(self):
        tmp_path = self.tmp_path
        malformed = 'node_A stimulation'
        network = read_sif(tmp_path,
            f'# Reference PMID: 101\n{malformed}\n'
            'node_B stimulation node_C\n'
            '# Reference PMID: 303\nnode_C inhibition node_D\n')
        assert network.edges['References'].tolist() == ['SIF file', '303']
        assert len(network.nodes) == 3

    def test_empty_reference_comment_cannot_reuse_previous_annotation(self):
        tmp_path = self.tmp_path
        network = read_sif(tmp_path,
            '# Reference PMID: 101\n# Reference PMID: \n'
            'node_A stimulation node_B\n')
        assert network.edges.loc[0, 'References'] == 'SIF file'

    def test_annotated_duplicate_signs_merge_existing_evidence_semantics(self):
        tmp_path = self.tmp_path
        network = read_sif(tmp_path,
            '# Reference PMID: 101; 202\nnode_A stimulation node_B\n'
            '# Reference PMID: 202; 303\nnode_A inhibition node_B\n'
            '# Reference PMID: 101\nnode_A stimulation node_B\n')
        assert len(network.edges) == 1
        assert network.edges.loc[0, 'Effect'] == 'bimodal'
        assert network.edges.loc[0, 'References'] == '101; 202; 303'

    def test_unannotated_opposite_signs_keep_existing_fallback(self):
        tmp_path = self.tmp_path
        network = read_sif(tmp_path,
            'node_A stimulation node_B\nnode_A inhibition node_B\n')
        assert len(network.edges) == 1
        assert network.edges.loc[0, 'Effect'] == 'bimodal'
        assert network.edges.loc[0, 'References'] == 'SIF file'

    def test_optional_type_and_unknown_effect_are_unchanged(self):
        tmp_path = self.tmp_path
        network = read_sif(tmp_path,
            '# Reference PMID: 101\nnode_A unknown node_B binding\n'
            '# Reference PMID: 202\nnode_B form_complex node_C\n')
        assert network.edges['Effect'].tolist() == ['undefined', 'form complex']
        assert network.edges.loc[0, 'Type'] == 'binding'
        assert pd.isna(network.edges.loc[1, 'Type'])
        assert network.edges['References'].tolist() == ['101', '202']

    def test_dangling_reference_is_not_carried_into_another_import(self):
        tmp_path = self.tmp_path
        first = read_sif(tmp_path,
            'node_A stimulation node_B\n# Reference PMID: 999\n')
        second = read_sif(tmp_path, 'node_C inhibition node_D\n', 'next.sif')
        assert first.edges.loc[0, 'References'] == 'SIF file'
        assert second.edges.loc[0, 'References'] == 'SIF file'


if __name__ == "__main__":
    unittest.main()
