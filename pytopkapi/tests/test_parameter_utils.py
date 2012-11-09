"""Tests for the parameter_utils sub-package.

"""
import numpy as np
import networkx as nx

from pytopkapi.parameter_utils.create_file import _make_strahler_dicts
from pytopkapi.parameter_utils.create_file import strahler_stream_order

def test_strahler_simple():
    nodes = np.arange(5)
    edges = [(4,3), (3,1), (2,1), (1,0)]

    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    stream_orders = {}

    nodes_per_arc, arcs_per_node = _make_strahler_dicts(G)
    strahler_stream_order(0, 1, nodes_per_arc, arcs_per_node, stream_orders)

    solution = {0 : 2, 1 : 1, 2 : 1, 3 : 1}
    assert(stream_orders == solution)

def test_strahler_complex():
    nodes = np.arange(14)
    edges = [(13,11), (12,11), (11,9), (10,9),
             (9,1), (1,0), (8,6), (7,6), (6,4), (5,4), (4,2), (3,2), (2,1)]

    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    stream_orders = {}

    nodes_per_arc, arcs_per_node = _make_strahler_dicts(G)
    strahler_stream_order(0, 1, nodes_per_arc, arcs_per_node, stream_orders)

    solution = {0: 3, 1: 2, 2: 1, 3: 2, 4: 1,
                5: 2, 6: 1, 7: 1, 8: 2, 9: 1, 10: 2, 11: 1, 12: 1}
    assert(stream_orders == solution)
