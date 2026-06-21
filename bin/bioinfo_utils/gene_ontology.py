import json
import networkx as nx
import obonet

# from bioinfo_utils.util import go_not_use_path, go_basic_path


def gos_not_to_use(go_not_use_path):
    notuse_json = json.load(open(go_not_use_path, "r"))
    nodes = notuse_json["graphs"][0]["nodes"]
    ids = [nodes[i]["id"] for i in range(len(nodes))]
    goids = ["GO:" + x.split("_")[-1] for x in ids if "GO_" in x]
    return set(goids)


def load_go_graph(go_basic_path):
    graph = obonet.read_obo(go_basic_path)
    return graph


def load_go_graph_strict(go_basic_path):
    # Load the raw multi-edge graph
    raw_graph = obonet.read_obo(go_basic_path)

    # Create a new directed graph for strict True Path Rule propagation
    strict_graph = nx.DiGraph()

    # We only want to propagate over 'is_a' and 'part_of'
    valid_relations = {"is_a", "part_of"}

    # Add all nodes first to ensure we don't lose any disconnected terms
    strict_graph.add_nodes_from(raw_graph.nodes(data=True))

    # Iterate through edges and only keep the valid ones
    for u, v, key, data in raw_graph.edges(keys=True, data=True):
        if key in valid_relations:
            strict_graph.add_edge(u, v)

    return strict_graph


def expand_go_set(goid: str, go_graph: nx.MultiDiGraph, goes_to_not_use: set):
    all_gos = set()

    if goid in go_graph:
        parents = sorted(nx.descendants(go_graph, goid))
        for parent in parents:
            if not parent in goes_to_not_use:
                all_gos.add(parent)

        all_gos.add(goid)

    return sorted(all_gos)


def expand_go_set_down(goid: str, go_graph: nx.MultiDiGraph, goes_to_not_use: set):
    all_gos = set()

    if goid in go_graph:
        parents = sorted(nx.ancestors(go_graph, goid))
        for parent in parents:
            if not parent in goes_to_not_use:
                all_gos.add(parent)

        all_gos.add(goid)

    return sorted(all_gos)


from typing import List


def expand_go_list(
    go_list: List[str],
    go_graph: nx.MultiDiGraph,
    goes_to_not_use: set,
    is_negative: bool = False,
):
    expanded_gos = set()
    for goid in go_list:
        if is_negative:
            expanded_gos.update(expand_go_set_down(goid, go_graph, goes_to_not_use))
        else:
            expanded_gos.update(expand_go_set(goid, go_graph, goes_to_not_use))
    return sorted(expanded_gos)
