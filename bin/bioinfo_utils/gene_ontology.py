import json
import networkx as nx
import obonet
from typing import List

# from bioinfo_utils.util import go_not_use_path, go_basic_path

GO_ROOTS = {
    "MF": "GO:0003674",
    "CC": "GO:0005575",
    "BP": "GO:0008150",
}

NAMESPACES = {
    "biological_process": "BP",
    "molecular_function": "MF",
    "cellular_component": "CC",
}


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


def list_ontology_members(go_obo_path):
    """
    [...]
    [Term]
    id: GO:0000001
    name: mitochondrion inheritance
    namespace: biological_process
    def: [...]
    """
    go_lists = {namespace: set() for namespace in NAMESPACES.keys()}
    go_alt_ids = {}

    last_goid = None
    for rawline in open(go_obo_path, "r"):
        line = rawline.strip()
        if line.startswith("id: ") and "GO:" in line:
            last_goid = line.replace("id: ", "")
        elif line.startswith("namespace: ") and last_goid:
            namespace = line.replace("namespace: ", "")
            if namespace != "external":
                go_lists[namespace].add(last_goid)
        elif line.startswith("alt_id: "):
            alt_id = line.replace("alt_id: ", "")
            go_alt_ids[alt_id] = last_goid

    for namespace, ont in NAMESPACES.items():
        go_lists[ont] = go_lists[namespace]
        del go_lists[namespace]

    return go_lists, go_alt_ids


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
