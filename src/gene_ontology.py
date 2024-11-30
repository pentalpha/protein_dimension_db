import json
import networkx as nx
import obonet

from util import go_not_use_path, go_basic_path

def gos_not_to_use():
    notuse_json = json.load(open(go_not_use_path, 'r'))
    nodes = notuse_json['graphs'][0]['nodes']
    ids = [nodes[i]['id']
        for i in range(len(nodes))]
    goids = ['GO:'+x.split('_')[-1] for x in ids if 'GO_' in x]
    return set(goids)

def load_go_graph():
    graph = obonet.read_obo(go_basic_path)
    return graph

def expand_go_set(goid: str, go_graph: nx.MultiDiGraph, goes_to_not_use: set):
    all_gos = set()

    if goid in go_graph:
        parents = sorted(nx.descendants(go_graph, goid))
        for parent in parents:
            if not parent in goes_to_not_use:
                all_gos.add(parent)
    
        all_gos.add(goid)
    
    return sorted(all_gos)