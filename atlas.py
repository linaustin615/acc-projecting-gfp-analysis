"""Read aligned Allen IDs without replacing or re-registering the source masks."""
import json
import re
import numpy as np
from config import ROOT


def ontology():
    tree = json.loads((ROOT / 'atlas/allen_structure_graph.json').read_text())
    nodes, paths = {}, {}
    def visit(node, ancestors):
        nodes[node['id']] = node
        paths[node['id']] = ancestors + [node['id']]
        for child in node['children']:
            visit(child, paths[node['id']])
    visit(tree['msg'][0], [])
    return nodes, paths


def display_group(sid, nodes, paths):
    if sid == 0:
        return 0
    node = nodes[sid]
    parent = nodes.get(node['parent_structure_id'])
    is_layer = bool(re.search(r'\blayers?\b|,\s*[1-6](?:/[23]|[ab])?$', node['name'], re.I))
    if parent and is_layer and 688 in paths[sid] and node['acronym'].startswith(parent['acronym']):
        suffix = node['acronym'][len(parent['acronym']):]
        if suffix in ['1', '2', '2/3', '3', '4', '5', '6', '6a', '6b']:
            return parent['id']
    return sid


def lookup(points, labels):
    points = np.asarray(points, dtype=float).reshape(-1, 2)
    if not np.isfinite(points).all() or not np.all((points >= 0) & (points < labels.shape)):
        raise ValueError('Detection coordinates fall outside the aligned label image')
    pixels = np.minimum(np.rint(points).astype(int), np.array(labels.shape)-1)
    return labels[pixels[:, 0], pixels[:, 1]].astype(np.int64)
