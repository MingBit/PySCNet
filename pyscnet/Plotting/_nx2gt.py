import graph_tool as gt

def get_prop_type(value, key=None):
    if isinstance(key, str):
        key = key.encode('utf-8', errors='replace')

    if isinstance(value, bool):
        tname = 'bool'
    elif isinstance(value, int):
        tname = 'float'
        value = float(value)
    elif isinstance(value, float):
        tname = 'float'
    elif isinstance(value, str):
        tname = 'string'
        value = value.encode('utf-8', errors='replace')
    elif isinstance(value, dict):
        tname = 'object'
    else:
        tname = 'string'
        value = str(value)

    try:
        key = key.decode('utf-8')
    except AttributeError:
        pass

    return tname, value, key

def nx2gt(nxG):
    gtG = gt.Graph(directed=nxG.is_directed())

    for key, value in list(nxG.graph.items()):
        tname, value, key = get_prop_type(value, key)
        prop = gtG.new_graph_property(tname)
        gtG.graph_properties[key] = prop
        gtG.graph_properties[key] = value

    nprops = set()
    for node, data in nxG.nodes(data=True):
        for key, val in list(data.items()):
            if key in nprops:
                continue
            tname, _, key = get_prop_type(val, key)
            prop = gtG.new_vertex_property(tname)
            gtG.vertex_properties[key] = prop
            nprops.add(key)

    gtG.vertex_properties['id'] = gtG.new_vertex_property('string')

    eprops = set()
    for src, dst, data in nxG.edges(data=True):
        for key, val in list(data.items()):
            if key in eprops:
                continue
            tname, _, key = get_prop_type(val, key)
            prop = gtG.new_edge_property(tname)
            gtG.edge_properties[key] = prop
            eprops.add(key)

    vertices = {}
    for node, data in nxG.nodes(data=True):
        v = gtG.add_vertex()
        vertices[node] = v
        data['id'] = str(node)
        for key, value in list(data.items()):
            gtG.vp[key][v] = value

    for src, dst, data in nxG.edges(data=True):
        e = gtG.add_edge(vertices[src], vertices[dst])
        for key, value in list(data.items()):
            gtG.ep[key][e] = value

    return gtG
