import osmnx as ox
import pandas as pd
import networkx as nx
import math


places = ["Winchester, MA", "Arlington, MA", "Woburn, MA", "Medford, MA", "Burlington, MA"]
dist_factor_vals = []

for p in places:

    G = ox.graph_from_place(query=p,
                            network_type="drive",
                            simplify=True,
                            retain_all=False)
    G = ox.project_graph(G, to_crs="EPSG:26986")

    labels = []
    for n in G.nodes(data=True):
        labels.append(n[0])

    #print(labels)
    dist = nx.all_pairs_dijkstra_path_length(G, weight="length")
    ratios = []

    for src, targets in dist:
        #print(G.nodes[src])
        src_x = G.nodes[src]['x']
        src_y = G.nodes[src]['y']
        for dst, d in targets.items():
            dst_x = G.nodes[dst]['x']
            dst_y = G.nodes[dst]['y']

            euc = math.sqrt((src_x-dst_x)**2 + (src_y-dst_y)**2)

            if (d == 0) | (euc == 0):
                continue
            ratios.append(d / euc)
            #print(d)
            pass

    dist_factor = sum(ratios) / len(ratios)
    dist_factor_vals.append(dist_factor)
    print(f"Place: {p}, Dist_Factor: {dist_factor}")

print(f"Average Dist Factor: {sum(dist_factor_vals)/len(dist_factor_vals)}")
