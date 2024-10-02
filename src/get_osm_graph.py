import argparse
import osmnx as ox
import pandas as pd

from shapely.geometry import LineString, Point
from shapely.ops import substring


def refine_graph(G, len_limit):
    repeat = True
    idx = 1
    while repeat == True:
        #print("Starting Refinement Iter")
        # Find Edges to reduce
        split_edges = []
        for (src, dst, data) in G.edges(data=True):
            
            if data["length"] > len_limit:
                split_edges.append((src,dst))


        #Reduce Edges
        repeat, idx = shorten_edges(G, split_edges, idx, len_limit)
        

def shorten_edges(G, split_edges, node_idx, len_limit):
    repeat = False
    mid = node_idx
    for (src, dst) in split_edges:
        
        # Get Edge Geometry
        data = G[src][dst][0]
        if "geometry" in data:
            geo = data["geometry"]
        else:
            geo = LineString([Point(G.nodes[src]["x"], G.nodes[src]["y"]), Point(G.nodes[dst]["x"], G.nodes[dst]["y"])])


        # Compute new line segments
        s1_geo = substring(geo, 0, 0.5, normalized=True)
        s2_geo = substring(geo, 0.5, 1, normalized=True)
        
        

        s2_proj = ox.projection.project_geometry(s2_geo, crs="EPSG:26986", to_latlong=True)

        # Add new node
        G.add_node(mid, x=s2_geo.coords[0][0], y=s2_geo.coords[0][1], street_count=2, lon=s2_proj[0].coords[0][0], lat=s2_proj[0].coords[0][1])
        
        # Add new edges
        G.add_edge(src, mid, **{**data, 'length': s1_geo.length, 'geometry': s1_geo})
        G.add_edge(mid, dst, **{**data, 'length': s2_geo.length, 'geometry': s2_geo})

        # Check if length is over limit
        if s1_geo.length > len_limit:
            repeat = True
        

        # Remove orginal edge
        G.remove_edge(src, dst)

        mid = mid + 1
    return repeat, mid
 
# Initialize parser
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("--file", 
                    help = "Path to output file",
                    type = str,
                    default= "OSMGraph")
parser.add_argument("--travel_speed",
                    help= "Travel speed in [m/s] for edge weight computation",
                    type = float,
                    default=1.25)
parser.add_argument("--place_str",
                    help= "String to pass into the OSMnx functions",
                    type= str,
                    default= "Winchester, MA")
parser.add_argument("--filters",
                    help= "Custom filters to pass into OSMnx",
                    type= str)
parser.add_argument("--CRS",
                    help = "CRS to project graph into. Example: EPSG:26986",
                    type = str,
                    default= "EPSG:4326")
parser.add_argument("--resolution",
                    help = "Maximum edge length in [m]",
                    type = int,
                    default = -1)
 
# Read arguments from command line
args = parser.parse_args()

walk_filters = (
            f'["highway"]["area"!~"yes"]["access"!~"private"]'
            f'["highway"!~"abandoned|bus_guideway|construction|cycleway|motor|no|planned|platform|proposed|raceway|razed"]'
            f'["foot"!~"no"]["service"!~"private"]'
            f'["highway"!~"service"]'
        )
## Load Graph and Reproject
G = ox.graph_from_place(query=args.place_str,
                        simplify=False,
                        custom_filter=args.filters,
                        retain_all=False)
G = ox.project_graph(G, to_crs=args.CRS)


## Simplify Graph
G = ox.simplify_graph(G,
        remove_rings=True,
        track_merged=False
        )

## Improve Resoluiton
if args.resolution != -1:
    refine_graph(G, args.resolution)


## Add Edge weights
for edge in G.edges(data=True):
    edge[2]["weight"] = edge[2]["length"] / args.travel_speed 
    edge[2]["mode"] = "walk"

## Save Output

# Create dataframe of the node data
node_data = []
i = 1
for node in G.nodes(data=True):
    rowDict = {}
    rowDict["index"] = i
    rowDict["nx_idx"] = node[0]
    rowDict["x"] = node[1]["x"]
    rowDict["y"] = node[1]["y"]
    rowDict["lat"] = node[1]["lat"]
    rowDict["lon"] = node[1]["lon"]
    node_data.append(rowDict)
    i = i + 1
node_df = pd.DataFrame(node_data)

# Create dataframe of edge data
edge_data = []
i = 1
for edge in G.edges(data=True):
    try:
        w = edge[2]["weight"]
    except:
        w = -1

    rowDict = {}
    rowDict["index"] = i
    rowDict["src"] = edge[0]
    rowDict["dst"] = edge[1]
    rowDict["weight"] = w
    rowDict["length"] = edge[2]["length"]
    edge_data.append(rowDict)
    i = i + 1
edge_df = pd.DataFrame(edge_data)

# Replace the OSMnx labels with the indexs
edge_df["src"] = edge_df["src"].replace(node_df["nx_idx"].to_list(), node_df["index"].to_list())
edge_df["dst"] = edge_df["dst"].replace(node_df["nx_idx"].to_list(), node_df["index"].to_list())

# Remove OSMnx indexs
node_df.drop(labels="nx_idx", axis=1, inplace=True)

# Find duplicate edges and remove - Better to check if reverse edge exists when building the combined graph?
n = edge_df.shape[0]
rmv_edges = []
for i in  range(n):
    src = edge_df.src[i]
    dst = edge_df.dst[i]

    for j in range(i+1,n,1):
        src2 = edge_df.src[j]
        dst2 = edge_df.dst[j]

        if src2 == dst and dst2 == src:
            rmv_edges.append(j)

edge_df = edge_df.drop(rmv_edges)
edge_df['index'] = [*range(1,edge_df.shape[0]+1,1)]

# Write to CSV files
node_df.to_csv(args.file + ".nodes.csv", index=False)
edge_df.to_csv(args.file + ".edges.csv", index=False)
