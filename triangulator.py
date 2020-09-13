import pandas as pd
import numpy as np
import networkx as nx

# Description: take a complete distance matrix that may not satisfy the triangle inequality
#   and replace all arc lengths by their shortest paths (which ensures triangle inequality).
# Note: it may not be appropriate to apply this code separately to a time and distance matrix,
#   as the shortest paths that are found may differ. A better alternative would be to modify this code,
#   find a shortest path in time or distance and then update the other accordingly (possibly violating the triangle inequality).

# Read matrix from file
# Matrix has row and column headers, which are equal and in the same order
matrix = pd.read_csv("matrix.csv")
matrix['stop_id'] = matrix['stop_id'].astype(int).astype(str)
matrix.set_index('stop_id', inplace=True)
matrix.columns = matrix.columns.astype(float).astype(int).astype(str)
assert all(matrix.index == matrix.columns)

# Shrink matrix for testing purposes
matrix = matrix.loc[['905879','905309','900432'], ['905879','905309','900432']] # TODO: REMOVE LINE

# Prepare list of edges and weights
edges = matrix.apply(lambda x: [(x.name, to, x[to]) for to in x.index], axis=1)
edges = edges.sum()

# Build graph
graph = nx.DiGraph()
graph.add_nodes_from(matrix.index)
graph.add_weighted_edges_from(edges)

# Solve all-pair shortest path
lengths = nx.floyd_warshall(graph)

# Convert lengths to a DataFrame matrix
length_matrix = pd.DataFrame(index=matrix.index, columns=matrix.columns)
for origin in lengths.keys():
    length_matrix.loc[origin,:] = [lengths[origin][destination] for destination in matrix.columns]

# matrix: the original matrix as a DataFrame
# length_matrix: matrix with shortest distances as a DataFrame
# If matrix - length_matrix > 0 for any cell, then matrix did not satify the triangle inequality

# No export is implemented