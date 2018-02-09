import sys
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import cluster_by_partitioning#, cluster_hierarchically
from .metrics import calc_silhouette, seq_scan, compute_similarity

# Some quick stuff to make sure the program is called correctly
#if len(sys.argv) < 4:
#    print("Usage: python -m hw2skeleton [-P| -H] <pdb directory> <output file>")
#    sys.exit(0)

#active_sites = read_active_sites(sys.argv[2])

# Choose clustering algorithm
#if sys.argv[1][0:2] == '-P':
#    print("Clustering using Partitioning method")
#    clustering = cluster_by_partitioning(active_sites)
#    write_clustering(sys.argv[3], clustering)

#if sys.argv[1][0:2] == '-H':
#    print("Clustering using hierarchical method")
#    clusterings = cluster_hierarchically(active_sites)
#    write_mult_clusterings(sys.argv[3], clusterings)




# # # # # # # # # # # # test hierarchical clustering # # # # # # # # # # # # # #
from collections import defaultdict
import random
import numpy as np

# Read in all the active sites to cluster
path = "/Users/student/Documents/BMI206/bmi203-2/data/"
active_sites = read_active_sites(path)

# Assign each site to its own cluster
clusters = []
sim_mat = np.zeros((len(active_sites),len(active_sites)))
print("Dimensions of matrix: ", sim_mat.shape)

for i, site_a in enumerate(active_sites):
    clusters.append([site_a])

    # Calculate the similarity between each "cluster"
    for j, site_b in enumerate(active_sites):
        sim_mat[i,j] = compute_similarity(site_a, site_b)

# Update diagonal to 0s
np.fill_diagonal(sim_mat, 0)

# Print preliminary clustering and similarity matrix - looks right I think
print("\nInitial clusters (should be every ActiveSite for himself): ", clusters)
print("\nInitial similarity matrix: ", sim_mat)

# While there are more than one clusters, iteratively add the most similar clusters
clusterings = []
scores = []
while len(clusters) > 1:

    # Find the maximum average similarity between all clusters
    print("\nHighest similarity value: ", np.max(sim_mat))

    # For each row in the similarity matrix, select the index of max similarity
    # Note: if more than two clusters share the same max similarity, only the first is chosen
    # (this algorithm is subject to produce different results based on initial conditions)

    # Decay this into the x and y axes, which will be the indexes to combine in clusters - looks good
    print("Focal index: ", np.argmax(sim_mat))
    row = np.argmax(sim_mat) // len(clusters)
    column = np.argmax(sim_mat) % len(clusters)
    print("Clusters to combine: %d, %d" %(row, column))

    # Update the clustering such that the row-th and column-th cluster are combined
    new_cluster = []
    new_cluster.append(clusters[row])
    new_cluster.append(clusters[column])
    new_cluster = [item for sublist in new_cluster for item in sublist]
    for i in sorted([row, column], reverse=True):
        del clusters[i]
    print("New cluster: ", new_cluster)

    # Update the similarity matrix
    # remove rows & columns of index row, column
    for i in sorted([row, column], reverse=True):
        sim_mat = np.delete(sim_mat, i, 0)
        sim_mat = np.delete(sim_mat, i, 1)

    # append new row & column to end of new cluster
    sim_mat_row = []
    for cluster in clusters:

        print("current focal comparison cluster: ", cluster)

        # calculate the average distance between this cluster and every other cluster
        avg_sim = []
        for site_a in list(cluster):
            sims = [] # holds the new-cluster similarity scores for site_a
            for site_b in new_cluster:
                print("Sites: ", site_a, " ", site_b)
                sims.append(compute_similarity(site_a, site_b)) # equals similarity between sites a and b
            avg_sim.append(sum(sims)/len(sims)) # equals average similarity between site a and all sites in new cluster
        sim_mat_row.append(sum(avg_sim)/len(avg_sim)) #equals average similarity between all sites in old cluster and all sites in new cluster

    # Append the new similarities to the similarity matrix
    np_row = np.array([sim_mat_row])

    # Append new column and row
    print("Shape of vector: ", np_row.shape)
    print("Shape of matrix: ", sim_mat.shape)

    sim_mat = np.concatenate((sim_mat,np_row), axis=0)
    np_row = np.append(np_row,0)
    print("New shape of matrix: ", sim_mat.shape)
    print("New shape of vector: ", np_row[:, None].shape)
    sim_mat = np.concatenate((sim_mat,np_row[:, None]), axis=1)

    print("New shape of matrix: ", sim_mat.shape)

    # Append the new clustering to the clusterings list & score
    clusters.append(new_cluster)
    print("\nClusters being appended into clusterings: ", clusters)
    clusterings.append(clusters)

    print("\nInput into calc_silhouette: ", clusters)
    scores.append(calc_silhouette(clusters))
