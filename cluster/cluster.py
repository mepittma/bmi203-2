from .utils import Atom, Residue, ActiveSite
from .metrics import calc_silhouette, seq_scan, compute_similarity, update_sim
from collections import defaultdict
import random
import numpy as np

def cluster_by_partitioning(active_sites):#, p):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """

    p = 15.0 # Predetermined "ideal percentile cutoff" for this dataset

    # Covering algorithm - modified from Threshold Bootstrap, Prosperi et al. 2010
    # (no bootstrapping of aligned residues)

    # 1. Initialization - shuffle the order of active sites
    random.shuffle(active_sites)

    # 2. Threshold assessment - calculate a prior distribution of distances, selected
    # randomly from the pairwise distances between active sites (n-squared)
    dist = []
    for i in range(0,len(active_sites) ** 2):

        # select two random active sites
        sites = random.sample(active_sites, 2)

        # calculate similarity
        sim = compute_similarity(sites[0], sites[1])

        # add to distribution
        dist.append(sim)

    # calculate threshold distance value (lower 10th percentile of similarity distribution)
    t = np.percentile(dist,p)

    # 3. Cluster assignment - seq through the active site sequences and cluster
    # based on similarity and threshold cutoff
    clustering = seq_scan(active_sites, defaultdict(list), t)

    # 4. Refinement: attempt to recluster the 10% smallest clusterings

    # Find the smallest 10% of clusters
    size_dist = []
    for cluster, sites in clustering.items():
        size_dist.append(len(sites))

    small_t = np.percentile(dist,10)
    small_list = list([ cluster for cluster, sites in clustering.items() if len(sites) <= small_t])

    # Extract sites from smallest clusters
    small_sites = []
    for small_c in small_list: #grabs cluster
        for active_site in clustering[small_c]: #grabs each active site in the cluster
             small_sites.append(active_site)

    # Remove smallest clusters from clustering - create new dict with only big_clusters
    big_dict = defaultdict(list,{cluster: sites for cluster, sites in clustering.items()
                 if cluster not in small_list})

    # Rerun through sequence scan, try to reassign these active sites
    final_clustering = seq_scan(small_sites, big_dict, t)
    # Recreate as a list of lists
    final_clust_list = list([sites for cluster, sites in final_clustering.items()])

    return final_clust_list


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Average Linkage clustering, an agglomerative method
    #https://home.deib.polimi.it/matteucc/Clustering/tutorial_html/hierarchical.html

    # Create initial similarity matrix, list of clusters (currently contain one ActSite each)
    clusters = []
    sim_mat = np.zeros((len(active_sites),len(active_sites)))

    # Loop through active sites, individually appending to list of clusters
    for i, site_a in enumerate(active_sites):
        clusters.append([site_a])

        # Calculate the similarity between each "cluster"
        for j, site_b in enumerate(active_sites):
            sim_mat[i,j] = compute_similarity(site_a, site_b)

    # Update diagonal to 0s
    np.fill_diagonal(sim_mat, 0)

    # Arrays to store clusterings and their respective scores
    clusterings = []
    scores = []

    # While the number of clusters is greater than 1, iteratively find the most
    # similar clusters and merge them. Update the similarity matrix.
    while len(clusters) > 1:

        # For each row in the similarity matrix, select the index of max similarity
        # Note: if more than two clusters share the same max similarity, only the first is chosen
        # (thus, this algorithm is subject to produce different results based on initial conditions)

        # Decay max index into the x and y axes, which will be the indexes to combine in clusters
        row = np.argmax(sim_mat) // len(clusters)
        column = np.argmax(sim_mat) % len(clusters)

        # Update the clustering such that the row-th and column-th cluster are combined
        new_cluster = []
        new_cluster.append(clusters[row])
        new_cluster.append(clusters[column])
        new_cluster = [item for sublist in new_cluster for item in sublist]
        for i in sorted([row, column], reverse=True):
            del clusters[i]

        # Update the similarity matrix
        # remove rows & columns of index row, column - make sure higher numbers deleted first
        for i in sorted([row, column], reverse=True):
            sim_mat = np.delete(sim_mat, i, 0)
            sim_mat = np.delete(sim_mat, i, 1)

        sim_mat = update_sim(clusters, new_cluster, sim_mat)

        # Append the new clustering to the clusterings list & score
        clusters.append(new_cluster)
        clusterings.append(clusters[:]) # why does python do this to me

        scores.append(calc_silhouette(clusters))

    # Only return cluster with highest silhouette score - choosing most "divisive" (first seen with max value)
    scores.index(max(scores))

    return clusterings#[scores.index(max(scores))] #uncomment to observe highest-scoring clustering (most numerous)

    # Comment 142 and Uncomment below to observe biggest highest-scoring clusters
    #rev_clust = list(reversed(clusterings))
    #rev_score = list(reversed(scores))
    #return rev_clust[rev_score.index(max(rev_score))]
