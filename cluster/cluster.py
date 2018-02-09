from .utils import Atom, Residue, ActiveSite
from .metrics import calc_silhouette, seq_scan, compute_similarity
from collections import defaultdict
import random
import numpy as np

def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
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
    t = np.percentile(dist,10)

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

    # Assign each site to its own cluster
    clusters = []
    sim_mat = np.zeros((len(clusters),len(clusters)))

    for site_a, i in enumerate(active_sites):
        clusters[i] = list(site_a)

        # Calculate the similarity between each "cluster"
        for site_b, j in enumerate(active_sites):
            sim_mat[i,j] = compute_similarity(site_a, site_b)

    # Print preliminary clustering and similarity matrix
    print("Initial clusters (should be every ActiveSite for himself): ", clusters)
    print("")


    # While there are more than one clusters, iteratively add the most similar clusters
    while len(clusters) > 1:

        # Find the minimum average distance between all clusters
        np.matrix()



    return clustering
