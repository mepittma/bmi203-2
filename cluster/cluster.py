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
    for cluster in clustering:
        size_dist.append(len(cluster))

    smalls = np.percentile(dist,10)

    # Extract sites from smallest clusters
    small_sites = []
    for small in smalls:
        small_sites.append(value(small))

    # Remove smallest clusters from clustering
    big_clusters = list(set(clustering)^set(smalls))

    # Rerun through sequence scan, try to reassign these active sites
    final_clustering = seq_scan(small_sites, big_clusters, t)

    return final_clustering


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Average Linkage clustering, an agglomerative method
    #https://home.deib.polimi.it/matteucc/Clustering/tutorial_html/hierarchical.html

    # Get list of names of active sites


    # 1. Calculate pairwise similarity for active sites
    mat = []
    for site_i in active_sites:
        j = []
        for site_j in active_sites:
            sim = compute_similarity(site_i, site_j)
            j.append(sim)
        mat.append(j)

    # 2. Find the most-similar pair in the matrix


    return clustering
