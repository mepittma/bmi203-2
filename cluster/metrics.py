from .io import read_active_sites, write_clustering, write_mult_clusterings
from .utils import Atom, Residue, ActiveSite
from collections import defaultdict
import numpy as np


def calc_silhouette(clustering):
    """
    This function accepts a clustering assignment (list of lists of clusters) and
    returns a score accounting for how similar objects in a cluster are to each other
    (cohesion) and how distant they are to other clusters (separation).

    Adapted from https://cs.fit.edu/~pkc/classes/ml-internet/silhouette.pdf

    Wow this is nasty! Surely there's a less complex solution?
    """

    # If it's just one giant cluster, return NA
    if len(clustering) == 1:
        return 0 # Assuming we want more than one single cluster

    avg_sil = [] # List that will contain silhouette score for each cluster
    # Works for cluster 0, cluster 1 having some issues
    for i, cluster in enumerate(clustering):

        # Create copy of clustering to safely remove current cluster from separation calculation - looks right
        other_clusters = list(clustering)
        #print("Updated clustering (should be same as input): ", other_clusters)
        other_clusters.pop(i)
        #print("\nUpdated other_clusters (should be all clusters but i): ", other_clusters)

        # Select each point in the cluster and calculate its cohesion with all
        # other points in the cluster
        silhouette = []
        for j, site_a in enumerate(cluster):

            cohesion = []

            # Select other sites in the cluster to calculate cohesion
            other_sites = list(cluster)
            other_sites.pop(j)

            # If the cluster is a singleton, cohesion = 1
            if len(other_sites) == 0:
                cohesion.append(1)

            # Looks like this is working for j = 0
            for site_b in other_sites:
                cohesion.append(compute_similarity(site_a, site_b))
            avg_cohesion = sum(cohesion)/len(cohesion)
            #print("\nAverage cohesion: ", avg_cohesion, "\n")

            # Calculate site_a's average distance to all other clusters - looks right if there's no overlap?
            sep_list = []
            for OC in other_clusters:
                sep_holder = []
                for site_c in OC:
                    sep_holder.append(compute_similarity(site_a, site_c))   #site-by-site similarity
                sep_list.append(sum(sep_holder)/len(sep_holder))            #average similarity from site to clusters
                #print("Calculated separation values - site-by-site: ", sep_holder)
            #print("Calculated separation values - averaged by cluster", sep_list,)

            # Find the maximum similarity from the focal site to another cluster
            separation = min(sep_list)
            #print("Minimum separation from other clusters: ", separation)

            # Calculate the silhouette coefficient for focal site - working
            silhouette.append((avg_cohesion - separation)/max(avg_cohesion, separation))
            #print("silhouette score for active site %r: %r\n\n" %(site_a.name, silhouette[j]))

        # Calculate the average silhouette coefficient for focal cluster
        avg_sil.append(sum(silhouette)/len(silhouette))
        #print("silhouette score for cluster %r: %r\n\n" %(i, avg_sil[i]))

    # Return the average silhouette coefficient for the entire clustering
    return sum(avg_sil)/len(avg_sil)

def seq_scan(active_sites, clusters, t):
    """
    Scans through active site sequences, calculates their similarity, and clusters
    based on similarity to existing sequences in the cluster. If the similarity does
    not meet a previously calculated threshold for any of the existing clusters, then it
    starts a new one.
    """
    for site in active_sites:
        assigned = False
        #print("\n\nAssigning %r to a cluster..." %site.name)

        # try to assign to an existing cluster
        for c in clusters:
            #print("Cluster number: ",c)

            # Find distribution of pairwise distances between already-assigned sites to the query site
            sim_dist = []
            for seen_site in clusters[c]:
                sim_dist.append(compute_similarity(site, seen_site))

            # If median similarity of the query sequence is greater than the threshold
            # cutoff, add it to the cluster; else restart the loop
            if (np.median(sim_dist) >= t):
                clusters[c].append(site)
                assigned = True
                #print("Assigned to cluster %r.\n" %c)
                break

            else:
                #print("Site %r is too different from cluster %r." %(site, c))
                continue

        # If the site wasn't assigned, create a new cluster with the site
        if (assigned == False):
            c = len(clusters) + 1
            #print("Creating new cluster %r.\n" %c)
            clusters[c].append(site)

    return clusters

def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """

    # # # # # # # # # # # # JACCARD SIMILARITY # # # # # # # # # # # #
    similarity = 0.0

    # Only keep the residue names, not their numbers
    a_residues = [str(x).split(' ')[0] for x in site_a.residues]
    b_residues = [str(x).split(' ')[0] for x in site_b.residues]

    # Simple Jaccard similarity
    intersection_cardinality = len(set.intersection(*[set(a_residues), set(b_residues)]))
    union_cardinality = len(set.union(*[set(a_residues), set(b_residues)]))
    similarity = intersection_cardinality/float(union_cardinality)

    return similarity

    # # # # # # # # # # # # CHARGE/SIZE SIMILARITY # # # # # # # # # # # #

def update_sim(clusters, new_cluster, sim_mat):
    """
    Inputs:
    clusters: a list of old clusters to find similarities to the newly-agglomerated cluster
    new_cluster: recently agglomerated cluster of active sites for which we are calculating new distances
    sim_mat: a numpy matrix containing average distances between clusters - constituents of new_cluster already deleted

    Output: an updated matrix in which the most similar clusters from the input
    have been combined into a single cluster
    """

    sim_mat_row = []
    for cluster in clusters:

        # calculate the average distance between this cluster and every other cluster
        avg_sim = []
        for site_a in list(cluster):
            sims = [] # holds the new-cluster similarity scores for site_a
            for site_b in new_cluster:
                #print("Sites: ", site_a, " ", site_b)
                sims.append(compute_similarity(site_a, site_b)) # equals similarity between sites a and b
            avg_sim.append(sum(sims)/len(sims)) # equals average similarity between site a and all sites in new cluster
        sim_mat_row.append(sum(avg_sim)/len(avg_sim)) #equals average similarity between all sites in old cluster and all sites in new cluster

    # Append the new similarities to the similarity matrix
    np_row = np.array([sim_mat_row])

    # Append new column and row
    sim_mat = np.concatenate((sim_mat,np_row), axis=0)
    np_row = np.append(np_row,0)

    sim_mat = np.concatenate((sim_mat,np_row[:, None]), axis=1)
    return sim_mat

def compute_rand(active_sites, h_clust, p_clust):
    """
    Inputs:
        active_sites: a list of ActiveSite instances that were clustered in two different ways
        h_clust: clustering (list of lists of ActiveSites, where each sublist is a distinct cluster)
        p_clust: alternative clustering of same format

    Output: Rand index score indicating the similarity of the two clusterings
    (0 indicates that the clusterings have no pairs in common; 1 indicates that the clusters are functionally the same)
    """

    # Rand index
    # Defined here: https://en.wikipedia.org/wiki/Rand_index
    a = 0   #The number of pairs of elements in active_sites that are in the same cluster in X and in the same cluster in Y
    b = 0   #The number of pairs of elements in active_sites that are in different clusters in both clusterings
    c = 0   #The number of pairs of elements in active_sites that are in the same in X but different in Y
    d = 0   #The number of pairs of elements in active_sites that are in different clusters in X, but same in Y

    # Select all pairs of elements in active site list
    for i, site_a in enumerate(active_sites):

        # Remove focal site from consideration
        other_sites = list(active_sites)
        other_sites.pop(i)

        # For every other site in the dataset, find a,b,c, and d
        for j, site_b in enumerate(other_sites):

            h_status = ""
            p_status = ""

            # Are a and b in the same cluster in h_clust?
            for cluster in h_clust:
                if site_a in cluster:
                    if site_b in cluster:
                        h_status = "SAME"
                    else:
                        h_status = "DIFF"

            # Are they in the same cluster in partitioning?
            for cluster in p_clust:
                if site_a in cluster:
                    if site_b in cluster:
                        p_status = "SAME"
                    else:
                        p_status = "DIFF"

            # Determine which category (a,b,c,d) the pair belongs to
            if h_status == "SAME":
                if p_status == "SAME":
                    a+=1
                elif p_status == "DIFF":
                    c+=1
                else:
                    print("Hmm this is sus")
            elif h_status == "DIFF":
                if p_status == "SAME":
                    d+=1
                elif p_status == "DIFF":
                    b+=1
                else:
                    print("Hmm this is sus")

    # Calculate the Rand index
    R = (a+b)/(a+b+c+d)
    print("Rand index (similarity of my two clusterings): ", R)

    return R
