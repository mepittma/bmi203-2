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

    print("\nHere's the input: ", clustering, "\n")

    avg_sil = [] # List that will contain silhouette score for each cluster
    # Works for cluster 0, cluster 1 having some issues
    for i, cluster in enumerate(clustering):

        print("\n\nCalculating silhouette score for cluster ", i, "\n")

        # Create copy of clustering to safely remove current cluster from separation calculation - looks right
        other_clusters = list(clustering)
        print("Updated clustering (should be same as input): ", other_clusters)
        other_clusters.pop(i)
        print("\nUpdated other_clusters (should be all clusters but i): ", other_clusters)

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
            print("\nAverage cohesion: ", avg_cohesion, "\n")

            # Calculate site_a's average distance to all other clusters - looks right if there's no overlap?
            sep_list = []
            for OC in other_clusters:
                sep_holder = []
                for site_c in OC:
                    sep_holder.append(compute_similarity(site_a, site_c))   #site-by-site similarity
                sep_list.append(sum(sep_holder)/len(sep_holder))            #average similarity from site to clusters
                print("Calculated separation values - site-by-site: ", sep_holder)
            print("Calculated separation values - averaged by cluster", sep_list,)

            # Find the maximum similarity from the focal site to another cluster
            separation = min(sep_list)
            print("Minimum separation from other clusters: ", separation)

            # Calculate the silhouette coefficient for focal site - working
            silhouette.append((avg_cohesion - separation)/max(avg_cohesion, separation))
            print("silhouette score for active site %r: %r\n\n" %(site_a.name, silhouette[j]))

        # Calculate the average silhouette coefficient for focal cluster
        avg_sil.append(sum(silhouette)/len(silhouette))
        print("silhouette score for cluster %r: %r\n\n" %(i, avg_sil[i]))

    # Return the average silhouette coefficient for the entire clustering
    print("Total cluster silhouette score: ", sum(avg_sil)/len(avg_sil))
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

    similarity = 0.0

    # Only keep the residue names, not their numbers
    a_residues = [str(x).split(' ')[0] for x in site_a.residues]
    b_residues = [str(x).split(' ')[0] for x in site_b.residues]

    # Simple Jaccard similarity
    intersection = set.intersection(*[set(a_residues), set(b_residues)])
    union = set.union(*[set(a_residues), set(b_residues)])

    # Hard-code BLOSUM matrix to penalize the appearance of x instead of y
    #BLOSUM = {}

    similarity = len(intersection)/len(union)

    return similarity
