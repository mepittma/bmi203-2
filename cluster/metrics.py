from .io import read_active_sites, write_clustering, write_mult_clusterings
from .utils import Atom, Residue, ActiveSite
from collections import defaultdict
import numpy as np


def calc_silhouette(clustering):
    """
    This function accepts a clustering assignment (list of lists of clusters) and
    returns a score accounting for how similar objects in a cluster are to each other
    (cohesion) and how distant they are to other clusters (separation).
    """

    #for cluster in clustering:

        # Calculate the cohesion among sites in that cluster
        #for site in cluster:


def seq_scan(active_sites, clusters, t):
    """
    Scans through active site sequences, calculates their similarity, and clusters
    based on similarity to existing sequences in the cluster. If the similarity does
    not meet a previously calculated threshold for any of the existing clusters, then it
    starts a new one.
    """
    for site in active_sites:
        assigned = False

        # try to assign to an existing cluster
        for c in clusters:

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
            else:
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

    print("A residues: ", site_a.residues)
    print("B residues: ", site_b.residues)

    # Only keep the residue names, not their numbers
    a_residues = [str(x).split(' ')[0] for x in site_a.residues]
    b_residues = [str(x).split(' ')[0] for x in site_b.residues]

    print("A residues: ", a_residues)
    print("B residues: ", b_residues)

    # Simple Jaccard similarity
    intersection = set.intersection(*[set(a_residues), set(b_residues)])
    union = set.union(*[set(a_residues), set(b_residues)])

    # Hard-code BLOSUM matrix to penalize the appearance of x instead of y
    #BLOSUM = {}

    similarity = len(intersection)/len(union)

    return similarity
