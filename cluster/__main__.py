import sys
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import cluster_by_partitioning, cluster_hierarchically
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




# # # # # # # # # # # # test simple similarity index # # # # # # # # # # # # # #

from .io import read_active_site
import os, random
import numpy as np
from collections import defaultdict

# Test Jaccard similiarity index on two random files
path = "/Users/student/Documents/BMI206/bmi203-2/data/"
actsite1 = read_active_site(path + random.choice(os.listdir(path)))
actsite2 = read_active_site(path + random.choice(os.listdir(path)))
print("Similarity: ", (compute_similarity(actsite1, actsite2)))

# Read in all the active sites and partition them
active_sites = read_active_sites(path)

# # # # # # # # SUBSET OF PARTITION FUNCTION # # # # # # # #
# 1. Initialization - shuffle the order of active sites
random.shuffle(active_sites)

# 2. Threshold assessment - calculate a prior distribution of distances, selected
# randomly from the pairwise distances between active sites (n-squared)
dist = []
for i in range(0,len(active_sites) ** 2):

    # select two random active sites
    sites = random.sample(active_sites, 2)

    # calculate similarity, print statements to confirm it's happening
    sim = compute_similarity(sites[0], sites[1])
    #print ("\nJaccard similarity between %r and %r: %r" %(sites[0].name, sites[1].name, sim))

    # add to distribution
    dist.append(sim)

# Print the distribution just for funsies
test_dist = [1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000]
print("Test threshold: ", np.percentile(test_dist,10))

# calculate threshold distance value (lower 10th percentile of similarity distribution)
# any two comparisons below this similarity score will be considered "too distant" to cluster together
t = np.percentile(dist,10)
print("Threshold value: ", t)

# 3. Cluster assignment - seq through the active site sequences and cluster
# based on similarity and threshold cutoff
clusters = defaultdict(list)
for site in active_sites:
    assigned = False
    #print("\n\nAssigning %r to a cluster..." %site.name)

    # try to assign to an existing cluster
    for c in clusters:
        #print("Key: ",c)

        #print("Type of cluster object: ", type(clusters[c]))
        #print("Type of key: ", type(c))

        # Find distribution of pairwise distances between already-assigned sites to the query site
        sim_dist = []
        for seen_site in clusters[c]:
            sim_dist.append(compute_similarity(site, seen_site))

        # If median similarity of the query sequence is greater than the threshold
        # cutoff, add it to the cluster; else restart the loop
        #print("Median distance between %r and cluster %r: %r" %(site.name, c, np.median(sim_dist)))
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

# Now check that metrics.py is working right
clustering = seq_scan(active_sites, defaultdict(list), t) # success!

print("Here's what came out of native code: ", clusters)
print("Here's what came out of the function: ", clustering)

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
# # # # # # # # SUBSET OF PARTITION FUNCTION # # # # # # # #

#cluster_by_partitioning(active_sites)


# Read in all the active sites and hierarchically cluster them
