import sys
import random
from .io import read_active_sites, write_clustering, write_mult_clusterings, read_active_site
from .cluster import cluster_by_partitioning, cluster_hierarchically
from .metrics import calc_silhouette, seq_scan, compute_similarity, compute_rand

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
    print("Usage: python -m cluster [-P| -H] <pdb directory> <output file>")
    sys.exit(0)

active_sites = read_active_sites(sys.argv[2])

# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    clustering = cluster_by_partitioning(active_sites)
    write_clustering(sys.argv[3], clustering)

if sys.argv[1][0:2] == '-H':
    print("Clustering using hierarchical method")
    clusterings = cluster_hierarchically(active_sites)
    write_mult_clusterings(sys.argv[3], clusterings)

# # # # # # # # SCRATCHWORK # # # # # # # #

### Create edge file for input into Cytoscape
# This is a file with three columns:
#   column 1 = ActiveSiteA
#   column 2 = ActiveSiteB
#   column 3 = edge weight (similarity score)
#path = "/Users/student/Documents/BMI206/bmi203-2/data/"
#active_sites = read_active_sites(path)
#with open('PDB_edges.tsv','w') as f:
#    for site_a in active_sites:
#        for site_b in active_sites:
#            sim = compute_similarity(site_a, site_b)
#            line = "%r\t%r\t%f\n" %(site_a.name, site_b.name, sim)
#            print(line)
#            print(line, file = f)


### Calculate ideal percentile cutoff score for this dataset (for partitioning)
#import numpy as np
#for percentile in np.arange(0,20,0.5):
#    avg = []
#    print("Percentile: ", percentile)
#    for i in range(0,10):
#        clustering = cluster_by_partitioning(active_sites, percentile)
#        avg.append(calc_silhouette(clustering))
#    print("Average score for ", percentile, "percentile cutoff: ", sum(avg)/len(avg))

### For ten trials, run two different instances of the partitioning function, calculate Rand
#rands = []
#for i in range(0,10):
#    clust_a = cluster_by_partitioning(active_sites)
#    clust_b = cluster_by_partitioning(active_sites)
#    rands.append(compute_rand(active_sites, clust_a, clust_b))

#print("Average Rand Index: ", sum(rands)/len(rands))
#print("Standard Deviation: ", np.std(rands, axis = 0))

### Create Cytoscape node file for the partitioning algorithm
#clustering = cluster_by_partitioning(active_sites)
#with open('partitioning_node_assignments.tsv','w') as f:
#    print("Active_Site\tCluster\n", file = f)
#    for i, cluster in enumerate(clustering):
#        for active_site in cluster:
#            line = "%r\t%f\n" %(active_site.name, i)
#            print(line, file = f)

### Create Cytoscape node file for the h_clust algorithm
### Run twice, the second time uncommenting the reverse sort command
#clustering = cluster_hierarchically(active_sites)
#with open('hierarchical_node_assignments.tsv','w') as f:
#with open('rev_hierarchical_node_assignments.tsv','w') as f:
#    print("Active_Site\tCluster\n", file = f)
#    for i, cluster in enumerate(clustering):
#        for active_site in cluster:
#            line = "%r\t%f\n" %(active_site.name, i)
#            print(line, file = f)

### Find silhouette scores for best p_clust and h_clust
### For ten trials, run two different instances of the partitioning function, calculate Rand
#h = []
#p = []
#for i in range(0,10):
#    clust_p = cluster_by_partitioning(active_sites)
#    p.append(calc_silhouette(clust_p))
#    clust_h = cluster_hierarchically(active_sites)
#    h.append(calc_silhouette(clust_h))

#print("Average Rand Index Hier: ", sum(h)/len(h))
#print("Standard Deviation Hier: ", np.std(h, axis = 0))

#print("Average Rand Index Part: ", sum(p)/len(p))
#print("Standard Deviation Part: ", np.std(p, axis = 0))

# Get Rand Index comparing
#clust_a = cluster_by_partitioning(active_sites)
#clust_b = cluster_hierarchically(active_sites)
#print(compute_rand(active_sites, clust_a, clust_b))
