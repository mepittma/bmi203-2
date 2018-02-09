import sys
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import cluster_by_partitioning, cluster_hierarchically
from .metrics import calc_silhouette, seq_scan, compute_similarity, compute_rand

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




# # # # # # # # # # # # test rand index # # # # # # # # # # # # # #

# Read in all the active sites to cluster
path = "/Users/student/Documents/BMI206/bmi203-2/data/"
active_sites = read_active_sites(path)

p_clust = cluster_by_partitioning(active_sites)
print("Partitioning cluster: ", p_clust)

h_clust = cluster_hierarchically(active_sites)
print("Hierarchical cluster: ", h_clust)

print("Rand index of two clusterings: ", compute_rand(active_sites, p_clust, h_clust))
