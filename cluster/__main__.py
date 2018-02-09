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




# # # # # # # # # # # # test rand index # # # # # # # # # # # # # #

# Read in all the active sites to cluster
path = "/Users/student/Documents/BMI206/bmi203-2/data/"
active_sites = read_active_sites(path)

p_clust = cluster_by_partitioning(active_sites)
print(p_clust)

h_clust = cluster_hierarchically(active_sites)
print(h_clust)

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
