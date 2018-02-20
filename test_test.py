from cluster import cluster
from cluster import io
from cluster import utils
import os
import random

# tractable subset
pdb_ids = [276, 4629, 10701]

active_sites = []
for id in pdb_ids:
    filepath = os.path.join("data", "%i.pdb"%id)
    active_sites.append(io.read_active_site(filepath))

print(active_sites)

# Since covering algorithm is subject to initial conditions, set random seed
random.seed(1)

answer = [[utils.ActiveSite("4629"),utils.ActiveSite("10701"),utils.ActiveSite("276")]]

# update this assertion
print("Answer being tested: ", answer)
print("Function returns: ", cluster.cluster_by_partitioning(active_sites))
print("Type of first element in answer: ", type(answer[0]))
print("Type of first element in first list: ", type(answer[0][0]))
print("Type of first element in partitioning: ", type(cluster.cluster_by_partitioning(active_sites)[0]))
print("Type of first element in first list of partitioning: ", type(cluster.cluster_by_partitioning(active_sites)[0][0]))
assert cluster.cluster_by_partitioning(active_sites) == [[utils.ActiveSite("4629"),utils.ActiveSite("10701"),utils.ActiveSite("276")]]
