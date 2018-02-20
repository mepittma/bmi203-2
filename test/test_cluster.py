from cluster import cluster
from cluster import io
from cluster import utils
import os
import random

def test_similarity():
    filename_a = os.path.join("data", "276.pdb")
    filename_b = os.path.join("data", "4629.pdb")

    activesite_a = io.read_active_site(filename_a)
    activesite_b = io.read_active_site(filename_b)

    # update this assertion
    assert cluster.compute_similarity(activesite_a, activesite_b) == float(1/8)

    # Calculation: number in common = {His Asp} vs {Asp Thr Arg Ser Lys Tyr Asn}
    # Common = 1, Total unique = 8

def test_partition_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # Since covering algorithm is subject to initial conditions, set random seed
    random.seed(1)

    # update this assertion
    assert cluster.cluster_by_partitioning(active_sites) == [[utils.ActiveSite("4629"),utils.ActiveSite("10701"),utils.ActiveSite("276")]]

def test_hierarchical_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    assert cluster.cluster_hierarchically(active_sites) == [[[utils.ActiveSite("276")], [utils.ActiveSite("4629"), utils.ActiveSite("10701")]], [[utils.ActiveSite("276"), utils.ActiveSite("4629"), utils.ActiveSite("10701")]]]
