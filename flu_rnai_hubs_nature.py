""" Are the 295 genes from the chandra 2010
    nature paper enriched in HPRD hubs?
"""
import utils_graph, utils_stats

rnai_file = '../Thesis/Data/Network/Flu/nature2010/rnai_hits'
hubs_file = '../Thesis/Data/Hubs2/HPRD.entrez.expand.hubs20'
net_file = '../Thesis/Data/Network/Human/HPRD/hprd_new.intr.ls.entrez'
party_file = '../Thesis/Data/Hubs2/2.party'
date_file = '../Thesis/Data/Hubs2/2.date'

rnai_genes = set(utils_graph.getNodes(rnai_file))
network_genes = set(utils_graph.getNodes(net_file))
hubs = set(utils_graph.getNodes(hubs_file))
party = set(utils_graph.getNodes(party_file))
date = set(utils_graph.getNodes(date_file))

print len(party & rnai_genes), len(date & rnai_genes)


