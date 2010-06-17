"""At this point, I am questioning the
   conclusion that RNAi data is enriched
   in hubs. I have shown that HCV RNAi
   hits that are genes that are bad
   for the virus (probably immune genes)
   are not enriched in hubs.

   I want to test the same thing for flu.
   This is complicated by the Z-scores
   that are provided in the Cell article
   w/ Aviv Regev as co author. I'll have
   to pick an arbitrary cutoff, and look
   at hub enrichment compared to all hits
   in the file.  All hits is not the best
   way to go; really I need to ~1700 genes
   that were tested & could be knocked
   down w/o killing the cell, but these
   are not available.
"""
import utils_stats, utils_graph

flu_rnai_file = '../Thesis/Data/Network/Flu/cell09/all_rnai'
hubs_file = '../Thesis/Data/Hubs2/HPRD.entrez.expand.hubs20'
net_file = '../Thesis/Data/Network/Human/HPRD/hprd_new.intr.ls.entrez'

network_genes = set(utils_graph.getNodes(net_file))
hubs = set(utils_graph.getNodes(hubs_file))

all_rnai = {}
replication_rnai = {}
with open(flu_rnai_file) as f:
    for line in f:
        [entrez, delNS1, 
         vRNA, replication] = [float(x) for x in line.strip().split('\t')]
        ID = str(int(entrez))
        if vRNA > float(1):
            replication_rnai[ID] = True
        all_rnai[ID] = True

bg = set(network_genes & set(all_rnai.keys()))
rep_set = set(replication_rnai.keys())
print 'background', len(bg)
print 'hubs in background', len(bg & hubs)
print 'replication bg', len(rep_set & bg)
print 'replication hubs', len(hubs & rep_set)

rep_nonHubs = set(set(rep_set & bg) - hubs)
bg_hubs = set(bg & hubs) - rep_set
bg_nonHubs = set(bg - hubs) - rep_nonHubs


print len(bg_hubs), len(bg_nonHubs)
print len(hubs & rep_set), len(rep_nonHubs)
print utils_stats.fisher_positive_pval([len(bg_hubs), len(bg_nonHubs)],
                                       [len(hubs & rep_set),
                                        len(rep_nonHubs)])
