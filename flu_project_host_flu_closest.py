"""Host and flu share few ELM sequences,
   so I'm making clusters of ELM sequences
   to be able to make comparisons between
   host and flu.

   I want to see how the closest distances
   are distributed.
"""
import Bio.Cluster
import Levenshtein, sys
sys.path.append('../flu/flELM/')
import itertools, sys, os, utils, random, global_settings, numpy
from collections import defaultdict

def get_closest_distances(flu_strings, host_strings):
    """For each flu sequence, find the closet distance to a human string"""

    min_distances = []

    for flu_s in flu_strings:
        dis = Levenshtein.distance(flu_s, 
                                   host_strings[0])
        for host_s in host_strings[1:]:
            d = Levenshtein.distance(flu_s, 
                                     host_s)
            dis = min([dis, d])
        print dis

def count_0s(ls):
    count = 0
    for item in ls:
        if not item:
            count += 1
    return count

def count_flu(protein2counts, all_elmSeqs):
    """Given hits from get_flu_counts, return ELMseq counts"""
    
    counts = defaultdict(utils.init_zero)
    for protein in protein2counts:
        for seq in protein2counts[protein]:
            for elmSeq in protein2counts[protein][seq]:
                elm, sequence = elmSeq.split(':')
                counts[elmSeq] += protein2counts[protein][seq][elmSeq]
                all_elmSeqs[elmSeq] = True
    return counts

def get_flu_counts(afile, proteins):
    """Make protein_name -> seq_name -> elm_seq_counts"""

    counts = {}
    with open(afile) as f:
        for line in f:
            (protein, st, stp,
             elm, seq, junk) = line.strip().split('\t')
            name = protein.split('.')[-1]
            elmSeq = elm + ':' + seq
            if name in proteins:
                if name not in counts:
                    counts[name] = {}
                if protein not in counts[name]:
                    counts[name][protein] = {}
                if elmSeq not in counts:
                    counts[name][protein][elmSeq] = 0
                counts[name][protein][elmSeq] += 1
    return counts

def mk_vec(counts, all_elmSeqs):
    """mk long vector of ELM:seq counts for this host's counts"""
    
    vec = []
    for elmseq in all_elmSeqs:
        if elmseq in counts:
            vec.append(counts[elmseq])
        else:
            vec.append(float(0))
    return vec

def mk_count_vecs(counts, all_elmSeqs):
    """mk long vector of ELM:seq counts for all hosts"""

    vecs = {}
    for host in counts:
        vecs[host] = mk_vec(counts[host],
                            all_elmSeqs)
    return vecs

def mk_count_dists(vecs):
    """change count vectors into distributions"""

    dists = {}
    for host in vecs:
        dists[host] = utils.getDistFromCount(vecs[host])
    return dists

hosts = ('H_sapiens', 'Gallus_gallus')
flus = ('human', 'chicken')
proteins = ('hemagglutinin', 'neuraminidase', 'nucleocapsid protein',
            'matrix protein 1', 'nonstructural protein 1', 'matrix protein 2',
            'nonstructural protein 2', 'polymerase PA', 'polymerase PB2',
            'polymerase PB1', 'PB1-F2 protein')

# count elm:seq occurence
flu_counts = {}
pre_flu_counts = {}
host_counts = {}
all_elmSeqs = {}

for flu in flus:
    pre_flu_counts[flu] = get_flu_counts('../flu/flELM/results/' + flu + '.H5N1.elms', 
                                         proteins)

flu_counts['human'] = count_flu(pre_flu_counts['human'], all_elmSeqs)
flu_counts['chicken'] = count_flu(pre_flu_counts['chicken'], all_elmSeqs)

host_seqs = []
for host in hosts:
    host_counts[host] = defaultdict(utils.init_zero)
    with open('../flu/flELM/results/roundup_all/elmdict_' + host + '.init') as f:
        for line in f:
            (elm, seq, count, fq) = line.strip().split('\t')
            elmSeq = elm + ':' + seq
            host_counts[host][elmSeq] += int(count)
            host_seqs.append(seq)

flu_seqs = []
for flu in flu_counts:
    for elmSeq in flu_counts[flu]:
        elm, seq = elmSeq.split(':')
        flu_seqs.append(seq)

get_closest_distances(flu_seqs, host_seqs)



