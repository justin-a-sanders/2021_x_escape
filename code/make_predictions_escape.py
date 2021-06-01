import os
import csv
import json
import sys
import pysam
os.environ["CUDA_VISIBLE_DEVICES"] = '-1' ### run on CPU
import tensorflow as tf
print(tf.__version__)
if tf.__version__[0] == '1':
    tf.compat.v1.enable_eager_execution()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cooltools.lib.numutils import set_diag
from basenji import dataset, dna_io, seqnn

def from_upper_triu(vector_repr, matrix_len, num_diags):
    z = np.zeros((matrix_len,matrix_len))
    triu_tup = np.triu_indices(matrix_len,num_diags)
    z[triu_tup] = vector_repr
    for i in range(-num_diags+1,num_diags):
        set_diag(z, np.nan, i)
    return z + z.T

fasta_open = pysam.Fastafile(sys.argv[1])
escape_genes = sys.argv[2]

### load params, specify model ###
model_dir = './'
params_file = model_dir+'params.json'
model_file  = model_dir+'model_best.h5'
with open(params_file) as params_open:
    params = json.load(params_open)
    params_model = params['model']
    params_train = params['train']

seqnn_model = seqnn.SeqNN(params_model)

seqnn_model.restore(model_file)
print('successfully loaded')

data_dir = './data/'

hic_targets = pd.read_csv(data_dir+'/targets.txt',sep='\t')
hic_file_dict_num = dict(zip(hic_targets['index'].values, hic_targets['file'].values) )
hic_file_dict     = dict(zip(hic_targets['identifier'].values, hic_targets['file'].values) )
hic_num_to_name_dict = dict(zip(hic_targets['index'].values, hic_targets['identifier'].values) )

# read data parameters
data_stats_file = '%s/statistics.json' % data_dir
with open(data_stats_file) as data_stats_open:
    data_stats = json.load(data_stats_open)
seq_length = data_stats['seq_length']
target_length = data_stats['target_length']
hic_diags =  data_stats['diagonal_offset']
target_crop = data_stats['crop_bp'] // data_stats['pool_width']
target_length1 = data_stats['seq_length'] // data_stats['pool_width']

#------------------------------------------------------------------------------#

sequences = pd.read_csv(data_dir+'sequences.bed',sep='\t',  names=['chr','start','stop','type'])
seqs_per_tf_default = 256
test_tr_num = 0
sequences_test = sequences.iloc[  sequences['type'].values=='test']
sequences_test.reset_index(inplace=True, drop=True)

#------------------------------------------------------------------------------#

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])

import subprocess
if not os.path.isfile('./data/hg38.ml.fa'):
    print('downloading hg38.ml.fa')
    subprocess.call('curl -o ./data/hg38.ml.fa.gz https://storage.googleapis.com/basenji_barnyard/hg38.ml.fa.gz', shell=True)
    subprocess.call('gunzip ./data/hg38.ml.fa.gz', shell=True)

target_length1_cropped = target_length1 - 2*target_crop
print('flattened representation length:', target_length)
print('symmetrix matrix size:', '('+str(target_length1_cropped)+','+str(target_length1_cropped)+')')



################################################################################
                            # Run Model Here #
################################################################################

# fasta_open = pysam.Fastafile('./data/hg38.ml.fa')  # human
# fasta_open = pysam.Fastafile('./data/chrX.SNPs_introduced_hap1.fa')  # GM12878 parental allele 1
# fasta_open = pysam.Fastafile('./data/chrX.SNPs_introduced_hap2.fa')  # GM12878 parental allele 2
# fasta_open = pysam.Fastafile('./data/chrX.fa') # mouse

genes = []
with open(escape_genes) as fd:
    rd = csv.reader(fd, delimiter="\t")
    for row in rd:
        genes.append((row[3], int(row[1]), int(row[2])))

gene_to_motifs = {}
with open("peak_motifs.bed") as fd:
    rd = csv.reader(fd, delimiter="\t")
    for row in rd:
        pos = int(row[1])
        for name, start, end in genes:
            if pos > start - 100000 and pos < end + 100000:
                # print(name, int(row[1]), int(row[2]))
                if name in gene_to_motifs:
                    gene_to_motifs[name].append((int(row[1]), int(row[2])))
                else:
                    gene_to_motifs[name] = [(int(row[1]), int(row[2]))]

target_index = 2
for name, start, end in genes:
    print(name)
    os.makedirs("predictions/" + name, exist_ok=True)
    seq_start = start - 450000 - 2**19
    seq_end = seq_start + 2*2**20
    seq = fasta_open.fetch('chrX', seq_start, seq_end).upper()
    s1 = dna_io.dna_1hot(seq[:2**20])
    s2 = dna_io.dna_1hot(seq[2**19:2**19+2**20])
    s3 = dna_io.dna_1hot(seq[2**20:])
    wt_pred_from_seq1 = seqnn_model.model.predict(np.expand_dims(s1,0))
    wt_pred_from_seq2 = seqnn_model.model.predict(np.expand_dims(s2,0))
    wt_pred_from_seq3 = seqnn_model.model.predict(np.expand_dims(s3,0))
    wt_pred_from_seq1 = from_upper_triu(wt_pred_from_seq1[:,:,target_index], target_length1_cropped, hic_diags)
    wt_pred_from_seq2 = from_upper_triu(wt_pred_from_seq2[:,:,target_index], target_length1_cropped, hic_diags)
    wt_pred_from_seq3 = from_upper_triu(wt_pred_from_seq3[:,:,target_index], target_length1_cropped, hic_diags)
    wt_pred_from_seq = np.zeros((960,960))
    wt_pred_from_seq[0:448,0:448] = wt_pred_from_seq1
    wt_pred_from_seq[512:,512:] = wt_pred_from_seq3
    wt_pred_from_seq[256:704,256:704] = wt_pred_from_seq2
    np.save("predictions/" + name + "/" + name + "_unchanged", wt_pred_from_seq)

    for motif_start, motif_end in gene_to_motifs[name]:
        # predict effect of inverision
        motif_end += 100
        motif_start -= 100
        motif_len = motif_end - motif_start
        rev_seq = fasta_open.fetch('chrX', seq_start, seq_end).upper()
        #print("    ", motif_start)
        #print("    ", rev_seq[motif_start-seq_start:motif_start-seq_start+motif_len], reverse_complement(rev_seq[motif_start-seq_start:motif_start-seq_start+motif_len]))
        rev_seq = rev_seq[:motif_start-seq_start] + \
                    reverse_complement(rev_seq[motif_start-seq_start:motif_start-seq_start+motif_len]) + \
                    rev_seq[motif_start-seq_start+motif_len:]
        s1 = dna_io.dna_1hot(rev_seq[:2**20])
        s2 = dna_io.dna_1hot(rev_seq[2**19:2**19+2**20])
        s3 = dna_io.dna_1hot(rev_seq[2**20:])
        wt_pred_from_seq1 = seqnn_model.model.predict(np.expand_dims(s1,0))
        wt_pred_from_seq2 = seqnn_model.model.predict(np.expand_dims(s2,0))
        wt_pred_from_seq3 = seqnn_model.model.predict(np.expand_dims(s3,0))
        wt_pred_from_seq1 = from_upper_triu(wt_pred_from_seq1[:,:,target_index], target_length1_cropped, hic_diags)
        wt_pred_from_seq2 = from_upper_triu(wt_pred_from_seq2[:,:,target_index], target_length1_cropped, hic_diags)
        wt_pred_from_seq3 = from_upper_triu(wt_pred_from_seq3[:,:,target_index], target_length1_cropped, hic_diags)
        wt_pred_from_seq = np.zeros((960,960))
        wt_pred_from_seq[0:448,0:448] = wt_pred_from_seq1
        wt_pred_from_seq[512:,512:] = wt_pred_from_seq3
        wt_pred_from_seq[256:704,256:704] = wt_pred_from_seq2
        np.save("predictions/" + name + "/" + name + "_" + str(motif_start) + "_reversed", wt_pred_from_seq)

        #predict effect of deletion
        del_seq = fasta_open.fetch('chrX', seq_start, seq_end + motif_len).upper()
        del_seq = del_seq[:motif_start-seq_start] + \
                    del_seq[motif_start-seq_start+motif_len:]
        s1 = dna_io.dna_1hot(del_seq[:2**20])
        s2 = dna_io.dna_1hot(del_seq[2**19:2**19+2**20])
        s3 = dna_io.dna_1hot(del_seq[2**20:])
        wt_pred_from_seq1 = seqnn_model.model.predict(np.expand_dims(s1,0))
        wt_pred_from_seq2 = seqnn_model.model.predict(np.expand_dims(s2,0))
        wt_pred_from_seq3 = seqnn_model.model.predict(np.expand_dims(s3,0))
        wt_pred_from_seq1 = from_upper_triu(wt_pred_from_seq1[:,:,target_index], target_length1_cropped, hic_diags)
        wt_pred_from_seq2 = from_upper_triu(wt_pred_from_seq2[:,:,target_index], target_length1_cropped, hic_diags)
        wt_pred_from_seq3 = from_upper_triu(wt_pred_from_seq3[:,:,target_index], target_length1_cropped, hic_diags)
        wt_pred_from_seq = np.zeros((960,960))
        wt_pred_from_seq[0:448,0:448] = wt_pred_from_seq1
        wt_pred_from_seq[512:,512:] = wt_pred_from_seq3
        wt_pred_from_seq[256:704,256:704] = wt_pred_from_seq2
        np.save("predictions/" + name + "/" + name + "_" + str(motif_start) + "_deletion", wt_pred_from_seq)
