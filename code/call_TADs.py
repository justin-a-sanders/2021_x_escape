import os
import csv
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
import math
from scipy import signal

escape_genes = sys.argv[1]
peak_motifs = sys.argv[2]

def insulation_scores(mat):
    mat = 2**mat
    window_size = 100
    scores = []
    for i in range(window_size, len(mat) - window_size):
        diamond = mat[i-window_size:i, i:i+window_size]
        scores.append(np.mean(diamond))
    mean_score = np.nanmean(scores)
    scores = [math.log2(s/mean_score) for s in scores]
    return ([None] * window_size) + scores + ([None] * window_size)

vmin=-1; vmax=1

genes = []
with open(escape_genes) as fd:
    rd = csv.reader(fd, delimiter="\t")
    for row in rd:
        genes.append((row[3], int(row[1]), int(row[2]), row[5]))


gene_to_motifs = {}
with open(peak_motifs) as fd:
    rd = csv.reader(fd, delimiter="\t")
    for row in rd:
        pos = int(row[1])
        for name, start, end, _ in genes:
            if pos > start - 100000 and pos < end + 100000:
                if name in gene_to_motifs:
                    gene_to_motifs[name].append((int(row[1]), int(row[2])))
                else:
                    gene_to_motifs[name] = [(int(row[1]), int(row[2]))]

prelude = """<html>
<head>
<style>
h2 {text-align: center;}
h3 {text-align: center;}
tr {text-align: center;}
p {padding-left: 25px; padding-bottom:35px}
table {margin-left: auto;  margin-right: auto;}
</style>
</head>
<body>"""

row = """<tr>
<td><h4>{0}</h4></td>
<td><img src='{1}' height=110 width=540></td>
<td><img src='{2}' height=110 width=540></td>
</tr>
"""

of = open("tad_predictions.html", "w")
of.write(prelude)
for name, start, end, orientation in genes:
    min1, max1, min2, max2 = -0.5, 0.5, -0.15, 0.15
    print(name)
    os.makedirs("tads/" + name, exist_ok=True)
    unchanged = np.load("predictions/" + name + "/" + name + "_unchanged.npy")
    f, ax2 = plt.subplots(figsize=(13,6))
    gs = grd.GridSpec(2, 3, height_ratios=[8,1], width_ratios=[12,1.8, 12], wspace=0.18)
    ax = plt.subplot(gs[0])
    p = ax.matshow(unchanged, cmap= 'RdBu_r', vmax=vmax, vmin=vmin, aspect='auto')
    plt.ylabel('chrX: ' + str(start-450000 - 500000) + "-" + str(start-450000-500000+2*2**20))
    plt.title(name + ' unchanged')
    gene_start_bin = int(450000/2048) + 256 - 32
    gene_end_bin = int((450000 + end - start)/2048) + 256 - 32
    if orientation == '+':
        plt.xticks(ticks=[gene_start_bin, gene_end_bin], labels=["TSS     ", "     TES"])
    else:
        plt.xticks(ticks=[gene_start_bin, gene_end_bin], labels=["TES     ", "     TSS"])

    colorAx = plt.subplot(gs[1])
    cb = plt.colorbar(p, cax = colorAx)

    ax = plt.subplot(gs[2])
    p = ax.matshow(unchanged[gene_start_bin-50:gene_end_bin+50,gene_start_bin-50:gene_end_bin+50], cmap= 'RdBu_r', vmax=vmax, vmin=vmin, aspect='auto')
    plt.title("Zoomed in " + name + ' unchanged')
    ax.yaxis.tick_right()
    adj_gene_end_bin = gene_end_bin - gene_start_bin + 50
    adj_gene_start_bin = 50
    if orientation == '+':
        plt.xticks(ticks=[adj_gene_start_bin, adj_gene_end_bin], labels=["TSS     ", "     TES"])
    else:
        plt.xticks(ticks=[adj_gene_start_bin, adj_gene_end_bin], labels=["TES     ", "     TSS"])

    ax2 = plt.subplot(gs[3])
    insulation = insulation_scores(np.nan_to_num(unchanged))
    boundaries = signal.find_peaks([-i for i in insulation[100:-100]], prominence = 0.1)
    ax2.plot(insulation)
    for idx in boundaries[0]+100:
        ax2.plot([idx,idx], [min1,max1], c='#2ca02c')
    plt.xlim(0,len(insulation))
    plt.ylim(min1,max1)
    plt.savefig("tads/" + name + "/" + name + "_unchanged_TADs_mat.png", dpi=800,bbox_inches='tight')
    #plt.clf()
    unchanged_insulation = insulation

    header = "\n <h2> " + name + """ </h2>
    <table> <tr> <td> """ + "<img src='" + "tads/" + name + "/" + name + "_unchanged_TADs_mat.png" + "' height=530 width=870>" + """</td> </tr> </table>
    <table> <tr>
    <td><h3>Motif change</h3></td>
    <td><h3>TAD boundaries</h3></td>
    <td><h3>Delta Insulation Score</h3></td>
    </tr> \n"""
    of.write(header)

    # gs = grd.GridSpec(len(gene_to_motifs[name]) * 2 + 1, 1)
    # gs.update(wspace = 1.5, hspace = 1.5)
    # ax2 = plt.subplot(gs[0])
    # insulation = insulation_scores(np.nan_to_num(unchanged))
    # boundaries = signal.find_peaks([-i for i in insulation[100:-100]], prominence = 0.1)
    # ax2.plot(insulation)
    # for idx in boundaries[0]+100:
    #     ax2.plot([idx,idx], [min(insulation[100:-100]),max(insulation[100:-100])], c='#2ca02c')
    # plt.xlim(0,len(insulation))
    # plt.title('unchanged', fontdict={'fontsize':6})
    # ax2.set_xticks([], [])
    # ax2.set_yticks([], [])

    f, ax2 = plt.subplots(figsize=(5,0.85))
    ax2.plot(insulation)
    for idx in boundaries[0]+100:
        ax2.plot([idx,idx], [min1,max1], c='#2ca02c')
    plt.xlim(0,len(insulation))
    plt.ylim(min1,max1)
    plt.title(name + ' unchanged', fontdict={'fontsize':7})
    ax2.set_xticks([], [])
    #ax2.set_yticks([], [])
    f.tight_layout()
    if orientation == '+':
        plt.xticks(ticks=[gene_start_bin, gene_end_bin], labels=["TSS     ", "     TES"], fontsize=7)
    else:
        plt.xticks(ticks=[gene_start_bin, gene_end_bin], labels=["TES     ", "     TSS"], fontsize=7)
    plt.savefig("tads/" + name + "/" + name + "_unchanged_TADs.png", dpi=800,bbox_inches='tight')
    #plt.clf()

    f, ax2 = plt.subplots(figsize=(5,0.85))
    ax2.plot([None] * 100 + [a-b for (a,b) in zip(insulation,unchanged_insulation) if (a != None and b != None)])
    plt.xlim(0,len(insulation))
    plt.ylim(min2,max2)
    plt.title(name + ' unchanged difference', fontdict={'fontsize':7})
    ax2.set_xticks([], [])
    f.tight_layout()
    if orientation == '+':
        plt.xticks(ticks=[gene_start_bin, gene_end_bin], labels=["TSS     ", "     TES"], fontsize=7)
    else:
        plt.xticks(ticks=[gene_start_bin, gene_end_bin], labels=["TES     ", "     TSS"], fontsize=7)
    plt.savefig("tads/" + name + "/" + name + "_unchanged_TADs_difference.png", dpi=800,bbox_inches='tight')
    plt.clf()

    of.write(row.format("unchanged", "tads/" + name + "/" + name + "_unchanged_TADs.png", "tads/" + name + "/" + name + "_unchanged_TADs_difference.png"))


    for i, (motif_start, motif_end) in enumerate(gene_to_motifs[name]):
        deleted = np.load("predictions/" + name + "/" + name + "_" + str(motif_start) + "_deletion.npy")
        reversed = np.load("predictions/" + name + "/" + name + "_" + str(motif_start) + "_reversed.npy")

        f, ax2 = plt.subplots(figsize=(5,1))
        insulation = insulation_scores(np.nan_to_num(deleted))
        boundaries = signal.find_peaks([-i for i in insulation[100:-100]], prominence = 0.1)
        ax2.plot(insulation)
        for idx in boundaries[0]+100:
            ax2.plot([idx,idx], [min1,max1], c='#2ca02c')
        ax2.plot([int((motif_start-(start-450000))/2048) + 256 - 32,int((motif_start-(start-450000))/2048) + 256 - 32], [min1,max1], c='#242322', linewidth=1)
        plt.xlim(0,len(insulation))
        plt.ylim(min1,max1)
        plt.title(name + " " + str(motif_start) + ' deleted', fontdict={'fontsize':7})
        if orientation == '+':
            plt.xticks(ticks=[gene_start_bin, gene_end_bin, int((motif_start-(start-450000))/2048) + 256 - 32], labels=["TSS     ", "     TES",""], fontsize=7)
        else:
            plt.xticks(ticks=[gene_start_bin, gene_end_bin, int((motif_start-(start-450000))/2048) + 256 - 32], labels=["TES     ", "     TSS",""], fontsize=7)
        f.tight_layout()
        plt.savefig("tads/" + name + "/" + name + "_" + str(motif_start) + "_del_TADs.png", dpi=800,bbox_inches='tight')
        #plt.clf()

        f, ax2 = plt.subplots(figsize=(5,1))
        insulation = insulation_scores(np.nan_to_num(deleted))
        boundaries = signal.find_peaks([-i for i in insulation[100:-100]], prominence = 0.1)
        ax2.plot([None] * 100 + [a-b for (a,b) in zip(insulation,unchanged_insulation) if (a != None and b != None)])
        ax2.plot([int((motif_start-(start-450000))/2048) + 256 - 32,int((motif_start-(start-450000))/2048) + 256 - 32], [min1,max1], c='#242322', linewidth=1)
        plt.xlim(0,len(insulation))
        plt.ylim(min2,max2)
        plt.title(name + " " + str(motif_start) + ' deleted difference', fontdict={'fontsize':7})
        if orientation == '+':
            plt.xticks(ticks=[gene_start_bin, gene_end_bin, int((motif_start-(start-450000))/2048) + 256 - 32], labels=["TSS     ", "     TES",""], fontsize=7)
        else:
            plt.xticks(ticks=[gene_start_bin, gene_end_bin, int((motif_start-(start-450000))/2048) + 256 - 32], labels=["TES     ", "     TSS",""], fontsize=7)
        #ax2.set_yticks([], [])
        f.tight_layout()
        plt.savefig("tads/" + name + "/" + name + "_" + str(motif_start) + "_del_TADs_difference.png", dpi=800,bbox_inches='tight')
        plt.clf()

        of.write(row.format(str(motif_start) + " deleted", "tads/" + name + "/" + name + "_" + str(motif_start) + "_del_TADs.png", "tads/" + name + "/" + name + "_" + str(motif_start) + "_del_TADs_difference.png"))

        f, ax2 = plt.subplots(figsize=(5,1))
        insulation = insulation_scores(np.nan_to_num(reversed))
        boundaries = signal.find_peaks([-i for i in insulation[100:-100]], prominence = 0.1)
        ax2.plot(insulation)
        for idx in boundaries[0]+100:
            ax2.plot([idx,idx], [min1,max1], c='#2ca02c')
        ax2.plot([int((motif_start-(start-450000))/2048) + 256 - 32,int((motif_start-(start-450000))/2048) + 256 - 32], [min1,max1], c='#242322', linewidth=1)
        plt.xlim(0,len(insulation))
        plt.ylim(min1,max1)
        plt.title(name + " " + str(motif_start) + ' reversed', fontdict={'fontsize':7})
        if orientation == '+':
            plt.xticks(ticks=[gene_start_bin, gene_end_bin, int((motif_start-(start-450000))/2048) + 256 - 32], labels=["TSS     ", "     TES",""], fontsize=7)
        else:
            plt.xticks(ticks=[gene_start_bin, gene_end_bin, int((motif_start-(start-450000))/2048) + 256 - 32], labels=["TES     ", "     TSS",""], fontsize=7)
        # ax2.set_yticks([], [])
        f.tight_layout()
        plt.savefig("tads/" + name + "/" + name + "_" + str(motif_start) + "_rev_TADs.png", dpi=800,bbox_inches='tight')
        plt.clf()

        f, ax2 = plt.subplots(figsize=(5,1))
        insulation = insulation_scores(np.nan_to_num(reversed))
        boundaries = signal.find_peaks([-i for i in insulation[100:-100]], prominence = 0.1)
        ax2.plot([None] * 100 + [a-b for (a,b) in zip(insulation,unchanged_insulation) if (a != None and b != None)])
        ax2.plot([int((motif_start-(start-450000))/2048) + 256 - 32,int((motif_start-(start-450000))/2048) + 256 - 32], [min1,max1], c='#242322', linewidth=1)
        plt.xlim(0,len(insulation))
        plt.ylim(min2,max2)
        plt.title(name + " " + str(motif_start) + ' reversed difference', fontdict={'fontsize':7})
        if orientation == '+':
            plt.xticks(ticks=[gene_start_bin, gene_end_bin, int((motif_start-(start-450000))/2048) + 256 - 32], labels=["TSS     ", "     TES",""], fontsize=7)
        else:
            plt.xticks(ticks=[gene_start_bin, gene_end_bin, int((motif_start-(start-450000))/2048) + 256 - 32], labels=["TES     ", "     TSS",""], fontsize=7)
        #ax2.set_yticks([], [])
        f.tight_layout()
        plt.savefig("tads/" + name + "/" + name + "_" + str(motif_start) + "_rev_TADs_difference.png", dpi=800,bbox_inches='tight')
        plt.clf()

        of.write(row.format(str(motif_start) + " reversed", "tads/" + name + "/" + name + "_" + str(motif_start) + "_rev_TADs.png", "tads/" + name + "/" + name + "_" + str(motif_start) + "_rev_TADs_difference.png"))

    of.write("</table> \n")
of.close()
