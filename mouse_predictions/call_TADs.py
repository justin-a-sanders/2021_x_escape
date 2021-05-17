import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
import math
from scipy import signal

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

genes = [ #human
    ("Ddx3x", 41333308, 41364472, -0.5, 0.5, -0.15, 0.15),
    ("Kdm6a", 44873182, 45112779, -0.5, 0.3, -0.2, 0.2),
    ("Eif2s3", 24054946,24078810, -0.5, 0.5, -0.1, 0.1),
    ("Xist", 73817774, 73852754, -0.5, 0.5, -0.4, 0.4),
    ("Pbdc1", 76173040, 76178314, -0.25, 0.25, -0.15, 0.15),
    ("Kdm5c", 53176277, 53225422, -0.6, 0.5, -0.15, 0.15),
    ("Mid1", 10445310, 10833683, -0.5, 0.5, -0.2, 0.2)
]

genes = [ #mouse
    ("Ddx3x", 13280970, 13294052, -0.5, 0.5, -0.15, 0.15),
    ("Kdm6a", 18162575, 18279936, -0.5, 0.5, -0.15, 0.15),
    ("Eif2s3x", 94188709, 94212651, -0.5, 0.5, -0.15, 0.15),
    ("Xist", 103460373, 103483233, -0.5, 0.5, -0.15, 0.15),
    ("Pbdc1", 105079756, 105117090, -0.5, 0.5, -0.15, 0.15),
    ("Kdm5c", 152233020, 152274535, -0.5, 0.5, -0.15, 0.15),
    ("Mid1", 169685199, 169990798, -0.5, 0.5, -0.15, 0.15),
    ("Car5b", 163976822, 164027997, -0.5, 0.5, -0.15, 0.15)
]


gene_to_motifs = {}
with open("peak_motifs_mm10.bed") as fd:
    rd = csv.reader(fd, delimiter="\t")
    for row in rd:
        pos = int(row[1])
        for name, start, end, _, _, _, _ in genes:
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

of = open("tad_predictions_mouse.html", "w")
of.write(prelude)
for name, start, end, min1, max1, min2, max2 in genes:
    print(name)
    os.makedirs("tads/" + name, exist_ok=True)
    unchanged = np.load("predictions/" + name + "/" + name + "_unchanged.npy")
    f, ax2 = plt.subplots(figsize=(4.5,6))
    gs = grd.GridSpec(2, 2, height_ratios=[8,1], width_ratios=[12,1], wspace=0.05)
    ax = plt.subplot(gs[0])
    p = ax.matshow(unchanged, cmap= 'RdBu_r', vmax=vmax, vmin=vmin, aspect='auto')
    plt.ylabel('chrX: ' + str(start-450000 - 500000) + "-" + str(start-450000-500000+2*2**20))
    plt.title(name + ' unchanged')
    gene_start_bin = int(450000/2048) + 256 - 32
    gene_end_bin = int((450000 + end - start)/2048) + 256 - 32
    plt.xticks(ticks=[gene_start_bin, gene_end_bin], labels=["s", "e"])

    colorAx = plt.subplot(gs[1])
    cb = plt.colorbar(p, cax = colorAx)

    ax2 = plt.subplot(gs[2])
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
    <table> <tr> <td> """ + "<img src='" + "tads/" + name + "/" + name + "_unchanged_TADs_mat.png" + "' height=530 width=510>" + """</td> </tr> </table>
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
    plt.title('unchanged', fontdict={'fontsize':7})
    ax2.set_xticks([], [])
    #ax2.set_yticks([], [])
    f.tight_layout()
    plt.savefig("tads/" + name + "/" + name + "_unchanged_TADs.png", dpi=800,bbox_inches='tight')
    #plt.clf()

    f, ax2 = plt.subplots(figsize=(5,0.85))
    ax2.plot([None] * 100 + [a-b for (a,b) in zip(insulation,unchanged_insulation) if (a != None and b != None)])
    plt.xlim(0,len(insulation))
    plt.ylim(min2,max2)
    plt.title('unchanged difference', fontdict={'fontsize':7})
    ax2.set_xticks([], [])
    #ax2.set_yticks([], [])
    f.tight_layout()
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
        plt.xlim(0,len(insulation))
        plt.ylim(min1,max1)
        plt.title(str(motif_start) + ' deleted', fontdict={'fontsize':7})
        ax2.set_xticks([], [])
        ax2.set_xticks([int((motif_start-(start-450000))/2048) + 256 - 32])
        ax2.set_xticklabels(['m'], fontdict={'fontsize':6})
        #ax2.set_yticks([], [])
        f.tight_layout()
        plt.savefig("tads/" + name + "/" + name + "_" + str(motif_start) + "_del_TADs.png", dpi=800,bbox_inches='tight')
        #plt.clf()

        f, ax2 = plt.subplots(figsize=(5,1))
        insulation = insulation_scores(np.nan_to_num(deleted))
        boundaries = signal.find_peaks([-i for i in insulation[100:-100]], prominence = 0.1)
        ax2.plot([None] * 100 + [a-b for (a,b) in zip(insulation,unchanged_insulation) if (a != None and b != None)])
        plt.xlim(0,len(insulation))
        plt.ylim(min2,max2)
        plt.title(str(motif_start) + ' deleted difference', fontdict={'fontsize':7})
        ax2.set_xticks([], [])
        ax2.set_xticks([int((motif_start-(start-450000))/2048) + 256 - 32])
        ax2.set_xticklabels(['m'], fontdict={'fontsize':6})
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
        plt.xlim(0,len(insulation))
        plt.ylim(min1,max1)
        plt.title(str(motif_start) + ' reversed', fontdict={'fontsize':7})
        ax2.set_xticks([], [])
        ax2.set_xticks([int((motif_start-(start-450000))/2048) + 256 - 32])
        ax2.set_xticklabels(['m'], fontdict={'fontsize':6})
        # ax2.set_yticks([], [])
        f.tight_layout()
        plt.savefig("tads/" + name + "/" + name + "_" + str(motif_start) + "_rev_TADs.png", dpi=800,bbox_inches='tight')
        plt.clf()

        f, ax2 = plt.subplots(figsize=(5,1))
        insulation = insulation_scores(np.nan_to_num(reversed))
        boundaries = signal.find_peaks([-i for i in insulation[100:-100]], prominence = 0.1)
        ax2.plot([None] * 100 + [a-b for (a,b) in zip(insulation,unchanged_insulation) if (a != None and b != None)])
        plt.xlim(0,len(insulation))
        plt.ylim(min2,max2)
        plt.title(str(motif_start) + ' reversed difference', fontdict={'fontsize':7})
        ax2.set_xticks([], [])
        ax2.set_xticks([int((motif_start-(start-450000))/2048) + 256 - 32])
        ax2.set_xticklabels(['m'], fontdict={'fontsize':6})
        #ax2.set_yticks([], [])
        f.tight_layout()
        plt.savefig("tads/" + name + "/" + name + "_" + str(motif_start) + "_rev_TADs_difference.png", dpi=800,bbox_inches='tight')
        plt.clf()

        of.write(row.format(str(motif_start) + " reversed", "tads/" + name + "/" + name + "_" + str(motif_start) + "_rev_TADs.png", "tads/" + name + "/" + name + "_" + str(motif_start) + "_rev_TADs_difference.png"))

    of.write("</table> \n")
of.close()
