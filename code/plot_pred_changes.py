import os
import csv
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.decomposition import PCA

escape_genes = sys.argv[1]
peak_motifs = sys.argv[2]

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
        for name, start, end in genes:
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
<td><h4>{8}</h4></td>
<td><img src='{0}' height=240 width=320></td>
<td><img src='{1}' height=240 width=320></td>
<td><img src='{2}' height=240 width=320></td>
<td><img src='{3}' height=240 width=320></td>
</tr>

<tr>
<td></td>
<td><p>Pearson: {4:.3f} , Spearman: {5:.3f}</p></td>
<td></td>
<td><p>Pearson: {6:.3f} , Spearman: {7:.3f}</p></td>
<td></td>
</tr>"""

f = open("motif_predictions.html", "w")
f.write(prelude)
for name, start, end in genes:
    print(name)
    os.makedirs("plots/" + name, exist_ok=True)
    unchanged = np.load("predictions/" + name + "/" + name + "_unchanged.npy")
    fig, ax = plt.subplots()
    im = plt.matshow(unchanged, fignum=False, cmap= 'RdBu_r', vmax=vmax, vmin=vmin)
    plt.colorbar(im, fraction=.04, pad = 0.05);
    plt.title(name + ' unchanged')
    plt.ylabel('chrX: ' + str(start-450000) + "-" + str(start-450000+2**20))
    gene_start_bin = int(450000/2048) + 256 - 32
    gene_end_bin = int((450000 + end - start)/2048) + 256 - 32
    plt.xticks(ticks=[gene_start_bin, gene_end_bin], labels=["s", "e"])
    fig.tight_layout()
    plt.savefig("plots/" + name + "/" + name + "_unchanged.png")
    plt.show()
    plt.clf()

    header = "\n <h2> " + name + """ </h2>
        <table> <tr>
        <td><h3>Motif</h3></td>
        <td><h3>Deletion</h3></td>
        <td><h3>Deletion Fold Change</h3></td>
        <td><h3>Inversion</h3></td>
        <td><h3>Inversion Fold Change</h3></td>
        </tr> \n"""
    f.write(header)

    for motif_start, motif_end in gene_to_motifs[name]:
        del_name = "plots/" + name + "/" + name + "_" + str(motif_start) + "_deleted.png"
        del_diff_name = "plots/" + name + "/" + name + "_" + str(motif_start) + "_deleted_log_diff.png"
        inv_name = "plots/" + name + "/" + name + "_" + str(motif_start) + "_reversed.png"
        inv_diff_name = "plots/" + name + "/" + name + "_" + str(motif_start) + "_reversed_log_diff.png"
        # plot predictions and log fold change from control for motif inversion
        reversed = np.load("predictions/" + name + "/" + name + "_" + str(motif_start) + "_reversed.npy")
        im = plt.matshow(reversed, fignum=False, cmap= 'RdBu_r', vmax=vmax, vmin=vmin)
        plt.colorbar(im, fraction=.04, pad = 0.05);
        plt.title(name + " " + str(motif_start) + ' reversed compliment')
        plt.ylabel('chrX: ' + str(start-450000) + "-" + str(start-450000+2**20))
        gene_start_bin = int(450000/2048) + 256 - 32
        gene_end_bin = int((450000 + end - start)/2048) + 256 - 32
        motif_bin = int((motif_start-(start-450000))/2048) + 256 - 32
        plt.xticks(ticks=[gene_start_bin, gene_end_bin, motif_bin], labels=["s", "e", "m"])
        fig.tight_layout()
        plt.savefig(inv_name)
        plt.clf()

        im = plt.matshow(reversed - unchanged, fignum=False, cmap= 'RdBu_r', vmax=vmax, vmin=vmin)
        plt.colorbar(im, fraction=.04, pad = 0.05);
        plt.title(name + " " + str(motif_start) + ' reversed compliment log fold change')
        plt.ylabel('chrX: ' + str(start-450000) + "-" + str(start-450000+2**20))
        plt.xticks(ticks=[gene_start_bin, gene_end_bin, motif_bin], labels=["s", "e", "m"])
        fig.tight_layout()
        plt.savefig(inv_diff_name)
        plt.clf()

        # do the same for motif deletion
        deleted = np.load("predictions/" + name + "/" + name + "_" + str(motif_start) + "_deletion.npy")
        im = plt.matshow(deleted, fignum=False, cmap= 'RdBu_r', vmax=vmax, vmin=vmin)
        plt.colorbar(im, fraction=.04, pad = 0.05);
        plt.title(name + " " + str(motif_start) + ' deleted')
        plt.ylabel('chrX: ' + str(start-450000) + "-" + str(start-450000+2**20))
        plt.xticks(ticks=[gene_start_bin, gene_end_bin, motif_bin], labels=["s", "e", "m"])
        fig.tight_layout()
        plt.savefig(del_name)
        plt.clf()

        im = plt.matshow(deleted - unchanged, fignum=False, cmap= 'RdBu_r', vmax=vmax, vmin=vmin)
        plt.colorbar(im, fraction=.04, pad = 0.05);
        plt.title(name + " " + str(motif_start) + ' deletion log fold change')
        plt.ylabel('chrX: ' + str(start-450000) + "-" + str(start-450000+2**20))
        plt.xticks(ticks=[gene_start_bin, gene_end_bin, motif_bin], labels=["s", "e", "m"])
        fig.tight_layout()
        plt.savefig(del_diff_name)
        plt.clf()

        unchanged = np.nan_to_num(unchanged)
        reversed = np.nan_to_num(reversed)
        deleted = np.nan_to_num(deleted)

        rev_correlation = np.corrcoef(np.ndarray.flatten(unchanged), np.ndarray.flatten(reversed))[0,1]
        del_correlation = np.corrcoef(np.ndarray.flatten(unchanged), np.ndarray.flatten(deleted))[0,1]
        rev_spearman = stats.spearmanr(np.ndarray.flatten(unchanged), np.ndarray.flatten(reversed))[0]
        del_spearman = stats.spearmanr(np.ndarray.flatten(unchanged), np.ndarray.flatten(deleted))[0]

        pca = PCA(n_components=1)
        unchanged_pc1 = pca.fit_transform(unchanged)
        rev_pc1 = pca.fit_transform(reversed)
        del_pc1 = pca.fit_transform(deleted)

        rev_pc_corr = np.corrcoef(np.ndarray.flatten(unchanged_pc1), np.ndarray.flatten(rev_pc1))[0,1]
        del_pc_corr = np.corrcoef(np.ndarray.flatten(unchanged_pc1), np.ndarray.flatten(del_pc1))[0,1]

        f.write(row.format(del_name, del_diff_name, inv_name, inv_diff_name, del_correlation, del_spearman, rev_correlation, rev_spearman, motif_start))
        print("   ", name, motif_start, rev_correlation, rev_spearman, del_correlation, del_spearman)
    f.write("</table> \n")
f.close()
