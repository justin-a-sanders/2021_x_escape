import numpy as np
import csv
import sys

all_peaks = set()
num_motifs = 0
motif_in = sys.argv[1]

b = True
with open("peak_motifs.bed", 'w') as out:
    with open(motif_in) as fd:
        rd = csv.reader(fd, delimiter="\t")
        for row in rd:
            if b:
                b = False
                continue
            if float(row[8]) <= 0.01:
                all_peaks.add(row[2])
                motif = row[2].split("_")
                motif[2] = str(int(motif[1]) + int(row[4]))
                motif[1] = str(int(motif[1]) + int(row[3]))
                motif[5] = row[5]
                if motif[5] == '+':
                    motif[8] = '0,0,200'
                else:
                    motif[8] = '200,0,0'
                # motif.append(row[9])
                out.write("\t".join(motif) + "\n")
                num_motifs += 1

print("Number of motifs found: ", num_motifs)
print("Distinct peaks with motifs: ", len(all_peaks))
print("Peaks without any motifs: ", 1469 - len(all_peaks))
