import numpy as np
import csv
import sys

lengths = []

from pyfaidx import Fasta
# chroms = Fasta('hg38.ml.fa') #human
# chroms = Fasta('chrX.fa') #mouse
chroms = Fasta(sys.argv[1])
peaks = sys.argv[2]
with open("peak_seqs.fasta", 'w') as out:
    with open(peaks) as fd:
        rd = csv.reader(fd, delimiter="\t")
        for row in rd:
            if row[0] == 'chrX':
                row = row[0:6] + row[1:3] + ["100,0,0"]
                out.write(">" + "_".join(row) + "\n")
                out.write(str(chroms['chrX'][int(row[1]):int(row[2])]) + "\n")
                lengths.append(int(row[2])-int(row[1]))

import matplotlib.pyplot as plt
print(len(lengths))
plt.hist(lengths, density=False, bins=30)
plt.ylabel('Count')
plt.xlabel('Peak length');
plt.show()
