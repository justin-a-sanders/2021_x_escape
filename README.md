# 2021_x_escape # 
### Code for analyzing tad boundaries around escape genes with Akita ###

The following is a guide to the input/output of all of the scripts in the order you should use them. Unfortunately, because some of the steps in the pipeline are manual I don't have a single script to run everything. 

The first few steps are getting all of the CTCF motifs that fall within ChIP peaks on chromosome X. I've already done this for human GM12878 cells and mouse cells using ChIP profiles from ENCODE and provided the results in the files peak_motifs.bed and peak_motifs_mm10.bed, so unless you want to use a new set of ChIP peaks or analyze a new celltype/species, you can skip ahead to the prediction step. 

**get_peak_seqs.py** - Takes in two arguments, a genomic Fasta file with a sequence for 'chrX' and a bed file of CTCF peaks in the cell type of interest. It produces as output the file peak_seqs.fasta which has the sequence underlying every peak on the X-chromosome. 

The next step is to upload this file to a motif scanning tool such as FIMO. In addition, you will need to upload a .meme file for the motif you want to scan for. The motif for CTCF, ctcf.meme, is included in the repository. Download the result as a .tsv file. 

**motifs_to_bed.py** - Takes as input the .tsv file of motifs identified by FIMO. Produces the file peak_motifs.bed, which has the locations and strandedness of every CTCF motif within a ChIP peak on chrX. 

SKIP TO HERE IF YOU WANT TO USE CURRENT MOUSE OR HUMAN MOTIFS

The next step requires Akita. I reccomend following the install guide [here](https://github.com/calico/basenji) to setup the local enviornment and download the model. Then, the easiest way to use it is probably to move the file make_predictions_escape.py to the /manuscripts/akita/ directory in the basenji git repo, along with the file peak_motifs.bed produced by the previous step.  

**make_predictions_escape.py** - Takes two inputs, the path to a FASTA file containing a sequence for 'chrX', and a bed file of regions/genes that you want to make predictions of. The escape genes we have been using thus far are provided in the files mouse_escape_genes.bed and human_escape_genes.bed. Produces as output .npy objects in /predictions/<gene_name>/* containing the predicted Hi-C matrices for each perturbation of each motif.

**plot_pred_changes.py** - Takes two inputs, the bed file of escape genes and the bed file of peak motifs. Additionally, make sure the /predictions/* folder produced by the last step is in the same directory. This script plots both the predicted Hi-C matrices and the differences from the control when each motif was perturbed, as well as the correlations between each perturbed and unperturbed matrix. Additionally, it lays out all these results in the file motif_predictions.html. 

**call_TADs.py** - Takes two inputs, the bed file of escape genes and the bed file of peak motifs. Additionally, again make sure the /predictions/* folder produced by the last step is in the same directory. This script calls TADs on each Hi-C matrix and plots the differences in insulation scores. Additionally, it lays out all the plots it generates in the file tad_predictions.html. 
