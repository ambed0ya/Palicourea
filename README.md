# Palicourea
Code for phylogenomics, climatic niche overlap, biogeography and character state reconstruction for Bedoya et al., In Review.

## Contents







Initial Analyses of target enrichment data (including orthology inference) follows Morales-Briones (Morales-Briones, D.F., B. Gehrke, H. Chien-Hsun Huang, A. Liston, M. Hong. H.E. Marx, D.C. Tank & Y. Yang. 2022. Analysis of paralogs in target enrichment data pinpoints multiple ancient polyploidy events in Alchemilla s.l. (Rosaceae). Systematic Biology 71(1):190–207) in the treatment of putative paralogs. Diego's scripts are indicated with **

## Filter loci: remove loci with low coverage across samples and total sequence length

After contig assembly with Hybpiper2, select loci with at least 50% sequence recovered in >20% sequences. Use [Loci_filtered.R](https://github.com/ambed0ya/Palicourea/blob/main/Loci_filtered.R "Loci_filtered.R script") (adapted from https://github.com/ajhelmstetter/afrodyn/blob/master/local_scripts/75_75.R) to identify loci with <50% sequence length and >80% missing data. Use the Rscript [identify_putative_paralogs.R](https://github.com/ambed0ya/Palicourea/blob/main/identify_putative_paralogs.R "identify_putative_paralogs.R script") to identify the targeted loci with >1 paralog warning from paralog_report.tsv (putative paralogs). Results from the latter are in 'putative_paralogs_list.csv'

Move files identified with [Loci_filtered.R](https://github.com/ambed0ya/Palicourea/blob/main/Loci_filtered.R "Loci_filtered.R script") and not identified as putative paralogs to 'filtered' folder (created inside Hybpiper output folder).

Similarly, from 'paralogs_no_chimeras', move files identified as putative paralogs and with at least 50% sequence recovered in >20 sequences to 'filtered' folder (create folder first inside of 'paralogs_no_chimeras')

## Generate and edit alignments
Create scripts for alignment with macse in 'filtered' folder inside Hybpiper2 folder (single copy loci) and in 'filtered' folder inside 'paralogs_no_chimeras' folder (putative paralogs). Output sequence files for putative paralogs from Spades first need to be edited (e.g., >HIL_parasitica_Jiminez2189 single_hit to >HIL_parasitica_Jiminez2189). Use for example, sed (e.g., `sed -i 's/ single_hit//g' *.fasta`).

`for filename in $(ls *.fasta); do echo macse -prog alignSequences -seq $filename -out_NT $(cut -d'.' -f1 <<<"$filename").NT.aln -out_AA $(cut -d'.' -f1 <<<"$filename").AA.aln > $(cut -d'.' -f1 <<<"$filename").sh; done`

Run in parallel
`parallel -j 32 bash ::: *.sh`

Edit alignments
`python remove_shifted_codons_from_macse.py . aln . nt` **

`python pxclsq_wrapper.py . 0.1 dna` **

Infer trees
`for filename in $(ls *NT.fs.aln-cln); do echo raxml-ng --all --msa $filename --model GTR+G --prefix out$filename --seed 2 --threads 2 --bs-metric fbp,tbe > $filename.sh; done`

Run in parallel (`parallel -j 20 bash ::: *.sh`)

##Processing putative paralogs
Change names of resulting FBP trees

`rename "s/.supportFBP/.tt/" *FBP`
`rename "s/out//" *.tt`
`rename "s/NT.fs.aln-cln.raxml.//" *.tt`

Create folder 'trees' and 'alignments' inside filtered paralogs_no_chimeras folder and move trees and alignments in there respectively. Go to each of those folders and use [name_changes.sh](https://github.com/ambed0ya/Palicourea/blob/main/name_changes.sh "name_changes.sh script") to change name structure in both alignments and trees for subsequent masking. (You'll need to edit the file to your species names; names@author but this may vary according to your dataset). Note: If preliminary analyses suggest that there are paraphyletic spp, you may need your name to include the collection number followed by @anything_you_want

Mask monophyletic and paraphyletic taxa:
`python mask_tips_by_taxonID_transcripts.py trees_folder alignments_folder y` (n if you don't want to mask paraphyletic taxa) **

Use treeshrink to detect and remove long branches that may be spurious branches
`python tree_shrink_wrapper.py . .tt.mm 0.01 treeshrink` (test different quantiles to see which one fit better your data; these are you final homolog trees) **

Rename all outputs in every folder with [files_folders_namechanges.sh](https://github.com/ambed0ya/Palicourea/blob/main/files_folders_name_changes.sh "files_folders_name_changes.sh script") so they can then be all moved to a new folder (MO) to be pruned (e.g. `mv -t MO /.mm`)

Run MO keeping clades with >20% total taxa using all outgroups
`python prune_paralogs_MO.py your_MO_folder .mm 23 out_20percent` ** You'll need to modify this script to specify your own ingroups and outgroups. Also, play around with the minimum number of taxa to indicate in the MO analyses (23 here). For this, create a matrix occupancy stats plot as follows:
`python ortholog_occupancy_stats.py out_20percent/` **

Create fasta files from the final tree files to re-run trees and use those input trees for downstream analyses (optional)
`python write_ortholog_fasta_from_multiple_aln.py alignments_folder out_20percent aln-cln .tre output_folder` **

Align and infer trees as described above

Collapse low-support branches with newick utilities (see example below)
`nw_ed  810-allexons.trees 'i & b<=10' o > 810-allexons-BS10.trees`
