# Palicourea
Code for phylogenomics, climatic niche overlap, biogeography and character ste reconstruction for Bedoya et al., in prep.

Initial Analyses of target enrichment data (including orthology inference) follows Morales-Briones (Morales-Briones, D.F., B. Gehrke, H. Chien-Hsun Huang, A. Liston, M. Hong. H.E. Marx, D.C. Tank & Y. Yang. 2022. Analysis of paralogs in target enrichment data pinpoints multiple ancient polyploidy events in Alchemilla s.l. (Rosaceae). Systematic Biology 71(1):190–207) in the treatment of putative paralogs. Diego's scripts are indicated with **

## Filter loci: remove loci with low coverage across samples and total sequence length

After contig assembly with Hybpiper2, select loci with at least 50% sequence recovered in >20% sequences. Use [Loci_filtered.R](https://github.com/ambed0ya/Palicourea/blob/main/Loci_filtered.R "Loci_filtered.R script") (adapted from https://github.com/ajhelmstetter/afrodyn/blob/master/local_scripts/75_75.R) to identify loci with <50% sequence length and >80% missing data. Use the Rscript [identify_putative_paralogs.R](https://github.com/ambed0ya/Palicourea/blob/main/identify_putative_paralogs.R "odentify_putative_paralogs.R script") to identify the targeted loci with >1 paralog warning from paralog_report.tsv (putative paralogs). Results are in 'putative_paralogs_list.csv'

Move files identified with [Loci_filtered.R](https://github.com/ambed0ya/Palicourea/blob/main/Loci_filtered.R "Loci_filtered.R script") and not identified as putative paralogs to 'filtered' folder (created inside Hybpiper output folder).

Similarly, from 'paralogs_no_chimeras', move files identified as putative paralogs and with at least 50% sequence recovered in >20 sequences to 'filtered' folder (create folder first inside of 'paralogs_no_chimeras')

## Generate and edit alignments
Create scripts for alignment with macse in 'filtered' folder inside Hybpiper2 folder (single copy loci) and in 'filtered' folder inside 'paralogs_no_chimeras' folder (putative paralogs). Output sequence files for putative loci from Spades first need to be edited (e.g., >HIL_parasitica_Jiminez2189 single_hit to >HIL_parasitica_Jiminez2189). Use for example, sed (e.g., `sed -i '/ single_hit//}’ *.fasta`).

`for filename in $(ls *.fasta); do echo macse -prog alignSequences -seq $filename -out_NT $(cut -d'.' -f1 <<<"$filename").NT.aln -out_AA $(cut -d'.' -f1 <<<"$filename").AA.aln > $(cut -d'.' -f1 <<<"$filename").sh; done`

Run in parallel
`parallel -j 32 bash ::: *.sh`

Edit alignments
`python remove_shifted_codons_from_macse.py . aln . nt` **

`python pxclsq_wrapper.py . 0.1 dna` **

Infer trees
`for filename in $(ls *NT.fs.aln-cln); do echo raxml-ng --all --msa $filename --model GTR+G --prefix out$filename --seed 2 --threads 2 --bs-metric fbp,tbe > $filename.sh; done`
Run in parallel (`parallel -j 20 bash ::: *.sh`) 