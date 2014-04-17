README.md

#Description of the files

###cds_length_distributions.R
* plots the distribution of the CDS lengths

###codon_usage_distributions.R
* plots the distribution of the codon usage for all the coding sequences

###codon_usage_relative_entropies.R
* plots the distribution of the KL distances for the codon usages in all genes

###compare_gene_rankings.R
* calculates a z score for each gene based on various factors (does same thing as rank_genes_zscores.pl)

###compare_support_for_isoforms.pl
* wrapper script for coordinate_extracter.pl and find_overlapping_features2.pl

###coordinate_extracter.pl
* extracts the lines in the gff file that are within the range of a gene's coordinates

###exon_intron_length_distribution.R
* plots the distribution of exon and intron lengths

###exon_intron_lengths_old.pl
* calculates the lengths of exons and introns - old version

###exon_intron_lengths.pl
* calculates the lengths of exons and introns - current working version 

###exon_length_distributions.R
* plots the distribution of exon lengths

###genes_zscore_rank_distribution.R
* plots the distribution of the ranking of the genes (z-score method)

###genetic_algorithm_example.pl
* example of a genetic algorithm using shapes and their areas 

###intron_length_distributions.R
* plots the distribution of intron lengths

###make_source_feature_file_multiple.pl
* makes a file containing only the multiple desired source/feature combinations 

###make_source_feature_file.pl
* makes a file that contains only the single desired source/feature combination

###prelim_ranking_scores.R
* plots the distribution of the preliminary ranking scores (scored by # of overlapping bases)

###rank_genes_zscore.pl
* calculates a z score for each gene based on various factors (does same thing as compare_gene_rankings.R)

###sanity_check_OLD_SCRIPT.pl
* old script, refer to start_stop_codon_usage_cds_length.pl

###split_gff_file.pl
* splits a gff file by chromosome (currently hard coded for *C. elegans*)

###start_stop_codon_distributions.R
* plots the distribution of the start and stop codon usage

###start_stop_codon_usage_cds_length.pl
* does too many things...
    * prints the start/stop codons and counts
    * check if start coordinate > stop coordinate
    * prints the codon usage counts for all the coding sequences
    * prints the CDS lengths
    * checks for frameshifts