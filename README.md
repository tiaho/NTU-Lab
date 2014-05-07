README.md
#Description of the files

###cds_codon_composition.pl
* calculates the composition of codons for all the coding sequences

###cds_codon_composition_distributions.R
* plots the distribution of the codon composition for all the coding sequences

###cds_codon_composition_relative_entropies.R
* plots the distribution of the KL distances for the codon composition in all genes

###cds_length_distributions.R
* plots the distribution of the CDS lengths

###cds_lengths.pl
* calculates the lengths of the coding sequences

###check_start_stop_coords.pl
* checks if the start coordinate comes after the stop coordinate

###compare_gene_rankings.R
* calculates a z score for each gene based on various factors (does same thing as rank_genes_zscores.pl)

###compare_support_for_isoforms.pl
* wrapper script for coordinate_extracter.pl and find_overlapping_features2.pl

###exon_intron_length_distribution.R
* plots the distribution of exon and intron lengths

###exon_lengths.pl
* calculates the lengths of exons

###exon_intron_lengths_old.pl
* calculates the lengths of exons and introns - old version

###exon_length_distributions.R
* plots the distribution of exon lengths

###frameshifts.pl
* checks for frameshifts

###genes_zscore_rank_distribution.R
* plots the distribution of the ranking of the genes (z-score method)

###genetic_algorithm_example.pl
* example of a genetic algorithm using shapes and their areas 

###intron_composition.pl
* calculates the composition of the introns

###intron_length_distributions.R
* plots the distribution of intron lengths
  
###intron_lengths.pl
* calculates the lengths of the introns 

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

###splice_site_composition.pl
* calculates the composition of the splice sites

###split_gff_file.pl
* splits a gff file by chromosome (currently hard coded for *C. elegans*). need to manually change file name each time!

###start_stop_codon_distributions.R
* plots the distribution of the start and stop codon usage

###start_stop_codon_composition.pl
* calculates the composition of start and stop codons

###start_stop_codon_usage_cds_length.pl
* does too many things...
    * prints the start/stop codons and counts
    * check if start coordinate > stop coordinate
    * prints the codon usage counts for all the coding sequences
    * prints the CDS lengths
    * checks for frameshifts