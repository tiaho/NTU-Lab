#!/usr/bin/perl
#
# start_stop_codon_usage_cds_lengths.pl
#
# uses source/feature: curated coding_exon
#
# prints the start/stop codons and counts
# check if start coordinate > stop coordinate
# prints the codon usage counts for the coding sequences
# prints the CDS lengths
# checks for frameshifts
# prints the codon usage counts for introns

use strict; use warnings;
use File::Slurp;
use Data::Dumper;

# the CDS ID of the sequence that we want the codon counts for
die "Enter the CDS ID of the desired gene or enter 'none'\n" unless @ARGV == 1;
my $desiredID = $ARGV[0];

# hashes
my (%total_start_codons, %total_stop_codons, %total_codon_usage, %gene_codon_usage, %gene_codon_count, %gene_codon_entropy);

# the input files
my @fasta_files = ("xx01", "xx02", "xx03", "xx04", "xx05", "xx06");
my @gff_files = ("chromI_curated_coding_exon.gff2", "chromII_curated_coding_exon.gff2", "chromIII_curated_coding_exon.gff2", "chromIV_curated_coding_exon.gff2", "chromV_curated_coding_exon.gff2", "chromX_curated_coding_exon.gff2");

# runs the script for each of the 6 chromosomes of C elegans
for (my $i = 0; $i <= 5; $i++){
	my ($start_codons_ref, $stop_codons_ref, $genome_codon_usage_ref, $cds_length_ref) = extract_data($fasta_files[$i], $gff_files[$i]);
# 	for my $key (keys $start_codons_ref){
# 		$total_start_codons{$key} += $$start_codons_ref{$key};
# 	}
# 	for my $key (keys $stop_codons_ref){
# 		$total_stop_codons{$key} += $$stop_codons_ref{$key};
# 	}
	for my $key (keys $genome_codon_usage_ref){ # codon usage frequency for genome
		$total_codon_usage{$key} += $$genome_codon_usage_ref{$key};
	}
}

# adds 1 to each codon count in genome, pseudocounting
for my $key (keys %total_codon_usage){
	$total_codon_usage{$key}++;
}

# adds 1 to each codon count in genes
for my $cdsID (keys %gene_codon_usage){
	for my $codon (keys %{ $gene_codon_usage{$cdsID} }){
		if ($codon eq "chromosome"){next;}
		else{$gene_codon_usage{$cdsID}{$codon}++;}
	}
}

# sums up the number of codons in the genome
my $genome_sum = 0;
for my $key (keys %total_codon_usage){
	$genome_sum += $total_codon_usage{$key};
}

# # prints the genome codon usage counts
# for my $key (keys %total_codon_usage){
# 	print"$key, $total_codon_usage{$key}, genome\n";
# }

# converts the codon usage counts for the genome to frequencies
for my $key (keys %total_codon_usage){ # $key = codon
	$total_codon_usage{$key} = $total_codon_usage{$key} / $genome_sum;
# 	print"$key, $total_codon_usage{$key}\n";
}

# print Dumper \%gene_codon_usage;

# calculates the total number of codons for each gene 
for my $cdsID (keys %gene_codon_usage){ # $key = CDS ID
	# sums up the number of codons per gene
	for my $codon (keys %{ $gene_codon_usage{$cdsID} }){ # $key2 = codon
		if ($codon eq "chromosome"){
			next;
		}else{
			if (exists $gene_codon_count{$cdsID}{num_of_codons}) {
				$gene_codon_count{$cdsID}{num_of_codons} += $gene_codon_usage{$cdsID}{$codon};
			}else{
				$gene_codon_count{$cdsID}{num_of_codons} = $gene_codon_usage{$codon};
			}
		}
	}
}

# # prints the codon usage frequencies for a certain cds ID
if ($desiredID eq "none"){
	# does nothing
}elsif (exists $gene_codon_usage{$desiredID}){
	for my $codon (keys %{ $gene_codon_usage{$desiredID} }){
		print"$codon, $gene_codon_usage{$desiredID}{$codon}, $desiredID\n";
	}
}else{
	print"Enter a valid CDS ID\n";
}

# calculates the entropy of each codon for each gene 
for my $cdsID (keys %gene_codon_usage){
	# converts counts to frequency and calculates part of the relative entropy (all except for sum)
	for my $codon (keys %{ $gene_codon_usage{$cdsID} }){
		if ($codon eq "chromosome"){
			next;
		}else{
			$gene_codon_usage{$cdsID}{$codon} = $gene_codon_usage{$cdsID}{$codon} / $gene_codon_count{$cdsID}{num_of_codons};
			my $P = $gene_codon_usage{$cdsID}{$codon} or die "$gene_codon_usage{$cdsID}{$codon}\n";
			my $Q = $total_codon_usage{$codon} or die "$total_codon_usage{$codon}\n";
			my $entropy = $P * log($P/$Q)/log(2);
			$gene_codon_entropy{$cdsID}{$codon}{entropy} = $entropy;
		}
	}
	
	# sums up the entropies to get the relative entropy
	my $relative_entropy = 0;
	for my $codon (keys %{ $gene_codon_entropy{$cdsID} }){
		if ($codon eq "chromosome"){next;}
		else{$relative_entropy += $gene_codon_entropy{$cdsID}{$codon}{entropy};}
	}
	$gene_codon_entropy{$cdsID}{relative_entropy} = $relative_entropy;
}

# # prints the cdsID and relative entropy
# for my $key (keys %gene_codon_entropy){
# 	print"$gene_codon_entropy{$key}{relative_entropy}, $gene_codon_usage{$key}{chromosome}, $key\n";
# }

# print Dumper \%gene_codon_entropy;
	
# prints the start codons and counts
# print "---Start codons---\n";
# for my $key (keys %total_start_codons){
#     print"$key, $total_start_codons{$key}\n";
# }

# prints the stop codons and counts
# print "\n---Stop codons---\n";
# for my $key (keys %total_stop_codons){
#     print"$key, $total_stop_codons{$key}\n";
# }

# prints the codon usage counts
# print "\n---Codon usage---\n";
# for my $key (keys %total_codon_usage){
#     print"$key, $total_codon_usage{$key}\n";
# }

#################
###SUBROUTINES###
#################

# extracts data from a chromosome
sub extract_data{
    my ($fasta, $gff) = @_;
    
    # hashes
	my (%cds_seq, %start_codons, %stop_codons, %genome_codon_usage);
	my $cds_seq_ref = \%cds_seq;
    my $start_codons_ref = \%start_codons;
    my $stop_codons_ref = \%stop_codons;
    my $genome_codon_usage_ref = \%genome_codon_usage;
    
    # reads in the dna sequence
    my $sequence = read_file($fasta);
    $sequence =~ s/>[A-Z]+//; # removes the header specifying the chromosome
    $sequence =~ s/\n//g; # removes all new lines

    # extracts start/stop codon positions
    open(my $in, "<$gff") or die "Cannot open $gff\n";
    while (my $line = <$in>){
		chomp $line;
	
		# extracts the coding sequence
		extract_cds($sequence, $line, $cds_seq_ref);
    
	# checks start and stop coordinates
# 	check_start_stop_coord($start_pos, $stop_pos, $line);
    }
    close $in;

    for my $key (keys %cds_seq){

		# extracts the start and stop codons
# 		extract_start_stop_codons($key, $cds_seq_ref, $start_codons_ref, $stop_codons_ref);
		
		# gets a count of the codon usage
		check_codon_usage($sequence, $key, $cds_seq_ref, $genome_codon_usage_ref);
		
		# gets a list of CDS lengths
# 		get_CDS_length($key, $cds_seq_ref);
		
		# checks for frameshifts
# 		check_frameshift($key, $cds_seq_ref); # no frameshifts
    }
    
	return($start_codons_ref, $stop_codons_ref, $genome_codon_usage_ref);
}

# extracts the coding sequence
sub extract_cds{
    my ($sequence, $line, $cds_seq_ref) = @_;
    my ($chromosome, $source, $feature, $start_pos, $stop_pos, undef, $strand, undef, $attributes) = split("\t", $line);
    
    # only wants the coordinates from curated coding_exon
    if ($source eq "curated" and $feature eq "coding_exon"){
		my ($cdsID) = ($attributes =~ m/CDS \"([A-Za-z0-9_]+\.\d+[a-z]?\.?\d?)/);
		$$cds_seq_ref{$cdsID}{strand} = $strand;
		$$cds_seq_ref{$cdsID}{chromosome} = $chromosome;
		$$cds_seq_ref{$cdsID}{line} = $line;

		my $coding_exon = substr($sequence, $start_pos - 1, $stop_pos + 1 - $start_pos); # -1 b/c perl string starts at 0, + 1 to take account of last base
		
		# reverse strand
		if ($strand eq "-"){$coding_exon =~ tr/atcg/tagc/;} # reverse it later, or else sequence will be out of order
		$$cds_seq_ref{$cdsID}{cds} .= $coding_exon;
# 		print"$cdsID, $strand, $$cds_seq_ref{$cdsID}{cds}\n";
    }
}

# checks start and stop coordinates
sub check_start_stop_coord{
    my ($start_pos, $stop_pos, $line) = @_;
    if ($start_pos == $stop_pos){
# 	print"Start coordinate is equal to the stop coordinate $start_pos = $stop_pos\n";
	print"$line\n";
#     }elsif ($start_pos > $stop_pos){
# 	print"Start coordinate is greater than the stop coordinate $start_pos > $stop_pos\n";
# 	print"$line\n";
#     }else{
# 	print"$start_pos < $stop_pos\n";
    }
}

# extracts the start and stop codons
sub extract_start_stop_codons{
    my ($key, $cds_seq_ref, $start_codons_ref, $stop_codons_ref) = @_;
#     print"$key\n";
    my ($start, $stop);
	if ($$cds_seq_ref{$key}{strand} eq "-"){
		$$cds_seq_ref{$key}{cds} = reverse($$cds_seq_ref{$key}{cds}); # the cds from extract_cds() was not reversed for - strand
		$$cds_seq_ref{$key}{rev} = "yes";
	}else{
		$$cds_seq_ref{$key}{rev} = "no";
	}
	$start = substr($$cds_seq_ref{$key}{cds}, 0, 3);
	$stop = substr($$cds_seq_ref{$key}{cds}, - 3);
	
# 	print"$start, $stop\n";
	
	# gets a count of each codon
	if (exists $$start_codons_ref{$start}){$$start_codons_ref{$start}++;}
	else{$$start_codons_ref{$start} = 1;}
	if (exists $$stop_codons_ref{$stop}){$$stop_codons_ref{$stop}++;}
	else{$$stop_codons_ref{$stop} = 1;}
    
	# prints out the line if start codon isn't ATG and if stop codon isn't TAA/TAG/TGA
        if ($start ne "atg"){print"start, $start, $$cds_seq_ref{$key}{chromosome}, $$cds_seq_ref{$key}{strand}, $$cds_seq_ref{$key}{rev}, $key, $$cds_seq_ref{$key}{cds}\n";}
        if ($stop ne "taa" and $stop ne "tag" and $stop ne "tga"){print"stop, $stop, $$cds_seq_ref{$key}{chromosome}, $$cds_seq_ref{$key}{strand}, $$cds_seq_ref{$key}{rev}, $key, $$cds_seq_ref{$key}{cds}\n";}
}

# gets a count of the codon usage
sub check_codon_usage{
    my ($sequence, $key, $cds_seq_ref, $genome_codon_usage_ref) = @_;
    if ($$cds_seq_ref{$key}{strand} eq "-"){
		$$cds_seq_ref{$key}{cds} = reverse($$cds_seq_ref{$key}{cds});
    }
	my $CDS = $$cds_seq_ref{$key}{cds};
	my $seq_length = length($CDS);
	$gene_codon_usage{$key}{chromosome} = $$cds_seq_ref{$key}{chromosome};
# 	print"$gene_codon_usage{$key}{chromosome}\n";
	
	# for whole genome
    for (my $i = 1; $i <= $seq_length; $i += 3){
		my $codon = substr($CDS, $i - 1, 3);
		if (exists $$genome_codon_usage_ref{$codon}) {$$genome_codon_usage_ref{$codon}++;}
		else{$$genome_codon_usage_ref{$codon} = 1;}
    }
    
    # for each gene
	for (my $i = 1; $i <= $seq_length; $i += 3){
		my $codon = substr($CDS, $i - 1, 3);
		if (exists $gene_codon_usage{$key}{$codon}) {$gene_codon_usage{$key}{$codon}++;}
		else {$gene_codon_usage{$key}{$codon} = 1;}
	}
}

# gets a list of CDS lengths
sub get_CDS_length{
    my ($key, $cds_seq_ref) = @_;
    my $length = length($$cds_seq_ref{$key}{cds});
    print"$length, $$cds_seq_ref{$key}{chromosome}, $key\n";
}

# checks for frameshifts
sub check_frameshift{
    my ($key, $cds_seq_ref) = @_;
    my $modulo = length($$cds_seq_ref{$key}{cds}) % 3;
    if ($modulo != 0){print "$modulo\n";}
#     else{print"okay\n";}
}
