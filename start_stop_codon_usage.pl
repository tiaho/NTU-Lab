#!/usr/bin/perl
#
# start_stop_codon_usage.pl
#
# calculates the usage of start and stop codons
# uses source/feature: curated coding_exon

use strict; use warnings;
use File::Slurp;

# wants to print start or stop codon counts?
die "Usage: start_stop_codon_usage.pl <start OR stop>\n" unless @ARGV == 1 and ($ARGV[0] eq "start" or $ARGV[0] eq "stop");
my $type = $ARGV[0];

# hashes
my (%total_start_codons, %total_stop_codons);

# the input files
my @fasta_files = ("xx01", "xx02", "xx03", "xx04", "xx05", "xx06");
my @gff_files = ("chromI_curated_coding_exon.gff2", "chromII_curated_coding_exon.gff2", "chromIII_curated_coding_exon.gff2", "chromIV_curated_coding_exon.gff2", "chromV_curated_coding_exon.gff2", "chromX_curated_coding_exon.gff2");

# runs the script for each of the 6 chromosomes of C elegans
for (my $i = 0; $i <= 1; $i++){
	my ($start_codons_ref, $stop_codons_ref) = extract_data($fasta_files[$i], $gff_files[$i]);
	for my $key (keys $start_codons_ref){
		$total_start_codons{$key} += $$start_codons_ref{$key};
	}
	for my $key (keys $stop_codons_ref){
		$total_stop_codons{$key} += $$stop_codons_ref{$key};
	}
}

# prints the start codons and counts
if ($type eq "start"){
	for my $key (keys %total_start_codons){
		print"$key, $total_start_codons{$key}\n";
	}
}

# prints the stop codons and counts
if ($type eq "stop"){
	for my $key (keys %total_stop_codons){
		print"$key, $total_stop_codons{$key}\n";
	}
}


#################
###SUBROUTINES###
#################

# extracts data from a chromosome
sub extract_data{
    my ($fasta, $gff) = @_;
    
    # hashes
	my (%cds_seq, %start_codons, %stop_codons);
	my $cds_seq_ref = \%cds_seq;
    my $start_codons_ref = \%start_codons;
    my $stop_codons_ref = \%stop_codons;
    
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
# 		check_start_stop_coord($start_pos, $stop_pos, $line);
    }
    close $in;

    for my $key (keys %cds_seq){
		# extracts the start and stop codons
		extract_start_stop_codons($key, $cds_seq_ref, $start_codons_ref, $stop_codons_ref);
    }
    
	return($start_codons_ref, $stop_codons_ref);
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