#!/usr/bin/perl
#
# frameshifts.pl
#
# checks for frameshifts
# uses source/feature: curated coding_exon

use strict; use warnings;
use File::Slurp;

# the input files
my @fasta_files = ("xx01", "xx02", "xx03", "xx04", "xx05", "xx06");
my @gff_files = ("chromI_curated_coding_exon.gff2", "chromII_curated_coding_exon.gff2", "chromIII_curated_coding_exon.gff2", "chromIV_curated_coding_exon.gff2", "chromV_curated_coding_exon.gff2", "chromX_curated_coding_exon.gff2");

# runs the script for each of the 6 chromosomes of C elegans
for (my $i = 0; $i <= 0; $i++){
	extract_data($fasta_files[$i], $gff_files[$i]);
}


#################
###SUBROUTINES###
#################

# extracts data from a chromosome
sub extract_data{
    my ($fasta, $gff) = @_;
    
    # hashes
	my %cds_seq;
	my $cds_seq_ref = \%cds_seq;
    
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
    }
    close $in;

    for my $key (keys %cds_seq){
		# checks for frameshifts
		check_frameshift($key, $cds_seq_ref); # no frameshifts
    }
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

# checks for frameshifts
sub check_frameshift{
    my ($key, $cds_seq_ref) = @_;
    my $modulo = length($$cds_seq_ref{$key}{cds}) % 3;
    if ($modulo != 0){print "$modulo\n";}
}