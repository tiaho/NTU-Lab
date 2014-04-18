#!/usr/bin/perl
#
# start_stop_coordinates.pl
#
# checks if the start coordinate comes after the stop coordinate
# uses source/feature: curated coding_exon


use strict; use warnings;
use File::Slurp;

# the input files
my @fasta_files = ("xx01", "xx02", "xx03", "xx04", "xx05", "xx06");
my @gff_files = ("chromI_curated_coding_exon.gff2", "chromII_curated_coding_exon.gff2", "chromIII_curated_coding_exon.gff2", "chromIV_curated_coding_exon.gff2", "chromV_curated_coding_exon.gff2", "chromX_curated_coding_exon.gff2");

# runs the script for each of the 6 chromosomes of C elegans
for (my $i = 0; $i <= 5; $i++){
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
		my ($chromosome, $source, $feature, $start_pos, $stop_pos, undef, $strand, undef, $attributes) = split("\t", $line);
    
		# checks start and stop coordinates
		check_start_stop_coord($start_pos, $stop_pos, $line);
    }
    close $in;
}

# checks start and stop coordinates
sub check_start_stop_coord{
    my ($start_pos, $stop_pos, $line) = @_;
    if ($start_pos == $stop_pos){
		print"Start coordinate is equal to the stop coordinate $start_pos = $stop_pos\n";
		print"$line\n";
    }elsif ($start_pos > $stop_pos){
		print"Start coordinate is greater than the stop coordinate $start_pos > $stop_pos\n";
		print"$line\n";
    }
}