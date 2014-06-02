#!/usr/bin/perl
#
# intron_lengths.pl
#
# calculates the lengths of the introns 
#
# uses these source/feature types: curated exon


use strict; use warnings;
use File::Slurp;
use Data::Dumper;

# hashes
my (%exon_coords, %intron_info);

# the input files
my @fasta_files = ("xx01", "xx02", "xx03", "xx04", "xx05", "xx06");
my @gff_files = ("chromI_curated_exon.gff2", "chromII_curated_exon.gff2", "chromIII_curated_exon.gff2", "chromIV_curated_exon.gff2", "chromV_curated_exon.gff2", "chromX_curated_exon.gff2");
my @sequences = my ($I, $II, $III, $IV, $V, $X);

# extracts curated exon coordinates and reads in the genomic seqence for each of the 6 chromosomes of C elegans
for (my $i = 0; $i <= 5; $i++){
	
	# extracts the coordinates of the exons for each gene
	extract_coords($gff_files[$i]);

	# reads in the fasta files into a string
	$sequences[$i] = read_sequence($fasta_files[$i]);
}
# print Dumper %exon_coords;

# gets the intron coordinates from the exon coordinates
for my $cdsID (keys %exon_coords){
	next if ($exon_coords{$cdsID}{count} == 1);
	for (my $i = 1; $i < $exon_coords{$cdsID}{count}; $i++){
		$intron_info{$cdsID}{$i}{start} = $exon_coords{$cdsID}{$i}{stop} + 1;
		$intron_info{$cdsID}{$i}{stop} = $exon_coords{$cdsID}{$i + 1}{start} - 1;
	}
	$intron_info{$cdsID}{strand} = $exon_coords{$cdsID}{strand};
	$intron_info{$cdsID}{chromosome} = $exon_coords{$cdsID}{chromosome};
	$intron_info{$cdsID}{count} = $exon_coords{$cdsID}{count} - 1;
}
# print Dumper %intron_info;

# calculates the lengths of the introns
for my $cdsID (keys %intron_info){
	for (my $i = 1; $i <= $intron_info{$cdsID}{count}; $i++){
		$intron_info{$cdsID}{lengths}{$i} = $intron_info{$cdsID}{$i}{stop} - $intron_info{$cdsID}{$i}{start} + 1;
	}
}
# print Dumper %intron_info;

# averages the lengths
for my $cdsID (keys %intron_info){
	my $sum = 0;
	for (my $i = 1; $i <= $intron_info{$cdsID}{count}; $i++){
		$sum += $intron_info{$cdsID}{lengths}{$i}; 
	}
	$intron_info{$cdsID}{lengths}{avg} = $sum / $intron_info{$cdsID}{count};
}
# print Dumper %intron_info;

# prints the cdsID, chromosome, and avg intron length
for my $cdsID (keys %intron_info) {
	my $length = log($intron_info{$cdsID}{lengths}{avg});
	print"$length, $intron_info{$cdsID}{chromosome}, $cdsID\n";
}

#################
###SUBROUTINES###
#################

# extracts the coordinates of exons for each gene
sub extract_coords{
	my $gff = $_[0];
    open(my $in, "<$gff") or die "Cannot open $gff\n";
    while (my $line = <$in>){
		chomp $line;
		my ($chromosome, $source, $feature, $start, $stop, undef, $strand, undef, $attributes) = split("\t", $line);
		my ($cdsID) = ($attributes =~ m/CDS \"([A-Za-z0-9_]+\.\d+[a-z]?\.?\d?)/);
		
		# assigns the info to a hash
		$exon_coords{$cdsID}{chromosome} = $chromosome;
		$exon_coords{$cdsID}{strand} = $strand;
		if (exists $exon_coords{$cdsID}{count}){$exon_coords{$cdsID}{count}++;}
		else {$exon_coords{$cdsID}{count} = 1;}
		my $count = $exon_coords{$cdsID}{count};
		$exon_coords{$cdsID}{$count}{start} = $start;
		$exon_coords{$cdsID}{$count}{stop} = $stop;
    }
    close $in;
}

# reads in the fasta files into a string
sub read_sequence{
	my $fasta = $_[0];
	
	# reads in the dna sequence
    my $sequence = read_file($fasta);
    $sequence =~ s/>[A-Z]+//; # removes the header specifying the chromosome
    $sequence =~ s/\n//g; # removes all new lines
    
    return($sequence);
}