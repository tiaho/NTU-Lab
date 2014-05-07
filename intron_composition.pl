#!/usr/bin/perl
#
# intron_composition.pl
#
# calculates the composition of the introns
# uses source/feature: curated exon

use strict; use warnings;
use File::Slurp;
use Data::Dumper;

### get the exon coordinates, extract exons, leaving only introns?
### compare them by freq of A, C, G, T


# hashes
my (%exon_coords, %intron_coords, %genome_intron_info);

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
		$intron_coords{$cdsID}{$i}{start} = $exon_coords{$cdsID}{$i}{stop} + 1;
		$intron_coords{$cdsID}{$i}{stop} = $exon_coords{$cdsID}{$i + 1}{start} - 1;
	}
	$intron_coords{$cdsID}{strand} = $exon_coords{$cdsID}{strand};
	$intron_coords{$cdsID}{chromosome} = $exon_coords{$cdsID}{chromosome};
	$intron_coords{$cdsID}{count} = $exon_coords{$cdsID}{count} - 1;
}
# print Dumper %intron_coords;

# extracts the intron sequences
for my $cdsID (keys %intron_coords){
	my $chrom_seq;
	
	# assigns the correct genomic sequence for each chromosome
	if ($intron_coords{$cdsID}{chromosome} eq "I") {$chrom_seq = $sequences[0];}
	elsif ($intron_coords{$cdsID}{chromosome} eq "II") {$chrom_seq = $sequences[1];}
	elsif ($intron_coords{$cdsID}{chromosome} eq "III") {$chrom_seq = $sequences[2];}
	elsif ($intron_coords{$cdsID}{chromosome} eq "IV") {$chrom_seq = $sequences[3];}
	elsif ($intron_coords{$cdsID}{chromosome} eq "V") {$chrom_seq = $sequences[4];}
	elsif ($intron_coords{$cdsID}{chromosome} eq "X") {$chrom_seq = $sequences[5];}
	for (my $i = 1; $i <= $intron_coords{$cdsID}{count}; $i++){
		extract_intron_seq($chrom_seq, $cdsID, $intron_coords{$cdsID}{strand}, $intron_coords{$cdsID}{$i}{start}, $intron_coords{$cdsID}{$i}{stop});
	}
	
}
# print Dumper %intron_coords;

# calculates the counts for each nucleotide for the whole genome
for my $cdsID (keys %intron_coords){
	for (my $i = 0; $i < length($intron_coords{$cdsID}{seq}); $i++){
		my $nuc = substr($intron_coords{$cdsID}{seq}, $i, 1);
		
		# stores the counts in a hash for each individual gene 
		if (exists $intron_coords{$cdsID}{counts}{$nuc}) {$intron_coords{$cdsID}{counts}{$nuc}++;}
		else {$intron_coords{$cdsID}{counts}{$nuc} = 1;}
		
		# stores the counts in a hash to calculate the genome wide frequency of nucleotides in introns
		if (exists $genome_intron_info{counts}{$nuc}) {$genome_intron_info{counts}{$nuc}++;}
		else {$genome_intron_info{counts}{$nuc} = 1;}
	}
}
# print Dumper %intron_coords;
# print Dumper %genome_intron_info;

# converts counts into frequencies (percentages) for the individual gene
for my $cdsID (keys %intron_coords){
	my $a = $intron_coords{$cdsID}{counts}{a};
	my $c = $intron_coords{$cdsID}{counts}{c};
	my $g = $intron_coords{$cdsID}{counts}{g};
	my $t = $intron_coords{$cdsID}{counts}{t};
	my $total = $a + $c + $g + $t;

	$intron_coords{$cdsID}{freq}{a} = $a/$total;
	$intron_coords{$cdsID}{freq}{c} = $c/$total;
	$intron_coords{$cdsID}{freq}{g} = $g/$total;
	$intron_coords{$cdsID}{freq}{t} = $t/$total;
}
# print Dumper %intron_coords;

# converts counts into frequencies (percentages) for the genome
my $a = $genome_intron_info{counts}{a};
my $c = $genome_intron_info{counts}{c};
my $g = $genome_intron_info{counts}{g};
my $t = $genome_intron_info{counts}{t};
my $total = $a + $c + $g + $t;

$genome_intron_info{freq}{a} = $a/$total;
$genome_intron_info{freq}{c} = $c/$total;
$genome_intron_info{freq}{g} = $g/$total;
$genome_intron_info{freq}{t} = $t/$total;

# print Dumper %genome_intron_info;

# finds the KL distance for the frequency of each nucleotide in a gene to the frequency it is found in all introns genome wide
my @nucleotides = ("a", "c", "g", "t");
for my $cdsID (keys %intron_coords){
	for my $nuc (@nucleotides){
		my $P = $intron_coords{$cdsID}{freq}{$nuc};
		my $Q = $genome_intron_info{freq}{$nuc};
		my $entropy = $P * log($P/$Q)/log(2);
		$intron_coords{$cdsID}{entropy}{$nuc} = $entropy;
	}
}
# print Dumper %intron_coords;

# sums the entropies for each gene, to get the relative entropy
for my $cdsID (keys %intron_coords){
	$intron_coords{$cdsID}{entropy}{total} = 0;
	for my $nuc (@nucleotides){
		$intron_coords{$cdsID}{entropy}{total} += $intron_coords{$cdsID}{entropy}{$nuc};
	}
}
# print Dumper %intron_coords;

# prints the cdsID, chromosome, and relative entropy
for my $cdsID (keys %intron_coords) {
	print"$intron_coords{$cdsID}{entropy}{total}, $intron_coords{$cdsID}{chromosome}, $cdsID\n";
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

# extracts the splice sites
sub extract_intron_seq{
	my ($sequence, $cdsID, $strand, $start, $stop) = @_;
	
	my $intron = substr($sequence, $start - 1, $stop - $start + 1); 
	if ($strand eq "+"){
		if (exists $intron_coords{$cdsID}{seq}) {$intron_coords{$cdsID}{seq} .= $intron;}
		else {$intron_coords{$cdsID}{seq} = $intron;}
	}
	if ($strand eq "-"){
		$intron =~ tr/acgt/tgca/;
		$intron = reverse($intron);
		if (exists $intron_coords{$cdsID}{seq}) {$intron_coords{$cdsID}{seq} .= $intron;}
		else {$intron_coords{$cdsID}{seq} = $intron;}
	}
}