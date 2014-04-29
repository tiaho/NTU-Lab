#!/usr/bin/perl
#
# splice_site_composition.pl
#
# calculates the composition of the splice sites
# uses source/feature: gene gene

use strict; use warnings;
use File::Slurp;
use Data::Dumper;

# keep all the curated exon start/stop coordinate in a hash, with the key = cdID and keep a counter variable for num of exons

# store each chromosome info in a different hash?
# extract the fastfa file for each choromosome into a different variable

# hashes
my (%exon_coords, %genome_splice_sites, %splice_sites);

# the input files
my @fasta_files = ("xx01", "xx02", "xx03", "xx04", "xx05", "xx06");
my @gff_files = ("chromI_curated_exon.gff2", "chromII_curated_exon.gff2", "chromIII_curated_exon.gff2", "chromIV_curated_exon.gff2", "chromV_curated_exon.gff2", "chromX_curated_exon.gff2");
my @sequences = my ($I, $II, $III, $IV, $V, $X);

# extracts curated exon coordinates and reads in the genomic seqence for each of the 6 chromosomes of C elegans
for (my $i = 0; $i <= 0; $i++){
	
	# extracts the coordinates of the exons for each gene
	extract_coords($gff_files[$i]);

	# reads in the fasta files into a string
	$sequences[$i] = read_sequence($fasta_files[$i]);
}
# print Dumper %exon_coords;

for my $cdsID (keys %exon_coords){
	my $genome_seq;

	# assigns the correct genomic sequence for each chromosome
	if ($exon_coords{$cdsID}{chromosome} eq "I") {$genome_seq = $sequences[0];}
	elsif ($exon_coords{$cdsID}{chromosome} eq "II") {$genome_seq = $sequences[1];}
	elsif ($exon_coords{$cdsID}{chromosome} eq "III") {$genome_seq = $sequences[2];}
	elsif ($exon_coords{$cdsID}{chromosome} eq "IV") {$genome_seq = $sequences[3];}
	elsif ($exon_coords{$cdsID}{chromosome} eq "V") {$genome_seq = $sequences[4];}
	elsif ($exon_coords{$cdsID}{chromosome} eq "X") {$genome_seq = $sequences[5];}

	for (my $i = 1; $i <= $exon_coords{$cdsID}{count}; $i++){
# 		extract_cds($genome_seq, $cdsID, $exon_coords{$cdsID}{strand}, $exon_coords{$cdsID}{$i}{start}, $exon_coords{$cdsID}{$i}{stop});

		# extracts splice sites
		my $position;
		if ($i == 1) {$position = "beg";}
		elsif ($i == $exon_coords{$cdsID}{count}) {$position = "end";}
		else {$position = "middle";}
		extract_splice_sites($position, $genome_seq, $cdsID, $exon_coords{$cdsID}{strand}, $exon_coords{$cdsID}{$i}{start}, $exon_coords{$cdsID}{$i}{stop});
	}
	
}

# reverses the sequences on the opposite strand
# for my $key (keys %exon_coords){
# 	if ($exon_coords{$key}{strand} eq "-"){
# 		$exon_coords{$key}{seq} =~ tr/atcg/tagc/;
# 		$exon_coords{$key}{seq} = reverse($exon_coords{$key}{seq});
# 	}
# }

print Dumper %genome_splice_sites;

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

# extracts the coding sequence
sub extract_cds{
    my ($sequence, $cdsID, $strand, $start, $stop) = @_;

	# extracts the sequence of the gene
	my $cds_seq = substr($sequence, $start - 1, $stop + 1 - $start); # -1 b/c perl string starts at 0, + 1 to take account of last base
	
	# reverse strand -- will reverse later, so order of the exons don't get messed up
# 	if ($strand eq "-"){
# 		$cds_seq =~ tr/atcg/tagc/;
# 		$cds_seq = reverse($cds_seq);
# 	}
	
	if (exists $exon_coords{$cdsID}{seq}) {$exon_coords{$cdsID}{seq} .= $cds_seq;}
	else {$exon_coords{$cdsID}{seq} = $cds_seq;}
}

# extracts the splice sites            ONLY LOOK AT POSITIVE STRANDS FIRST
sub extract_splice_sites{
	my ($position, $sequence, $cdsID, $strand, $start, $stop) = @_;
	
	my ($splice1, $splice2);
	
	if ($position eq "beg"){
		# extract 2 bases after the stop codon
		$splice1 = substr($sequence, $stop, 2); # wants to start at the base after the stop codon
		if ($strand eq "+"){
			if (exists $genome_splice_sites{donor}{$splice1}) {$genome_splice_sites{donor}{$splice1}++;}
			else {$genome_splice_sites{donor}{$splice1} = 1;}
		}
		if ($strand eq "-"){
			$splice1 =~ tr/acgt/tgca/;
			$splice1 = reverse($splice1);
			if (exists $genome_splice_sites{acceptor}{$splice1}) {$genome_splice_sites{acceptor}{$splice1}++;}
			else {$genome_splice_sites{acceptor}{$splice1} = 1;}
		}

# 		print"$splice1\n";
	} elsif ($position eq "end"){
		# extract 2 bases before the start codon
		$splice2 = substr($sequence, $start - 3, 2); # wants to start 2 bases before the start codon
		if ($strand eq "+"){
			if (exists $genome_splice_sites{acceptor}{$splice2}) {$genome_splice_sites{acceptor}{$splice2}++;}
			else {$genome_splice_sites{acceptor}{$splice2} = 1;}
		}
		if ($strand eq "-"){
			$splice2 =~ tr/acgt/tgca/;
			$splice2 = reverse($splice2);
			if (exists $genome_splice_sites{donor}{$splice2}) {$genome_splice_sites{donor}{$splice2}++;}
			else {$genome_splice_sites{donor}{$splice2} = 1;}
		}

# 		print"$splice2\n";
	} else{
		# extract both 2 bases before start codon and 2 bases after stop codon
		$splice1 = substr($sequence, $stop, 2);
		$splice2 = substr($sequence, $start - 3, 2);
		if ($strand eq "+"){
			if (exists $genome_splice_sites{donor}{$splice1}) {$genome_splice_sites{donor}{$splice1}++;}
			else {$genome_splice_sites{donor}{$splice1} = 1;}
			if (exists $genome_splice_sites{acceptor}{$splice2}) {$genome_splice_sites{acceptor}{$splice2}++;}
			else {$genome_splice_sites{acceptor}{$splice2} = 1;}
		}
		if ($strand eq "-"){
			$splice1 =~ tr/acgt/tgca/;
			$splice1 = reverse($splice1);
			$splice2 =~ tr/acgt/tgca/;
			$splice2 = reverse($splice2);
			if (exists $genome_splice_sites{acceptor}{$splice1}) {$genome_splice_sites{acceptor}{$splice1}++;}
			else {$genome_splice_sites{acceptor}{$splice1} = 1;}
			if (exists $genome_splice_sites{donor}{$splice2}) {$genome_splice_sites{donor}{$splice2}++;}
			else {$genome_splice_sites{donor}{$splice2} = 1;}
		}
# 		print"$splice1 \t $splice2\n";
	}
}

__END__
# determines intron composition
sub intron_comp{
	
}