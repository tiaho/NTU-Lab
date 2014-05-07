#!/usr/bin/perl
#
# splice_site_composition.pl
#
# calculates the composition of the splice sites
# uses source/feature: curated exon

use strict; use warnings;
use File::Slurp;
use Data::Dumper;

# calculate the frequency (percentage) of each letter at a given spot (from genomic data)
# calculate the frequency we'd expect to see each of the actual given combinations (for each donor/acceptor site)
# average that frequency for each CDS
# z-scores of how different each of them are from the mean frequency

# user either enters "donor" or "acceptor" to indicate which info they want
die "Usage: splice_site_composition.pl <donor OR acceptor>" unless (@ARGV == 1) and ( ($ARGV[0] eq "donor") or ($ARGV[0] eq "acceptor") );
my $desired_type = $ARGV[0];

# hashes
my (%exon_coords, %splice_sites, %splice_site_probability, %genome_splice_site_count, %genome_splice_site_freq);

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

# extracts the splice sites
for my $cdsID (keys %exon_coords){
	my $chrom_seq;

	if ($exon_coords{$cdsID}{count} == 1) {next;}
	
	# assigns the correct genomic sequence for each chromosome
	if ($exon_coords{$cdsID}{chromosome} eq "I") {$chrom_seq = $sequences[0];}
	elsif ($exon_coords{$cdsID}{chromosome} eq "II") {$chrom_seq = $sequences[1];}
	elsif ($exon_coords{$cdsID}{chromosome} eq "III") {$chrom_seq = $sequences[2];}
	elsif ($exon_coords{$cdsID}{chromosome} eq "IV") {$chrom_seq = $sequences[3];}
	elsif ($exon_coords{$cdsID}{chromosome} eq "V") {$chrom_seq = $sequences[4];}
	elsif ($exon_coords{$cdsID}{chromosome} eq "X") {$chrom_seq = $sequences[5];}

	for (my $i = 1; $i <= $exon_coords{$cdsID}{count}; $i++){
		# extracts splice sites
		my $position;
		if ($i == 1) {$position = "beg";}
		elsif ($i == $exon_coords{$cdsID}{count}) {$position = "end";}
		else {$position = "middle";}
		extract_splice_sites($position, $chrom_seq, $cdsID, $exon_coords{$cdsID}{strand}, $exon_coords{$cdsID}{$i}{start}, $exon_coords{$cdsID}{$i}{stop});
	}
	
}
# print Dumper %splice_sites;

# calculates the counts for each nucleotide for the whole genome
for my $cdsID (keys %splice_sites){
	for my $type (keys %{ $splice_sites{$cdsID} }){
		for my $site (keys %{ $splice_sites{$cdsID}{$type} }){
			for (my $i = 0; $i <= 4; $i++){
				my $nuc = substr($site, $i, 1);
				if (exists $genome_splice_site_count{$type}{$i + 1}{$nuc}) {$genome_splice_site_count{$type}{$i + 1}{$nuc} += $splice_sites{$cdsID}{$type}{$site};}
				else {$genome_splice_site_count{$type}{$i + 1}{$nuc}++;}
			}
		}
	}
}
# print Dumper %genome_splice_site_count;

# calculates the frequency of each nucleotide at each position for splice donors/acceptors
for my $type (keys %genome_splice_site_count){
	for (my $i = 1; $i <= 5; $i++){
		my $a = $genome_splice_site_count{$type}{$i}{a};
		my $c = $genome_splice_site_count{$type}{$i}{c};
		my $g = $genome_splice_site_count{$type}{$i}{g};
		my $t = $genome_splice_site_count{$type}{$i}{t};
		my $total = $a + $c + $g + $t;

		$genome_splice_site_freq{$type}{$i}{a} = $a/$total;
		$genome_splice_site_freq{$type}{$i}{c} = $c/$total;
		$genome_splice_site_freq{$type}{$i}{g} = $g/$total;
		$genome_splice_site_freq{$type}{$i}{t} = $t/$total;
	}
}
# print Dumper %genome_splice_site_freq;

# calculates the probability of finding a splice site
for my $cdsID (keys %splice_sites){
	for my $type (keys %{$splice_sites{$cdsID} }){
		for my $site (keys %{ $splice_sites{$cdsID}{$type} }){
			my (@nuc, @freq);
			for (my $i = 0; $i <= 4; $i++){
				$nuc[$i] = substr($site, $i, 1);
				$freq[$i] = $genome_splice_site_freq{$type}{$i + 1}{$nuc[$i]};
			}
			my $num_of_occurances = $splice_sites{$cdsID}{$type}{$site};
			my $site_freq = $freq[0] * $freq[1] * $freq[2] * $freq[3] * $freq[4] * $num_of_occurances; # all the individual frequencies for each site and the number of times that site occurs
			
			if (exists $splice_site_probability{$cdsID}{$type}{freq}) {$splice_site_probability{$cdsID}{$type}{freq} += $site_freq;}
			else {$splice_site_probability{$cdsID}{$type}{freq} = $site_freq;}
			if (exists $splice_site_probability{$cdsID}{$type}{count}) {$splice_site_probability{$cdsID}{$type}{count} += $num_of_occurances;}
			else {$splice_site_probability{$cdsID}{$type}{count} = $num_of_occurances;}
		}
	}
}
# print Dumper %splice_site_probability;

# calculates the average of the probabilities for each CDS
for my $cdsID (keys %splice_site_probability){
	for my $type (keys % {$splice_site_probability{$cdsID} }){
		my $freq = $splice_site_probability{$cdsID}{$type}{freq};
		my $count = $splice_site_probability{$cdsID}{$type}{count};
		$splice_site_probability{$cdsID}{$type}{avg} = $freq/$count;
	}
}
# print Dumper %splice_site_probability;

# prints the cdsID, chromosome, and avg frequency of finding the splice site for a certain user specified type (donor or acceptor)
for my $cdsID (keys %splice_site_probability) {
	my $chromosome = $exon_coords{$cdsID}{chromosome};
	my $avg = $splice_site_probability{$cdsID}{$desired_type}{avg};
	print"$avg, $chromosome, $cdsID\n";
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
sub extract_splice_sites{
	my ($position, $sequence, $cdsID, $strand, $start, $stop) = @_;
	
	my ($splice1, $splice2);
	
	if ($position eq "beg"){
		# extract 5 bases after the stop codon
		$splice1 = substr($sequence, $stop, 5); # wants to start at the base after the stop codon
		if ($strand eq "+"){
			# stores in a hash to keep count for each CDS
			if (exists $splice_sites{$cdsID}{donor}{$splice1}) {$splice_sites{$cdsID}{donor}{$splice1}++;}
			else {$splice_sites{$cdsID}{donor}{$splice1} = 1;}
		}
		if ($strand eq "-"){
			$splice1 =~ tr/acgt/tgca/;
			$splice1 = reverse($splice1);
			# stores in a hash to keep count for each CDS
			if (exists $splice_sites{$cdsID}{acceptor}{$splice1}) {$splice_sites{$cdsID}{acceptor}{$splice1}++;}
			else {$splice_sites{$cdsID}{acceptor}{$splice1} = 1;}
		}

# 		print"$splice1\n";
	} elsif ($position eq "end"){
		# extract 2 bases before the start codon
		$splice2 = substr($sequence, $start - 6, 5); # wants to start 2 bases before the start codon
		if ($strand eq "+"){
			# stores in a hash to keep count for each CDS
			if (exists $splice_sites{$cdsID}{acceptor}{$splice2}) {$splice_sites{$cdsID}{acceptor}{$splice2}++;}
			else {$splice_sites{$cdsID}{acceptor}{$splice2} = 1;}
		}
		if ($strand eq "-"){
			$splice2 =~ tr/acgt/tgca/;
			$splice2 = reverse($splice2);
			# stores in a hash to keep count for each CDS
			if (exists $splice_sites{$cdsID}{donor}{$splice2}) {$splice_sites{$cdsID}{donor}{$splice2}++;}
			else {$splice_sites{$cdsID}{donor}{$splice2} = 1;}
		}

# 		print"$splice2\n";
	} else{
		# extract both 2 bases before start codon and 2 bases after stop codon
		$splice1 = substr($sequence, $stop, 5);
		$splice2 = substr($sequence, $start - 6, 5);
		if ($strand eq "+"){
			# stores in a hash to keep count for each CDS
			if (exists $splice_sites{$cdsID}{donor}{$splice1}) {$splice_sites{$cdsID}{donor}{$splice1}++;}
			else {$splice_sites{$cdsID}{donor}{$splice1} = 1;}
			if (exists $splice_sites{$cdsID}{acceptor}{$splice2}) {$splice_sites{$cdsID}{acceptor}{$splice2}++;}
			else {$splice_sites{$cdsID}{acceptor}{$splice2} = 1;}
		}
		if ($strand eq "-"){
			$splice1 =~ tr/acgt/tgca/;
			$splice1 = reverse($splice1);
			$splice2 =~ tr/acgt/tgca/;
			$splice2 = reverse($splice2);
			# stores in a hash to keep count for each CDS
			if (exists $splice_sites{$cdsID}{acceptor}{$splice1}) {$splice_sites{$cdsID}{acceptor}{$splice1}++;}
			else {$splice_sites{$cdsID}{acceptor}{$splice1} = 1;}
			if (exists $splice_sites{$cdsID}{donor}{$splice2}) {$splice_sites{$cdsID}{donor}{$splice2}++;}
			else {$splice_sites{$cdsID}{donor}{$splice2} = 1;}
		}
# 		print"$splice1 \t $splice2\n";
	}
}