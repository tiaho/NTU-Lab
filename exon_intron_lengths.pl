#!/usr/bin/perl
#
# exon_intron_lengths.pl
#
# calculates the lengths of exons and introns - current working version 
#
# uses these source/feature types: curated exon
#
# prints the exon lengths
# prints the intron lengths


use strict; use warnings;
use File::Slurp;

# the input files
my @gff_files = ("chromI_curated_exon.gff2", "chromII_curated_exon.gff2", "chromIII_curated_exon.gff2", "chromIV_curated_exon.gff2", "chromV_curated_exon.gff2", "chromX_curated_exon.gff2");

# runs the script for each of the 6 chromosomes of C elegans
for (my $i = 0; $i <= 5; $i++){
	makes_file($gff_files[$i]);
	my $extracted_file_name = "extracted_" . $gff_files[$i];
	my $sorted_file_name = "sorted_" . $extracted_file_name;
	sorts_file($extracted_file_name);
	get_intron_exon_lengths($sorted_file_name);
	
# 	for my $line (@sorted_file){
# 		print"$line";
# 	}
}

#################
###SUBROUTINES###
#################

# makes files with only 6 columns (chromosome, source, feature, start, stop, id), with id being the extracted cds id
sub makes_file{
    my $gff = $_[0];
	my $outfile = "extracted_" . $gff;

    # extracts start/stop codon positions
    open(my $in, "<$gff") or die "Cannot open $gff\n";
    open(my $out, ">$outfile") or die "Cannot open $outfile\n";
    
    while (my $line = <$in>){
		chomp $line;
		my ($chromosome, $source, $feature, $start, $stop, undef, undef, undef, $attributes) = split("\t", $line);
		
		# extracts the ID
		my ($id) = ($attributes =~ m/CDS \"([A-Za-z0-9_]+\.[a-z0-9]+\.?[a-z0-9]*)/);
		
		# prints to outfile
		print $out "$chromosome $source $feature $start $stop $id\n";
    }
    close $in;
    close $out;
}

# sorts the files by 1st cds id then 2nd by start coordinate
sub sorts_file{
	my $file = $_[0];
	my $sorted_file_name = "sorted_" . $file;
	my $command = "sort -k6,6 -k4,4n $file > $sorted_file_name"; # sorts the start coords and then stop by reverse
	system($command) == 0 or die "Cannot run $command";
}

# gets the exon and intron lengths
sub get_intron_exon_lengths{
	my $in_file = $_[0];
	open(my $in, "<$in_file") or die "Cannot open $in_file\n";
	
	my @file = <$in>;
	my $file_length = scalar(@file);
	
	for(my $i = 0; $i <= $file_length - 1; $i++){ # -1 because array indices start at 0
		die "$i\n" unless chomp $file[$i];
		my ($chromosome, $source, $feature, $start, $stop, $id) = split(" ", $file[$i]);
		
		# calculates exon length
		my $exon_length = $stop + 1 - $start; # +1 because want to include both ends of the exon
# 		print"$exon_length, $chromosome, $id\n";
		
		# calculates intron length
		if ($i == 0){
			# do nothing
		}else{
			my (undef, undef, $prev_feature, $prev_start, $prev_stop, $prev_id) = split(" ", $file[$i - 1]);
			if ($id eq $prev_id){
				my $intron_length = $start - $prev_stop - 1; # -1 b/c doesn't include the exon start/end
				print"$intron_length, $chromosome, $id\n";
			}else{ # id and prev id aren't the same
				# do nothing
			}
		}
	}
}