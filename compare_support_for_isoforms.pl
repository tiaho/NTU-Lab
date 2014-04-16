#!/usr/bin/perl
#
# compare_support_for_isoforms.pl
#
# <runs the list of wormbase gene IDs, extracts lines in the file that are in the coordinate range of the gene,
#    then finds the coding exons of the gene and compares them for each isoform to see how they overlap with the other source/feature combos>
# Input: file with gene IDs, gff file, parameters file
# Output: 

use strict; use warnings FATAL => 'all';
use File::Copy;

die "Usage: <gene ID file> <GFF file> <parameters file> <number of genes to be run>" unless @ARGV == 4;

my ($gene_ID_file, $gff_file, $parameter_file, $n) = @ARGV;
my %gene_scores; # key: gene ID. value: score
my $count = 0;
my $abs_score_threshold = 50;

open(my $in, "<$gene_ID_file"); #opens the gene ID file
while(my $gene_ID = <$in>){
		chomp($gene_ID);

		my ($abs_score, $rel_score);
		my $gene_gff = $gene_ID."_0.gff";  #GFF file for the particular gene
		my $gene_gff_location = "GFF_files/".$gene_gff;
		
		# make GFF file if it doesn't exist or if it zero bytes
		if(not -e $gene_gff_location or -z $gene_gff_location){
				
		    # gets data for the gene that is within the coordinates of the gene
			my $command = "coordinate_extracter.pl $gff_file $gene_ID";
			system($command) == 0 or die "Cannot run $command\n";
			move($gene_gff, $gene_gff_location) or die "Failed to move file"; # moves output file into subdirectory
		}
		
	
#		print"$gene_ID\n";
		# gets the score for the gene
		my $command = "find_overlapping_features2.pl $gene_gff_location $gene_ID $parameter_file > output.tmp";
		system($command) == 0 or die "Cannot run $command"; # ranks the isoforms of that gene
		
		# captures the score and stores it in hash
		# consider using grep?
#		my $rel = `grep \"^Absolute\" $file`;
		open(my $in2, "<output.tmp"); # opens the temp file
		while(my $line = <$in2>){
			chomp($line);
			if($line =~ m/Absolute difference between scores of best\/worst isoform\: ([\d+\.]+)/){
				$abs_score = $1;
				#print"$abs_score\n";
			} elsif($line =~ m/Relative difference between scores of best\/worst isoform\: ([\d\.]+)/){
				$rel_score = $1;
			} else{
				# do nothing	
			}
		}
		close($in2);

# 		my $score = "$rel_score\t$abs_score";

		print "$gene_ID $rel_score $abs_score\n";

# 		if($rel_score >= 0.70 and $abs_score >= $abs_score_threshold){
# 		$gene_scores{$gene_ID} = $score;
		my $new_output = $gene_ID.".txt";
		my $new_output_location = "Output_files/".$new_output;
		rename("output.tmp", $new_output) or die "Cannot rename file";
		move($new_output, $new_output_location) or die "Failed to move file";
# 		} 
		
		if($count == int(0.10*$n)){
			warn("10% done\n"); 
		}elsif($count == int(0.20*$n)){
			warn("20% done\n"); 
		}elsif($count == int(0.30*$n)){
			warn("30% done\n"); 
		}elsif($count == int(0.40*$n)){
			warn("40% done\n"); 
		}elsif($count == int(0.50*$n)){
			warn("50% done\n"); 	
		}elsif($count == int(0.60*$n)){
			warn("60% done\n");
		}elsif($count == int(0.70*$n)){
			warn("70% done\n"); 
		}elsif($count == int(0.80*$n)){
			warn("80% done\n"); 	
		}elsif($count == int(0.90*$n)){
			warn("90% done\n"); 
		}
		
		$count++;
		last if ($count == $n);
}
close($in);

print "\n";

# prints the gene ID and score
# foreach my $key (sort{$gene_scores{$b} cmp $gene_scores{$a}} keys %gene_scores) {
# 	print "$key $gene_scores{$key}\n";
# }
