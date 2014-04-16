#!/usr/bin/perl -w
#
# coordinate_extracter.pl
# <Takes a wormbase ID, finds the gene gene entry for that ID then prints out all the lines of the file that are in that genes range of coordinates>
# Input: gff file and a wormbase ID and the number of bases that are allowed to be off but still be in the gff
# Output: all of the lines of the gff file that correlate to a specific gene.

use strict;
use warnings FATAL => 'all';
use feature ':5.10';

my($file, $id);

my $number_off = 0;

#allows the use of either 2 or three arguements. User can enter in #number_off or it is otherwise 0.
if(@ARGV == 2){
	($file, $id) = @ARGV;
}
elsif(@ARGV == 3){
	($file, $id, $number_off) = @ARGV;
}
else{
	die "specify two or three arguements: a file name, a wormbase id and optionally a number of bases that the gff file can contain that are off the gene.";
}


# need to track chromosome, start and stop coordinate for each gene
my $gene_chr;
my $gene_start;
my $gene_stop;


my $output_file = "${id}_${number_off}.gff";
if(-e $output_file and -s $output_file){
	warn "$output_file already exists, skipping\n";
	exit;
}



#go through the gff file the first time to find the coordinates of the gene gene entry
open(my $in, "< $file ") or die "Can't open the file";
while(my $line = <$in>){
	# skip comments
	next if ($line =~ m/^#/);
	
	my ($location, $source, $feature, $start, $stop, undef, undef, undef, $extra)  = split("\t", $line);
	if($source eq 'gene' and $feature eq 'gene' and $extra =~ m/$id/){
		$gene_start = $start - $number_off;
		$gene_stop  = $stop + $number_off;	
		$gene_chr = $location;
		last;
	}	
}
close($in);


#print all lines within the gffcoordinates

open(my $in2, "<", "$file") or die "Can't open the $file";
open(my $out, ">", $output_file) or die "Error: can't create $output_file";
while(my $line = <$in2>){

	# skip comments
	next if ($line =~ m/^#/);

	my ($location, undef, undef, $start, $stop, undef, undef, undef, undef)  = split("\t", $line);

	# not worth comparing coordiantes if the current feature is not from the same chromosome
	# as target gene
	next if ($location ne $gene_chr);

	if($start >= $gene_start and $stop <= $gene_stop or check_overlap($start,$stop,$gene_start, $gene_stop)){
		print $out "$line";
	}

	# take advantage of the fact that we have sorted GFF file by start coordiantes
	last if (($start > $gene_stop) and ($location eq $gene_chr));
}
close($in2);
close($out);


exit;


sub check_overlap{
#returns (0) if no overlap exists between the two objects
#returns (1) if there is any kind of overlap at even one point.

	my ($start1 , $stop1, $start2, $stop2) = @_;
#makes the first and second element of results a 0 if there is no overlap;
	if($start1 > $stop2 or $start2 > $stop1){
		return 0;
	}
	else{	
		return 1;
	}
}


