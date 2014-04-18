#!/usr/bin/perl
#
# make_source_feature_file.pl
# makes a file containing only the multiple desired source/feature combinations

use strict; use warnings;


my ($gff) = @ARGV;
die "<gff2>\n" unless @ARGV == 1;

# reads in the gff2 file
open(my $in, "<$gff") or die "Cannot open $gff\n";

my $file = "gene_Coding_transcript_exon.gff2";
# writes the matching source/feature lines out
open(my $out, ">$file") or die "Cannot open $file\n";
while (my $line = <$in>){
    chomp $line;
    if ($line =~ m/^##/) {next;}
    my (undef, $source, $feature, undef, undef, undef, undef, undef, undef) = split("\t", $line);
    if ((($source eq "gene") and ($feature eq "gene")) or (($source eq "Coding_transcript") and ($feature eq "exon"))){
		print $out "$line\n";
    }
}

close $in;
close $out;