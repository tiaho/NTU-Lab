#!/usr/bin/perl
#
# split_gff_file.pl
# splits a gff file by chromosome (currently hard coded for *C. elegans*). need to manually change file name each time!

use strict; use warnings;


# file to be split
my $file = $ARGV[0];
die "Please enter a gff file to be split\n" unless (@ARGV == 1);
open (my $in, "<$file") or die "Cannot open $file\n";

# gff file split by chromosome
open (my $out1, ">chromI_gene_gene.gff2");
open (my $out2, ">chromII_gene_gene.gff2");
open (my $out3, ">chromIII_gene_gene.gff2");
open (my $out4, ">chromIV_gene_gene.gff2");
open (my $out5, ">chromV_gene_gene.gff2");
open (my $out6, ">chromX_gene_gene.gff2");

while (my $line = <$in>){
    chomp $line;
    my ($chrom, undef) = split("\t", $line);
    if ($chrom eq "I") {print $out1 "$line\n";}
    elsif ($chrom eq "II") {print $out2 "$line\n";}
    elsif ($chrom eq "III") {print $out3 "$line\n";}
    elsif ($chrom eq "IV") {print $out4 "$line\n";}
    elsif ($chrom eq "V") {print $out5 "$line\n"}
    elsif ($chrom eq "X") {print $out6 "$line\n"}
}

close $in;
close $out1;
close $out2;
close $out3;
close $out4;
close $out5;
close $out6;
