#!/usr/bin/perl
#
# split_gff_file.pl
# splits the gff file by chromosome
# REMINDER change outfile names each time
use strict; use warnings;


# file to be split
my $file = $ARGV[0];
die "Please enter a gff file to be split\n" unless (@ARGV == 1);
open (my $in, "<$file") or die "Cannot open $file\n";

# gff file split by chromosome
open (my $out1, ">chromI_master.gff2");
open (my $out2, ">chromII_master.gff2");
open (my $out3, ">chromIII_master.gff2");
open (my $out4, ">chromIV_master.gff2");
open (my $out5, ">chromV_master.gff2");
open (my $out6, ">chromX_master.gff2");

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
