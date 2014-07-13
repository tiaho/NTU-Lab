#!/usr/bin/perl
#
# make_source_feature_file.pl
# makes a file that contains only the single desired source/feature combination

use strict; use warnings;


my ($gff, $desired_source, $desired_feature) = @ARGV;
die "<gff2> <source> <feature>\n" unless @ARGV == 3;

# reads in the gff2 file
open(my $in, "<$gff") or die "Cannot open $gff\n";

# resulting file name is source_feature.gff2
my $file = $desired_source . '_' . $desired_feature . '.gff2';

# writes the matching source/feature lines out
open(my $out, ">$file") or die "Cannot open $file\n";
while (my $line = <$in>){
    chomp $line;
    if ($line =~ m/^##/) {next;}
    my (undef, $source, $feature, undef, undef, undef, undef, undef, undef) = split("\t", $line);
    if (($source eq $desired_source) and ($feature eq $desired_feature)){
		print $out "$line\n";
    }
}

close $in;
close $out;