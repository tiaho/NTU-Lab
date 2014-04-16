#!/usr/bin/perl
#
# sanity_check.pl

# old script, refer to start_stop_codon_usage_cds_length.pl

# work on frameshifts
# how to determine codon bias with frameshifts?

use strict; use warnings;
use File::Slurp;

# hashes
my (%cds_seq, %start_codons, %stop_codons, %CDS_length, %codon_usage);

# runs the script for each of the 6 chromosomes of C elegans
# extract_data("xx01", "test.gff2");
extract_data("xx01", "chromI_coding_exon.gff2");
extract_data("xx02", "chromII_coding_exon.gff2");
extract_data("xx03", "chromIII_coding_exon.gff2");
extract_data("xx04", "chromIV_coding_exon.gff2");
extract_data("xx05", "chromV_coding_exon.gff2");
extract_data("xx06", "chromX_coding_exon.gff2");

# prints the start codons and counts
print "---Start codons---\n";
for my $key (keys %start_codons){
    print"$key, $start_codons{$key}\n";
}

# prints the stop codons and counts
# print "\n---Stop codons---\n";
# for my $key (keys %stop_codons){
#     print"$key, $stop_codons{$key}\n";
# }

# prints the codon length counts
# print "\n---Codon lengths---\n";
# for my $key (keys %CDS_length){
#     print"$key, $CDS_length{$key}\n";
# }

# # prints the codon usage counts
# print "\n---Codon usage---\n";
# for my $key (keys %codon_usage){
#     print"$key, $codon_usage{$key}\n";
# }

#################
###SUBROUTINES###
#################

# extracts data from a chromosome
sub extract_data{
    my ($fasta, $gff) = @_;
    
    # reads in the dna sequence
    my $sequence = read_file($fasta);
    $sequence =~ s/>[A-Z]+//; # removes the header specifying the chromosome
    $sequence =~ s/\n//g; # removes all new lines

    # extracts start/stop codon positions
    open(my $in, "<$gff") or die "Cannot open $gff\n";
    while (my $line = <$in>){
		chomp $line;
	
		# extracts the coding sequence
		extract_cds($sequence, $line);
    
	# checks start and stop coordinates
# 	check_start_stop_coord($start_pos, $stop_pos, $line);
    }
    close $in;

    for my $key (keys %cds_seq){
# 		print"$key\n";
# 		print"$key, $cds_seq{$key}{strand}\n";

		# extracts the start and stop codons
		extract_start_stop_codons($key);
	
		# gets a count of CDS lengths
# 		count_CDS_length($start_pos, $stop_pos);
    }

}

# extracts the coding sequence
sub extract_cds{
    my ($sequence, $line) = @_;
    my (undef, $source, $feature, $start_pos, $stop_pos, undef, $strand, undef, $attributes) = split("\t", $line);
    
    # only wants the coordinates from curated coding_exon
    if ($source eq "curated" and $feature eq "coding_exon"){
		my ($cdsID) = ($attributes =~ m/CDS \"([A-Za-z0-9_]+\.[0-9]+)/); # don't care about splice variants
		$cds_seq{$cdsID}{strand} = $strand;
		$cds_seq{$cdsID}{line} = $line;

		# forward strand
		if ($strand eq "+"){
			my $coding_exon = substr($sequence, $start_pos - 1, $stop_pos + 1 - $start_pos); # -1 b/c perl string starts at 0, + 1 to take account of last base
			$cds_seq{$cdsID}{cds} .= $coding_exon;
# 			print"$cdsID, $strand, $cds_seq{$cdsID}{cds}\n";
		# reverse strand
		}else{ # strand eq '-'
			my $coding_exon = substr($sequence, $start_pos - 1, $stop_pos + 1 - $start_pos); # -1 b/c perl string starts at 0, + 1 to take account of last base
			$coding_exon =~ tr/atcg/tagc/; # reverse it later, or else sequence will be out of order
			$cds_seq{$cdsID}{cds} .= $coding_exon;
# 			print"$cdsID, $strand, $cds_seq{$cdsID}{cds}\n";
		}
    }
}

# # checks start and stop coordinates
# sub check_start_stop_coord{
#     my ($start_pos, $stop_pos, $line) = @_;
#     if ($start_pos == $stop_pos){
# # 	print"Start coordinate is equal to the stop coordinate $start_pos = $stop_pos\n";
# 	print"$line\n";
# #     }elsif ($start_pos > $stop_pos){
# # 	print"Start coordinate is greater than the stop coordinate $start_pos > $stop_pos\n";
# # 	print"$line\n";
# #     }else{
# # 	print"$start_pos < $stop_pos\n";
#     }
# }

# extracts the start and stop codons
sub extract_start_stop_codons{
    my $key = $_[0];
#     print"$key\n";
    my ($start, $stop);
	if ($cds_seq{$key}{strand} eq "-") {$cds_seq{$key}{cds} = reverse($cds_seq{$key}{cds});}
	$start = substr($cds_seq{$key}{cds}, 0, 3);
	$stop = substr($cds_seq{$key}{cds}, -3);
	
# 	print"$start, $stop\n";
  
	
	# prints the lines that don't have expected start/stop codons
# 	if ($start ne "atg"){
# 		print"$start, $cds_seq{$key}{line}\n";
# # 		print"$length\n";
# 		print"$cds\n";
# 	}
# 	if ($stop ne "taa" and $stop ne "tag" and $stop ne "tga"){
# 		print"$stop, $cds_seq{$key}{line}\n";
# # 		print"$length\n";
# 		print"$cds\n";
# 	}
	
	# gets a count of each codon
	if (exists $start_codons{$start}){$start_codons{$start}++;}
	else{$start_codons{$start} = 1;}
	if (exists $stop_codons{$stop}){$stop_codons{$stop}++;}
	else{$stop_codons{$stop} = 1;}
    
	# prints out the line if start codon isn't ATG and if stop codon isn't TAA/TAG/TGA
#         if ($start ne "atg"){print"$cds_seq{$key}{line}\n";}
#         if ($stop ne "taa" and $stop ne "tag" and $stop ne "tga"){print"$cds_seq{$key}{line}\n";}

# 	# checks for frameshifts
# 	check_frameshift($key);

	# checks for codon bias
# 	check_codon_bias($start_pos, $stop_pos, $sequence, $key);
}

__END__

# # checks for frameshifts
# sub check_frameshift{
#     my ($key) = $_;
#     my $start_pos = $cds_seq{$key}{start};
#     my $stop_pos = $cds_seq{$key}{stop};
#     my $modulo = ($stop_pos + 1 - $start_pos) % 3;
#     if ($modulo != 0){print "$modulo\n";}
#     else{print"okay\n";}
#     
# }

# # checks for codon bias
# sub check_codon_bias{
#     my ($start_pos, $stop_pos, $sequence, $key) = @_;
#     my $seq_length = $stop_pos + 1 - $start_pos;
#     my $CDS = substr($sequence, $start_pos - 1, $seq_length);
#     if ($cds_seq{$key}{strand} eq "-"){
# 	$CDS =~ tr/atgc/tacg/;
# 	$CDS = reverse($CDS);
#     }
#     for (my $i = 1; $i <= $seq_length; $i += 3){
# 	my $codon = substr($CDS, $i - 1, 3);
# 	if (exists $codon_usage{$codon}){$codon_usage{$codon}++;}
# 	else{$codon_usage{$codon} = 1;}
#     }
# }

# gets a count of CDS lengths
sub count_CDS_length{
    my ($start_pos, $stop_pos) = @_;
    my $length = abs($stop_pos + 1 - $start_pos);
#     if (exists $CDS_length{$length}) {$CDS_length{$length}++;}
#     else{$CDS_length{$length} = 1;}
    print"$length\n";
}