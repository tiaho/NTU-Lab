#!/usr/bin/perl
#
# exon_intron_lengths.pl
#
#  uses these source/feature types: gene gene, curated exon
#
# prints the exon lengths
# prints the intron lengths


use strict; use warnings;
use File::Slurp;

# the input files
my @gff_files = ("chromI_gene_exon.gff2", "chromII_gene_exon.gff2", "chromIII_gene_exon.gff2", "chromIV_gene_exon.gff2", "chromV_gene_exon.gff2", "chromX_gene_exon.gff2");

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

# makes files with only 6 columns (chromosome, source, feature, start, stop, id), with id being the extracted cds/sequence id
sub makes_file{
    my $gff = $_[0];
	my $outfile = "extracted_" . $gff;

    # extracts start/stop codon positions
    open(my $in, "<$gff") or die "Cannot open $gff\n";
    open(my $out, ">$outfile") or die "Cannot open $outfile\n";
    
    while (my $line = <$in>){
		chomp $line;
		my ($chromosome, $source, $feature, $start, $stop, undef, undef, undef, $attributes) = split("\t", $line);
		
		my ($id, $biotype);
		# extracts the ID
		if ($source eq "gene" and $feature eq "gene"){
			($id) = ($attributes =~ m/Sequence_name \"([A-Za-z0-9_]+\.[a-z0-9]+)\"/);
			($biotype) = ($attributes =~ m/Biotype\"(\w+,*\w+)\"/);
		}else{
			($id) = ($attributes =~ m/CDS \"([A-Za-z0-9_]+\.[a-z0-9]+\.?[a-z0-9]*)/);
			$biotype = 0;
		} 
		
		# prints to outfile
		print $out "$chromosome $source $feature $start $stop $id $biotype\n";
    }
    close $in;
    close $out;
}

# sorts the files by 1st cds/sequence id then 2nd by start coord
sub sorts_file{
	my $file = $_[0];
	my $sorted_file_name = "sorted_" . $file;
	my $command = "sort -k6,6d -k4,4n -k5,5nr $file > $sorted_file_name"; # sorts the start coords and then stop by reverse
	system($command) == 0 or die "Cannot run $command";
}

# puts all gene gene start/stop coodinates and id into hash
sub get_intron_exon_lengths{
	my $in_file = $_[0];
	open(my $in, "<$in_file") or die "Cannot open $in_file\n";
	
	my @file = <$in>;
	my $file_length = scalar(@file);
	
	my ($holder_geneID, $holder_genestart, $holder_genestop);
	for(my $i = 0; $i <= $file_length - 1; $i++){ # -1 because array indices start at 0
		die "$i\n" unless chomp $file[$i];
		my ($chromosome, $source, $feature, $start, $stop, $id, $biotype) = split(" ", $file[$i]);
		
		if ($i == 0){
			if (($feature eq "gene") and ($biotype !~ m/RNA/)){
				$holder_geneID = $id;
				$holder_genestart = $start;
				$holder_genestop = $stop;
			}
		}else{
			my (undef, undef, $prev_feature, $prev_start, $prev_stop, $prev_id, $prev_biotype) = split(" ", $file[$i - 1]);
			if ($feature eq "gene"){
				if ($prev_feature eq "exon"){
					my $intron_length = $holder_genestop - $prev_stop - 1; # -1 b/c doesn't include the exon start/end
					print"$intron_length, $chromosome, $prev_id\n";
				}
				if ($biotype !~ m/RNA/){
					$holder_geneID = $id;
					$holder_genestart = $start;
					$holder_genestop = $stop;
					if ($prev_feature eq "gene"){
						my $intron_length = $prev_stop - $prev_start - 1; # -1 b/c doesn't include the exon start/end
						print"$intron_length, $chromosome, $prev_id\n";
					}
				}elsif ($biotype =~ m/RNA/){
					if ($prev_biotype =~m/RNA/){
						my $intron_length = $prev_stop - $prev_start - 1; # -1 b/c doesn't include the exon start/end
						print"$intron_length, $chromosome, $prev_id\n";
					}
				}
			}elsif ($feature eq "exon"){
				# get exon lengths
				my $exon_length = $stop + 1 - $start; # +1 because want to include both ends of the exon
# 				print"$exon_length, $chromosome, $id\n";
				if ($prev_feature eq "gene"){
					my $intron_length = $start - $holder_genestart - 1; # -1 b/c doesn't include the exon start/end
					print"$intron_length, $chromosome, $id\n";
				}elsif(($prev_feature eq "exon") and ($prev_id eq $id)){ # ($prev_id eq $id)
					my $intron_length = $start - $prev_stop - 1; # -1 b/c doesn't include the exon start/end
					print"$intron_length, $chromosome, $id\n";
				}elsif(($prev_feature eq "exon") and ($prev_id ne $id)){
					# intron between current exon start and gene start
					my $intron_length = $start - $holder_genestart - 1; # -1 b/c doesn't include the exon start/end
					print"$intron_length, $chromosome, $id\n";
					# intron between previous exon stop and gene stop
					$intron_length = $holder_genestop - $prev_stop - 1; # -1 b/c doesn't include the exon start/end
					print"$intron_length, $chromosome, $prev_id\n";
				}
			}
		}
	}
}

# gets the intron and exon lengths
# sub get_intron_exon_lengths{
# 	my $in_file = $_[0];
# 	open(my $in, "<$in_file") or die "Cannot open $in_file\n";
# 	
# 	my @file = <$in>;
# 	my $file_length = scalar(@file);
# 	
# 	my ($holder_geneID, $holder_genestart, $holder_genestop);
# 	for(my $i = 0; $i <= $file_length - 1; $i++){ # -1 because array indices start at 0
# 		die "$i\n" unless chomp $file[$i];
# 		my ($chromosome, $source, $feature, $start, $stop, $id, $biotype) = split(" ", $file[$i]);
# 		
# 		if ($i == 0){
# 			if (($feature eq "gene") and ($biotype !~ m/RNA/)){
# 				$holder_geneID = $id;
# 				$holder_genestart = $start;
# 				$holder_genestop = $stop;
# 			}
# 		}else{
# 			my (undef, undef, $prev_feature, $prev_start, $prev_stop, $prev_id, $prev_biotype) = split(" ", $file[$i - 1]);
# 			if ($feature eq "gene"){
# 				if ($prev_feature eq "exon"){
# 					my $intron_length = $holder_genestop - $prev_stop - 1; # -1 b/c doesn't include the exon start/end
# 					print"$intron_length, $chromosome, $prev_id\n";
# 				}
# 				if ($biotype !~ m/RNA/){
# 					$holder_geneID = $id;
# 					$holder_genestart = $start;
# 					$holder_genestop = $stop;
# 					if ($prev_feature eq "gene"){
# 						my $intron_length = $prev_stop - $prev_start - 1; # -1 b/c doesn't include the exon start/end
# 						print"$intron_length, $chromosome, $prev_id\n";
# 					}
# 				}elsif ($biotype =~ m/RNA/){
# 					if ($prev_biotype =~m/RNA/){
# 						my $intron_length = $prev_stop - $prev_start - 1; # -1 b/c doesn't include the exon start/end
# 						print"$intron_length, $chromosome, $prev_id\n";
# 					}
# 				}
# 			}elsif ($feature eq "exon"){
# 				# get exon lengths
# 				my $exon_length = $stop + 1 - $start; # +1 because want to include both ends of the exon
# # 				print"$exon_length, $chromosome, $id\n";
# 				if ($prev_feature eq "gene"){
# 					my $intron_length = $start - $holder_genestart - 1; # -1 b/c doesn't include the exon start/end
# 					print"$intron_length, $chromosome, $id\n";
# 				}elsif(($prev_feature eq "exon") and ($prev_id eq $id)){ # ($prev_id eq $id)
# 					my $intron_length = $start - $prev_stop - 1; # -1 b/c doesn't include the exon start/end
# 					print"$intron_length, $chromosome, $id\n";
# 				}elsif(($prev_feature eq "exon") and ($prev_id ne $id)){
# 					# intron between current exon start and gene start
# 					my $intron_length = $start - $holder_genestart - 1; # -1 b/c doesn't include the exon start/end
# 					print"$intron_length, $chromosome, $id\n";
# 					# intron between previous exon stop and gene stop
# 					$intron_length = $holder_genestop - $prev_stop - 1; # -1 b/c doesn't include the exon start/end
# 					print"$intron_length, $chromosome, $prev_id\n";
# 				}
# 			}
# 		}
# 	}
# }