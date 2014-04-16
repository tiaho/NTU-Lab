#!/usr/bin/perl
#
# find_overlapping_features.pl
#
# A script to take a GFF file for a single gene, and rank all isoforms for that 
# gene based on a user-specified list of desired GFF source/feature combinations
# 
#
# Author: Keith Bradnam, Genome Center, UC Davis
# based on code written by Ben Edwards
#
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.
#
# Version Number : 0.2

use strict;
use warnings FATAL => 'all';
use feature ':5.10';
use Data::Dumper;

die "Usage: <GFF file> <WormBase Gene ID> <parameter file>" unless (@ARGV == 3);
my ($gff_file, $gene_id, $param_file) = @ARGV;

# master multi-dimensional hash to hold (almost) all info from the gene's point of view
my %genes;

#        C27A2.2a  coding_exon            0
# $genes{$isoform}{$gff_feature}{exon_data}[$count]{start}
#                                                  {stop}
#                               {total_length}
#                               {results}{$target}{unique_hits}
#                                                 {hits_per_Kbp}
#                               {final_score}



# and another big multi-dimensional hash to store coordinates of all features of interest in a GFF file
my %features;
#          BLAT_EST_BEST  yk1051a06.3             0
# $features{$type}        {$feature_ID}[$count]{start}
# $features{$type}        {$feature_ID}[$count]{stop}


# need a hash to keep track of desired GFF source/feature combos to see if they match coding_exons or UTRs
# this gets populated from a file
my %desired_GFF_features;
load_desired_GFF_details();


# now need to get a list of all CDSs in the GFF file to determine the CDS identifiers 
# (will be used to match lines later on). Will cheat a bit bit and use command-line grep
my %CDSs;

foreach my $cds (`grep $gene_id $gff_file | grep curated | grep CDS`){
	chomp($cds);
	my ($location, $source, $feature, $start, $stop, undef, undef, undef, $comment)  = split("\t", $cds);
	my ($cds_id) = $comment =~ m/(\w+\.\d+[a-z]\.?\d?)/;
	$CDSs{$cds_id} = 1;
}




####################################################################################################
#
# 1st loop
#
# Parse GFF file to extract details of: coding exons, UTR regions, and CDS identifiers
#
####################################################################################################

open(my $gff, "< $gff_file ") or die "Can't open $gff_file\n";
while(my $line = <$gff>){
	chomp($line);

	my ($location, $source, $feature, $start, $stop, undef, undef, undef, $comment)  = split("\t", $line);

	# only want to look at Coding transcript lines which are coding_exons or UTRs
	next unless ($source eq 'Coding_transcript');
	next unless ($feature eq 'five_prime_UTR' or $feature eq 'three_prime_UTR' or $feature eq 'coding_exon');

	# extract Transcript ID
	my ($transcript_id) = $comment =~ m/Transcript \"(\w+\.\d*[a-z]?\.?\d*)\"/;

	# need to truncate $transcript_id when looking for a match to deal with
	# non-coding variants
	my $short_id = $transcript_id;
	if ($transcript_id =~ m/\d$/){
		$short_id =~ s/\.\d+$//; 
	}
	
	# does this transcript ID match the CDS IDs we extracted earlier?
	next unless (exists $CDSs{$short_id});

	# need to keep track of count of how many exon/UTR features we have seen for this isoform
	# if this is the first, set $count to 0;
	my $count = 0;
	if (defined ${genes{$transcript_id}{$feature}{exon_data}}){
		$count = scalar @{$genes{$transcript_id}{$feature}{exon_data}};
	} 
		
	# now add data to main hash
	$genes{$transcript_id}{$feature}{exon_data}[$count]{start}  = $start;
	$genes{$transcript_id}{$feature}{exon_data}[$count]{stop}   = $stop;

	# increment sum length of all coding exons or UTRs
	$genes{$transcript_id}{$feature}{total_length} += ($stop - $start + 1)

}
close($gff);




####################################################################################################
#
# 2nd loop
#
# Open gff file again. This time to capture details of all features of interest
#
####################################################################################################

open($gff, "< $gff_file ") or die "Can't open $gff_file\n";

while(my $line = <$gff>){
	chomp($line);
	my ($location, $source, $feature, $s_start, $s_stop, undef, undef, undef, $comment)  = split("\t", $line);

	# only want to consider GFF source/features that are in our list of desired features
	next unless (exists $desired_GFF_features{"$source\t$feature"});
	next unless ($comment =~ m/Target/);
	
	# extract Target sequence from $comment (if present)
	my $target_seq_id;
	if ($comment =~ m/Target \"([\w\d\.:-]+)\"/){
		$target_seq_id = $1;
	} 
	
	# only want to keep current feature if they overlap with any of the coding exons (or UTRs)
	# that we have already seen
	next unless (feature_overlaps_any_isoform($s_start, $s_stop));

	# need to keep track of how many different segments there are for this feature
	# segment = aligned region from parent feature
	my $count = 0;

	if(defined $features{$source}{$target_seq_id}){
		$count = @{$features{$source}{$target_seq_id}};
	}
	
	$features{$source}{$target_seq_id}[$count]{start} = $s_start;
	$features{$source}{$target_seq_id}[$count]{stop}  = $s_stop;
}
close($gff);




####################################################################################################
#
# 3rd loop
#
# now want to count which features are real hits to each isoform by looking at all
# 'segments' belonging to each feature
#
####################################################################################################

foreach my $gff_source (keys %features){	
#	print "Feature type: $gff_source\n";
	foreach my $feature (sort keys %{$features{$gff_source}}){
#		print "\t$feature\n";
		my $segment_count = @{$features{$gff_source}{$feature}};
#		print "\tSegments: $segment_count\n";
		# hash to store how many matches each feature has to each isoform
		my %isoform_matches;
		
		foreach my $isoform (keys %genes){
			$isoform_matches{$isoform} = 0;
#			print "\t\t$isoform\n";
			# now loop over all segments (aligned regions) for this feature
			for (my $i = 0; $i < $segment_count; $i++){
				my $s_start = $features{$gff_source}{$feature}[$i]{start};
				my $s_stop  = $features{$gff_source}{$feature}[$i]{stop};
				my $overlap = feature_overlaps_specific_isoform($s_start, $s_stop, $isoform);
#				print "\t\t$i) overlap = $overlap\n";
				$isoform_matches{$isoform} += $overlap;
			}
		}

		# what was the winning isoform?
		my @scores = sort {$b <=> $a} values %isoform_matches;
		my $top_score = $scores[0];
		
		foreach my $isoform (sort {$isoform_matches{$b} <=> $isoform_matches{$a}} keys %isoform_matches){
			my $score = $isoform_matches{$isoform};
#			print "\t$isoform\t$score";
			if ($score == $top_score){
				# this isoform is among the best fits to the feature data, 
				# so we can count this as a hit
				$genes{$isoform}{coding_exon}{results}{$gff_source}{unique_hits}++;								
#				print "*\n";
			} else{
#				print "\n";
			}
		}						
	}	
}


#print Dumper(\%genes);

########################
# final results
########################

# header row
print "Target_feature\t";
foreach my $isoform (sort keys %genes){	
	print "$isoform\t";
}
print "Ranking\n";


# hash to store counts of main ranking method
my %rank_counts;

# main results for each target feature
foreach my $key (sort keys %desired_GFF_features){
	my $target = $key;
	
	$target =~ tr /\t/:/;
	$target =~ s/:.*//;
	
	# don't want to consider coding_exon
	next if ($target eq 'Coding_transcript:coding_exon');

	print "$target\t";
	foreach my $isoform (sort keys %genes){		
		my $hits_per_Kbp = 0;	
		
		# if there were no hits, hits per Kbp will be zero. Otherwise, can calculate it...
		if (defined $genes{$isoform}{coding_exon}{results}{$target}{unique_hits}){
			my $hits       = $genes{$isoform}{coding_exon}{results}{$target}{unique_hits};
			my $sum_length = $genes{$isoform}{coding_exon}{total_length};
			$hits_per_Kbp  = sprintf("%.1f", ($hits / $sum_length) * 1000);		
		} 	

		
		# store details of hits per Kbp for this particular combination of isoform & target feature
		$genes{$isoform}{coding_exon}{results}{$target}{hits_per_Kbp} = $hits_per_Kbp;

		# sum across columns
		$genes{$isoform}{coding_exon}{final_score} += $hits_per_Kbp;
					
				
		print "$hits_per_Kbp\t";
	}

	# now calculate rankings
	my $rank = sort_hash($target);
	$rank_counts{$rank}++;
	print "$rank\n";
}



# print final sum scores
# will also want to grab scores of each isoform to work out extremes
my @absolute_scores;

print "Sum_scores:\t";


foreach my $isoform (sort keys %genes){		
	my $score = "$genes{$isoform}{coding_exon}{final_score}";
	print "$score\t";
	push(@absolute_scores, $score);
} 
print "\n\n";


# print counts of each ranking method
foreach my $ranking (sort {$rank_counts{$b} <=> $rank_counts{$a}} keys %rank_counts){		
	print "$ranking\t$rank_counts{$ranking}\n";
}

print "\n";

# now print out final differences between best and worst isoform
@absolute_scores = sort {$a <=> $b} @absolute_scores;

my ($max, $min) = ($absolute_scores[-1], $absolute_scores[0]); 

my $abs_diff = $max - $min;
my $rel_diff;


if($max == 0 and $min == 0){
	$rel_diff = 0;
	warn "POSSIBLE ERROR: $gene_id has 0 counts";
}else{
	$rel_diff = sprintf("%.2f", (1 - ($min / $max)));
}

print "Absolute difference between scores of best/worst isoform: $abs_diff\n";
print "Relative difference between scores of best/worst isoform: $rel_diff\n\n";


#print Dumper(\%genes);

exit;


#####################
#
# SUBROUTINEs
#
#####################

sub sort_hash{
 	my ($target)  = @_;

	# will build a new hash with scores of current feature
	my %hash;

	foreach my $isoform (sort keys %genes){		
		
		my $score = $genes{$isoform}{coding_exon}{results}{$target}{hits_per_Kbp};
		# trim isoform down to just the suffix
		$isoform =~ s/^.*\.[\d+]([a-z])/$1/;
		$hash{$isoform} = $score;
	}


	# need to keep track of when we are looking at tied values
	my $tied = 0;

	# need to build final output string	
	my $output_string = "";
	
	# need a way of sorting keys of tied values
	# to make sure we always see 'a(bc)de' and not sometimes 'a(cb)de'
	my @tied_values;
	
	# First sort main hash by values
	my @sorted_values = sort {$hash{$b} <=> $hash{$a}} (keys %hash);


	# now loop over array comparing pairs of adjacent values
	for (my $i = 0; $i < @sorted_values - 1; $i++){
	
		#print "$sorted_values[$i]\t$hash{$sorted_values[$i]}\n";			

		# if next value is identical, just add to array and set $tied flag
		if($hash{$sorted_values[$i]} == $hash{$sorted_values[$i+1]} and $hash{$sorted_values[$i]} != 0){
			push(@tied_values, $sorted_values[$i]);
			$tied = 1;			
		}
		elsif($hash{$sorted_values[$i]} == $hash{$sorted_values[$i+1]} and $hash{$sorted_values[$i]} == 0){
			push(@tied_values, $sorted_values[$i]);
			$tied = 2; #sets tied equal to two only when the two values are both 0's
			
		} else{
			# otherwise what we do depends on whether $tied is on or off
			if($tied == 1){
				push(@tied_values, $sorted_values[$i]);
				$output_string .= "(";
				$output_string .= join("", sort @tied_values);
				$output_string .= ")";
				@tied_values = ();
				$tied = 0;		
			}
			elsif($tied == 2 ) {
				push(@tied_values, $sorted_values[$i]);
				$output_string .= "[";
				$output_string .= join("", sort @tied_values);
				$output_string .= "]";
				@tied_values = ();
				$tied = 0;		
			}
			else{
				# if here, we have a regular not duplicate value that can be
				# added to output string
				$output_string .= "$sorted_values[$i]";
			}
		}
	}

	# Because above loop doesn't look at last value we have to treat
	# this separately 
	#print "$sorted_values[-1]\t$hash{$sorted_values[-1]}\n";
	
	if($tied == 0){
		$output_string .= "$sorted_values[-1]";
	} elsif($tied == 2) {
			push(@tied_values, $sorted_values[-1]);
			$output_string .= "[";
			$output_string .= join("", sort @tied_values);
			$output_string .= "]";
	} else{
			push(@tied_values, $sorted_values[-1]);
			$output_string .= "(";
			$output_string .= join("", sort @tied_values);
			$output_string .= ")";
	}

	#print "$output_string\n";	
	return $output_string;
}


# simple check to see coordinates of current feature (or part of feature)
# overlaps with previously seen coding exons coordinates (or UTRs)
sub feature_overlaps_any_isoform{
	my ($s_start, $s_stop) = @_;
	
	foreach my $isoform (keys %genes){
		foreach my $type (keys %{$genes{$isoform}}){
				for (my $i = 0; $i < @{$genes{$isoform}{$type}{exon_data}}; $i++){
				my $q_start  = $genes{$isoform}{$type}{exon_data}[$i]{start};
				my $q_stop   = $genes{$isoform}{$type}{exon_data}[$i]{stop};
	
				# skip if query and subject don't overlap
				if ($s_stop < $q_start or $s_start > $q_stop){
					next;
				} else{
					# we have an overlap, no need to look any further
					return(1);
				}
			}
		}
	}
	# no overlap
	return(0);
}

# a more complex check to see whether a segment from a feature overlaps any exons of a specific isoform
sub feature_overlaps_specific_isoform{
	my ($s_start, $s_stop, $isoform) = @_;
	my $type = 'coding_exon';
	
	for (my $i = 0; $i < @{$genes{$isoform}{$type}{exon_data}}; $i++){
		my $q_start  = $genes{$isoform}{$type}{exon_data}[$i]{start};
		my $q_stop   = $genes{$isoform}{$type}{exon_data}[$i]{stop};
		my $q_length = $q_stop - $q_start + 1;

		# skip if query and subject don't overlap
		if ($s_stop < $q_start or $s_start > $q_stop){
			next;
		} else{	
			# at this point we know that a specific segment of a feature (e.g. EST) overlaps an exon of 
			# a specific isoform. So don't need to continue checking against the remaining exons
			return(1);
		}
	}

	# no overlap
	return(0);
}


## simple hash to load details of all desired GFF source/feature combos from file into a hash
sub load_desired_GFF_details{
	open(my $params, "< $param_file") or die "Can't open $param_file";
	while(my $line = <$params>){
		chomp($line);
		next if ($line =~ m/^#/);
		$desired_GFF_features{$line} = 1;
	}
	close($params);
}


__END__

## now need to loop over all coding_exons/UTRs and then compare to all features to see what overlap
# Now we have a feature that we might be interested in, we need to compare
# coordinates to all coding exons and UTRs and look for overlaps
foreach my $isoform (keys %genes){	
	
	foreach my $type (keys %{$genes{$isoform}}){
		for (my $i = 0; $i < @{$genes{$isoform}{$type}{exon_data}}; $i++){
		
			my $q_start  = $genes{$isoform}{$type}{exon_data}[$i]{start};
			my $q_stop   = $genes{$isoform}{$type}{exon_data}[$i]{stop};
			my $q_length = $genes{$isoform}{$type}{exon_data}[$i]{'length'};
	
			# skip if query and subject don't overlap
			next if ($s_stop < $q_start or $s_start > $q_stop);
			# we have an overlap and can track the seq ID (if Target fields is present in $comment)
			if(defined $target_seq_id){
				# is this the first time we are seeing this target_seq_id?
				if(not defined ($matching_seq_IDs{"$isoform:$feature:$target_seq_id"})){
					$matching_seq_IDs{"$isoform:$feature:$target_seq_id"} = 1;
					$genes{$isoform}{$type}{unique_hits}++;								
				}
			}
			
			
			# if we haven't seen the current mix of query and subject we can
			# set up the sequence string to represent the current query feature
			# (coding exon or UTR). Use dashes to represent unmatched bases
			if (not defined $genes{$isoform}{$type}{exon_data}[$i]{matches}{$combined_key}{seq}){
				$genes{$isoform}{$type}{exon_data}[$i]{matches}{$combined_key}{seq} = '-' x $q_length;
				$genes{$isoform}{$type}{exon_data}[$i]{matches}{$combined_key}{hits} = 0;
			}

			
			# now convert chromosome coordinates into local match coordinates
			my ($m_start, $m_stop);
			
			if ($s_start < $q_start){
				$m_start = 0;
			} else {
				$m_start = $s_start - $q_start;
			}
	
			if ($s_stop > $q_stop){
				$m_stop = $q_length - 1;
			} else {
				$m_stop = $s_stop - $q_start;
			}
			my $m_length = $m_stop - $m_start + 1;

			
			# now substitute sequence of feature
			substr($genes{$isoform}{$type}{exon_data}[$i]{matches}{"$source:$feature"}{seq}, $m_start, $m_length) = ("o" x $m_length);	

			# can we increment the hit counter? I.e. was this a match? 
			# just look for any o's in the sequence string
			if($genes{$isoform}{$type}{exon_data}[$i]{matches}{$combined_key}{seq} =~ m/o/){
				$genes{$isoform}{$type}{exon_data}[$i]{matches}{$combined_key}{hits}++;
			}
			
#				my $hits = $genes{$isoform}{$type}{exon_data}[$i]{matches}{$combined_key}{hits};
#				print "\t$isoform\t$type\t$i\t$q_start\t$q_stop\t$q_length\t$m_start\t$m_stop\t$hits\t";
#				print "$genes{$isoform}{$type}{exon_data}[$i]{matches}{$combined_key}{seq}\n";		
		}
	}
}		
