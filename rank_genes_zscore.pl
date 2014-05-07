#!/usr/bin/perl
#
# rank_genes_zscore.pl
# calculates a z score for each gene based on various factors (does same thing as compare_gene_rankings.R)
#
# ranks genes by z score
# ranks based on cds length, exon length, intron length, codon usage, donor splice site usage, acceptor splice site usage, intron usage

use strict; use warnings;
use File::Slurp;
use Statistics::Zscore;
use Data::Dumper;

die "Enter a gene to get the score breakdown for or enter 'none'" unless @ARGV == 1;
my $gene = $ARGV[0];

# reads in the files
my @cds_length = read_file("cds_lengths.txt");
my @exon_length = read_file("exon_lengths.txt");
my @intron_length = read_file("intron_lengths.txt");
my @codon_usage = read_file("gene_codon_usage_relative_entropies.txt");
my @donor = read_file("splice_site_donor_usage.txt");
my @acceptor = read_file("splice_site_acceptor_usage.txt");
my @intron_usage = read_file("intron_usage.txt");

# adds the z scores to the respective arrays
@cds_length = add_zscore_to_array(@cds_length);
@exon_length = add_zscore_to_array(@exon_length);
@intron_length = add_zscore_to_array(@intron_length);
@codon_usage = add_zscore_to_array(@codon_usage);
@donor = add_zscore_to_array(@donor);
@acceptor = add_zscore_to_array(@acceptor);
@intron_usage = add_zscore_to_array(@intron_usage);

# print join("\n", @donor);

# adds the z scores together and puts the info for the genes into a hash
my %master;
add_zscores_together("cds", @cds_length);
add_zscores_together("exon_length", @exon_length);
add_zscores_together("intron_length", @intron_length);
add_zscores_together("codon", @codon_usage);
add_zscores_together("donor", @donor);
add_zscores_together("acceptor", @acceptor);
add_zscores_together("intron_usage", @intron_usage);

# calculates an average z score for each gene, and prints that result
# takes absolute values of scores b/c we don't want + and - zscores to even each other out
for my $key (keys %master){
	my $final_score = 0;
	my $count = 0;
	
	if (exists $master{$key}{cds}){
		$final_score += abs($master{$key}{cds}{zscore});
		$count++;
	}
	if (exists $master{$key}{exon_length}){
		if (exists $master{$key}{exon_length}{avg_zscore}){
			$final_score += abs($master{$key}{exon_length}{avg_zscore});
			$count++;
		}else{
			$final_score += abs($master{$key}{exon_length}{zscore});
			$count++;
		}
	}
	if (exists $master{$key}{intron_length}){
		if (exists $master{$key}{intron_length}{avg_zscore}){
			$final_score += abs($master{$key}{intron_length}{avg_zscore});
			$count++;
		}else{
			$final_score += abs($master{$key}{intron_length}{zscore});
			$count++;
		}
	}
	if (exists $master{$key}{codon}){
		$final_score += abs($master{$key}{codon}{zscore});
		$count++;
	}
	if (exists $master{$key}{donor}){
		$final_score += abs($master{$key}{donor}{zscore});
		$count++;
	}
	if (exists $master{$key}{acceptor}){
		$final_score += abs($master{$key}{acceptor}{zscore});
		$count++;
	}
	if (exists $master{$key}{intron_usage}){
		$final_score += abs($master{$key}{intron_usage}{zscore});
		$count++;
	}
	$master{$key}{final_score} = $final_score / $count;
	
	# prints result
	if ($gene eq "none"){
		print"$master{$key}{final_score}, $master{$key}{chromosome}, $key\n";
	}
}

if ($gene ne "none"){
	for my $key (keys %master){
		if ($key =~ m/$gene/){
			if (exists $master{$key}{cds}){
				print"$key, cds, $master{$key}{cds}{zscore}\n";
			}
			if (exists $master{$key}{exon}){
				print"$key, exon, $master{$key}{exon}{zscore}, $master{$key}{exon}{count}\n";
			}
			if (exists $master{$key}{intron}){
				print"$key, intron, $master{$key}{intron}{zscore}, $master{$key}{intron}{count}\n";
			}
			if (exists $master{$key}{codon}){
				print"$key, codon, $master{$key}{codon}{zscore}\n";
			}
			if (exists $master{$key}{donor}){
				print"$key, donor, $master{$key}{donor}{zscore}\n";
			}
			if (exists $master{$key}{acceptor}){
				print"$key, acceptor, $master{$key}{acceptor}{zscore}\n";
			}
			if (exists $master{$key}{intron_usage}){
				print"$key, intron_usage, $master{$key}{intron_usage}{zscore}\n";
			}
		}
	}
}

# print Dumper(\%master);
#################
###SUBROUTINES###
#################
sub add_zscore_to_array{
	my @array = @_;
	
	# calculate z score
	my $z = Statistics::Zscore->new;
	my $zscore_ref = $z->standardize( \@array);
	my @zscore = @$zscore_ref;

	# adds the z score to the array for each entry
	my $i = 0;
	for my $line (@array){
		chomp $line;
		$array[$i] = $line . ", " . $zscore[$i];
		$i++;
	}
	return(@array);
}

sub add_zscores_together{
	my ($type, @array) = @_;
	my $i = 0;
	my ($prev_gene, $prev_zscore);
	for my $line (@array){
		chomp $line;
		my ($length, $chromosome, $gene, $zscore) = split(", ", $line);
	# 	print"$line\n";
		if ($i == 0){
			$master{$gene}{chromosome} = $chromosome;
			$master{$gene}{$type}{zscore} = $zscore;
			$master{$gene}{$type}{count} = 1;
		}else{
			if ($gene ne $prev_gene){
				$master{$gene}{chromosome} = $chromosome;
				$master{$gene}{$type}{zscore} = $zscore;
				$master{$gene}{$type}{count} = 1;
			}elsif ($gene eq $prev_gene){
				$master{$gene}{chromosome} = $chromosome;
				$master{$gene}{$type}{zscore} += $zscore;
				$master{$gene}{$type}{count}++
			}
		}
		$prev_gene = $gene;
		$prev_zscore = $zscore;
		$i++;
	}
}