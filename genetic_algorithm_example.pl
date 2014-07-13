#!/usr/bin/perl
#
# genetic_algorithm_example.pl
#
# example of a genetic algorithm using shapes and their areas

use strict; use warnings;
use Math::Trig;


die "<generations> <max population> <% to keep> <macro mutation % chance> <micro mutation % chance> <% to keep safe from mutation> \n" unless (@ARGV == 6);
my ($gen, $max_population, $percent_to_keep, $macro_mut, $micro_mut, $no_mutation) = @ARGV;
 
die "Invalid value for number of generations\n" unless ($gen > 0);
die "Invalid value for max population\n" unless ($max_population > 0);
die "Invalid value for percent to keep\n" unless ($percent_to_keep > 0 and $percent_to_keep < 100);
die "Invalid value for macro mutation percent chance\n" unless ($macro_mut > 0 and $macro_mut < 100);
die "Invalid value for micro mutation percent chance\n" unless ($micro_mut > 0 and $micro_mut < 100);
die "Invalid value for percent to keep safe from mutation\n" unless ($no_mutation < (100-$percent_to_keep)/2);


# creates the population of shapes and their lengths
my %pop;

create_population();

for (my $g = 1; $g <= $gen; $g++){
    print"\n\nGeneration $g\n--------------\n";

    # delete mutants (except in 1st generation)
    unless ($g == 1){
	print "Delete mutants\n";
	delete_mutants();
    }

#     # calculates fitness for non mutants
#     print "Calculate fitness\n";
#     calc_fitness();

    # individuals sorted by fitness
    # deletes those that are not in the top x%
    print "Deletes unfit\n";
    delete_unfit();
    
    # repopulates the population
    print "Repopulate population\n";
    create_population();

#     calc_fitness();
    
    
    # mutation
    # the top y% of population safe from being mutated
    print "Mutate population\n";
    mutates_population();

}



###SUBROUTINES###

# chooses a shape
sub choose_shape{
    my @shapes = ("circle", "triangle", "square", "hexagon");
    my $rand = int(rand(@shapes));
    my $rand_shape = $shapes[$rand];
    return($rand_shape);
}

# chooses a length for a shape
sub choose_length{
    my $length = int(rand(10)) + 1;
    return($length);
}

# calculates the area for a shape
# need shape and length as inputs
sub calc_area{
    my ($shape, $length) = @_;
    my $area;
    if ($shape eq "circle"){
	$area = pi * $length * $length;
    }elsif ($shape eq "triangle"){
	$area = $length * $length * sqrt(3)/2 * 1/2;
    }elsif ($shape eq "square"){
	$area = $length * $length;
    }elsif ($shape eq "hexagon"){
	$area = $length * $length * sqrt(3)/2 * 1/2 * 6;
    }
    return($area);
}


# creates and/or repopulates the population
sub create_population{
    for (my $i = 1; $i <= $max_population; $i++){
	next unless (not exists $pop{$i}{inner_shape});
	my $shape1 = choose_shape();
	my $shape2 = choose_shape();
	redo if ($shape1 eq $shape2);
	
	my $length1 = choose_length();
	my $length2 = choose_length();
	my $inner_area = calc_area($shape1, $length1);
	my $outer_area = calc_area($shape2, $length2);
	redo if ($inner_area > $outer_area);
	my $fitness = $outer_area - $inner_area;
	
	# if we are here, everything is good
	$pop{$i}{inner_shape} = $shape1;
	$pop{$i}{outer_shape} = $shape2;
	$pop{$i}{inner_length} = $length1;
	$pop{$i}{outer_length} = $length1;
	$pop{$i}{inner_area} = $inner_area;
	$pop{$i}{outer_area} = $outer_area;
	$pop{$i}{fitness} = sprintf("%.2f", $fitness);
    }

    my $n = keys %pop;
}


# finds mutants. mutants if the 2 shapes are the same, or if inner area > outer area
sub delete_mutants{
    my $count = 0;
    for my $key (keys %pop){
	if ($pop{$key}{inner_shape} eq $pop{$key}{outer_shape}){
	    delete $pop{$key};
	    $count++;
	}elsif($pop{$key}{inner_area} > $pop{$key}{outer_area}){
	    delete $pop{$key};
	    $count++;
	}
    }
    print "\tDeleted $count mutants\n";
}

# # calculates fitness
# sub calc_fitness{
#     for my $key (keys %pop){
# 	$pop{$key}{fitness} = sprintf("%.2f", $pop{$key}{outer_area} - $pop{$key}{inner_area});
#     }
# }


# individuals sorted by fitness
# deletes those that are not in the top x%
sub delete_unfit{
    my $count = 1;
    my $keep = int($percent_to_keep * 0.01 * $max_population);
    print "\tKeeping the $keep fittest individuals out of $max_population\n";

    for my $key (sort {$pop{$a}{fitness} <=> $pop{$b}{fitness}} keys %pop){
	if ($count <= $keep){
	    print "\t$key) Fitness = $pop{$key}{fitness} ($pop{$key}{inner_shape} $pop{$key}{inner_length} $pop{$key}{outer_shape} $pop{$key}{outer_length})\n";
	}else{
	    delete $pop{$key};
	}
	$count++;
    }
}


# mutates population
sub mutates_population{
    my $safe = int($no_mutation * 0.01 * $max_population);
    print "Keep the top $safe safe from mutation\n";
    my $count = 0;
    
    for my $key (sort {$pop{$a}{fitness} <=> $pop{$b}{fitness}} keys %pop){
	$count++;
	
	# skips mutating the top % percent of the population
	next if ($count <= $safe);
	
	# is there a macro mutation for this individual's inner/outer shape?
	$pop{$key}{inner_shape} = choose_shape() if (rand(100) <= $macro_mut);
	$pop{$key}{outer_shape} = choose_shape() if (rand(100) <= $macro_mut);

	# is there a macro mutation for this individual's inner/outer length?
	$pop{$key}{inner_length} = choose_length() if (rand(100) <= $macro_mut);
	$pop{$key}{outer_length} = choose_length() if (rand(100) <= $macro_mut);


	# is there a micro mutation for this individual's inner/outer length?
	$pop{$key}{inner_length} += 0.1  if (rand(100) <= $micro_mut);
	$pop{$key}{inner_length} -= 0.1  if (rand(100) <= $micro_mut);
	$pop{$key}{inner_length} += 0.01 if (rand(100) <= $micro_mut);
	$pop{$key}{inner_length} -= 0.01 if (rand(100) <= $micro_mut);
	$pop{$key}{outer_length} += 0.1  if (rand(100) <= $micro_mut);
	$pop{$key}{outer_length} -= 0.1  if (rand(100) <= $micro_mut);
	$pop{$key}{outer_length} += 0.01 if (rand(100) <= $micro_mut);
	$pop{$key}{outer_length} -= 0.01 if (rand(100) <= $micro_mut);
    }
}



