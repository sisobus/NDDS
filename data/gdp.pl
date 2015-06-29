#!/usr/bin/perl

=pod

=head1 SYNOPSIS

This perl scripts generates N data points such that they have C inherent clusters.
Each data point has D dimensions and each dimension has Max size of A. The values N, C, D, A are passed as command line arguments.
You can also specify percentage of outliers to control facrtion of outliers.

    e.g. To generate 100 points (each having 10 dimensions of alphabet size  4) having 5 clusters with 7%outliers use:
        generate_data_points.pl -outlier_percentage 7 100 5 10 4 

    Use -perldoc option to see detailed help message.

=head1 OPTIONS

=over

=item -outlier_percentage

This option is used to specify percentage outliers. Default is 3%

=item -help

Displays a brief help message and exits

=item -perldoc

Displays a detailed documentation.

=back

=cut

use strict;
use warnings 'all';

use Getopt::Long;
use Pod::Usage;

my $percent_outliers = 3;

GetOptions('outlier_percentage=s' => \$percent_outliers, 'help' => \my $help, 'perldoc' => \my $perldoc);

if($perldoc)
{
    pod2usage(-verbose => 2);
}

if($help)
{
    pod2usage(-help => 1);
}

if (scalar(@ARGV) < 4)
{
	print("The script expects 4 arguments\n");
	print("Number of data points\n");
	print("Number of clusters\n");
	print("Number of dimenstions\n");
	print("Alphabet size\n");
	exit 1;
}

my $num_data_pt = $ARGV[0];
my $num_clusters = $ARGV[1];
my $num_dim = $ARGV[2];
# my $size_of_subset = $ARGV[3];            # upper bound for seed
# my $no_of_overlap = 0;
# my $dinominator = $ARGV[4];
my $alphabet_size = $ARGV[3];

#Preliminary validation
if($num_data_pt > ($alphabet_size**$num_dim))
{
	print("Number of data points is larger than the number of possible combinations\n");
	exit(2);
}

# print values of all the parameters
# print("$num_data_pt,$num_clusters,$num_dim,$alphabet_size;\n");
# Seed the random number generator
srand(time);

# An array of seeds for the data points.
my @seeds;
my @seed_ids;
# Generate seeds
my %used_dim;
for(my $i=0;$i<$num_clusters;$i++)
{
    for(my $j=0;$j<$num_dim;$j++)
    {
        my $el = 65 + int(rand($alphabet_size));
        while(defined(${used_dim{$j}}{$el}) && $num_clusters<=$alphabet_size)
        {
        	$el = 65 + int(rand($alphabet_size));
        }
        ${used_dim{$j}}{$el} = 1;
        if($el > 90)
        {
            $el += 6;
        }
        ${seeds[$i]}[$j] = sprintf("%c", $el);

        # ID of any data point has following format : Cluster#-Point#
        # All the points in the same cluster will have same Cluster#
    }
    $seed_ids[$i] = sprintf("%04d-%06d", $i, $i);
}

#my $datafile = "./data.txt";
my $datafile = $ARGV[4];
open(DFILE,">$datafile");

#Print a comment informing reader about parameters
print DFILE "$num_data_pt $num_dim $num_clusters \n";

for(my $i=0;$i<$num_clusters;$i++)
{
    print DFILE "$seed_ids[$i]:", join(",", @{$seeds[$i]}), "\n";
}

my $points_to_be_created = $num_data_pt - $num_clusters;
my %ids;
# Percentage of outliers
$percent_outliers /= 100;
my $threshold = $percent_outliers* $num_data_pt;

while($points_to_be_created > 0)
{
    my $next_seed = int(rand($num_clusters));
    my @next_data_point;
    my $id;
    for(my $i=0;$i<$num_dim;$i++)
    {
        $next_data_point[$i] = ${seeds[$next_seed]}[$i];
    }
    # Any point in the cluster can be different from seed point in at the most
    # 50% of dimensions.
    my $num_dim_to_change;
    if($points_to_be_created > $threshold)
    {
        $num_dim_to_change = 1 + int(rand(0.25*$num_dim));
        $id = sprintf("%04d-%06d", $next_seed, ($num_data_pt - $points_to_be_created))
    }
    else
    {
        # This is an outlier
        $num_dim_to_change = 1 + 0.35 * $num_dim + int(rand(0.2*$num_dim));
        $id = sprintf("%04d-#%05d", $next_seed, ($num_data_pt - $points_to_be_created))
    }
    for(my $i=0;$i<$num_dim_to_change;$i++)
    {
        # Randomly select a dimension and change it
        my $dim_to_change = int(rand($num_dim));
        my $val = 65 + int(rand($alphabet_size));
        if($val > 90)
        {
            $val += 6;
        }
        my $new_value = sprintf("%c", $val);
        # Make sure the dimension DOES get changed
        while($new_value eq $next_data_point[$dim_to_change])
        {
            $val = 65 + int(rand($alphabet_size));
            if($val > 90)
            {
                $val += 6;
            }
            $new_value = sprintf("%c", $val);
        }
        $next_data_point[$dim_to_change] = $new_value;
    }
    my $vector = join(",", @next_data_point);
    my $flag = 1;
    for(my $i=0;$i < scalar(@seeds);$i++)
    {
        my $seed = join(",", @{$seeds[$i]});
        if($vector eq $seed)
        {
            $flag = 0;
            last;
        }
    }
    if($flag && !defined($ids{$vector}))
    {
        #This point is new (i.e there is no duplicate)
        $ids{$vector} = $id;
        $points_to_be_created--;
    }
}
foreach my $key(keys(%ids))
{
    print DFILE "$ids{$key}:$key\n";
}

close(DFILE);
