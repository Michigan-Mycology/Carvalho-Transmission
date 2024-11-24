#!/usr/bin/perl

use List::Util qw(sum);

# This program takes a vcf file and calculates the allele frequencies for each strain for each chromosome

# open input and output files
open (INPUTFILE, $ARGV[0]);
open (OUTPUTALLELES, "> $ARGV[1]");

# Pass number of strains as the third argument after removing any line breaks
$numstrains=$ARGV[2];
$numstrains =~ s/\s+$//;

if (@ARGV < 2) {
    die "Improper number of arguments. Usage: perl extract_allele_frequencies.pl inputvcffile output_allele_file numstrains\n";
}

print "numstrains=$numstrains\n";

# Cols is the number of columns before the real data that are variable between vcf files
$cols = 9;

# maxdepth is the maximum fold over average depth for each strain to remove.
$maxdepth = 4;
# lowdepth is the same as max except on the lower end
$lowdepth = 4;
# GQmin is the value if Genotype is lower than this quality is filtered out
$GQmin = 50;

# Array of rolling sum of depths

my $depths;
@depths = (0) x ($numstrains+$cols);

# Array of rolling sum of positions with depths

my $sumpos;
@sumpos = (0) x ($numstrains+$cols);

# First loop through to get the average depth for determining what cutoff level to use for depths

while (<INPUTFILE>)
{
	unless ($_ =~ /##/)
	{
# This prints the strain names
		if ($_ =~ /CHROM/)
		{
			@genos = split(/\t/, $_);
			$genos[$numstrains+$cols-1] =~ s/\s+$//;
			for ($a=$cols;$a<($numstrains+$cols);$a++)
			{	
				$names[$a] = $genos[$a]; 
			}
		}	
		else
		{
			@genos = split(/\t/, $_);
			$genos[$numstrains+$cols] =~ s/\s+$//;
# Look to see SNPs that were filtered by GATK VariantFiltration
# Or if the second allele has more than 1 variant in which case there is a ','
			unless (($genos[6] =~ /$filter/)|($genos[4] =~ /\,/))
			{
				for ($b=$cols;$b<($numstrains+$cols);$b++)
				{
# Get depths as 3rd element of array
					@indgenos = split(/:/, $genos[$b]);
					unless ($indgenos[0] =~ /\./) 
					{
						$depths[$b] += $indgenos[2];
						$sumpos[$b]++;
					}
				}
			}
		}
	}
}

# Calculate depth means
for ($a=$cols;$a<($numstrains+$cols);$a++)
{	
	$meandepth[$a] = $depths[$a]/$sumpos[$a];
	print "meandepth[$names[$a]] = $meandepth[$a]\n";
}

close INPUTFILE;

# Then loop through again and print if passes QC

open (INPUTFILE, $ARGV[0]);

print OUTPUTALLELES "Contig\tPosition\t ";

while (<INPUTFILE>)
{
	unless ($_ =~ /##/)
	{
		if ($_ =~ /CHROM/)
		{
			@genos = split(/\t/, $_);
			$genos[$numstrains+$cols-1] =~ s/\s+$//;
			for ($a=$cols;$a<($numstrains+$cols);$a++)
			{	
# printing strain names
				print OUTPUTALLELES "$genos[$a]\t"
			}
		print OUTPUTALLELES "\n";
		}	
		else
		{
			@genos = split(/\t/, $_);
			$genos[$numstrains+$cols] =~ s/\s+$//;
# This prints the Contig number and location if the SNP passed
			$genos[0] =~ s/Supercontig_1.//g;
# Look to see SNPs that were filtered by GATK VariantFiltration
# Or if the second allele has more than 1 variant in which case there is a ','
			unless (($genos[6] =~ /$filter/)|($genos[4] =~ /\,/))
			{
				print OUTPUTALLELES "$genos[0]\t$genos[1]\t";
				$refallele=$genos[3];
				$altallele=$genos[4];
				for ($b=$cols;$b<($numstrains+$cols);$b++)
				{
					@indgenos = split(/:/, $genos[$b]);
					if ($indgenos[0] =~ /\./)
					{
# detected missing genotype, print NA;
						print OUTPUTALLELES "NA";
					}
					elsif (($indgenos[2] > ($maxdepth * $meandepth[$b])) | ($indgenos[2] < ($meandepth[$b] / $lowdepth)) | ($indgenos[3] < $GQmin))
					{
# Failed QC or depth, print NA
						$temp1 = $maxdepth * $meandepth[$b];
						$temp2 = $meandepth[$b]/$lowdepth;
						print "depth=$indgenos[2]\t qual=$indgenos[3] maxdepth*meandepth[$names[$b]]=$temp1 meandepth[$names[$b]]/lowdepth=$temp2 GQmin=$GQmin\n";
						print OUTPUTALLELES "NA";
					}
					elsif (($indgenos[0] eq '1|1') | ($indgenos[2] =~ /1\/1/))
# Homozygous, print NA
					{
						print OUTPUTALLELES "NA";
					}
					elsif ($indgenos[0] =~ /2/)
# Multiallelic, print NA (note shouldn't get to this point as should be filtered earlier)
					{
						print OUTPUTALLELES "NA";
					}				
					else 
# Should be either 0/1 or 0|1
					{
# Get depths of each allele
						@alleles = split(/,/, $indgenos[1]);
						printallelefreq($indgenos[0],$alleles[0],$alleles[1]);
					}
					unless ($b==$numstrains+$cols-1)
					{
						print OUTPUTALLELES "\t";
					}
				}
				print OUTPUTALLELES "\n";
			}
		}
	}
}
sub printallelefreq {
my ($genotype, $all1count, $all2count) = @_;

if (($genotype =~ /0\/1/) | ($genotype eq '0|1'))
{
	if ($all2count == 0) 
	{
		print OUTPUTALLELES "NA";
	}
	else
	{
		$freq = $all1count/($all1count + $all2count);
		printf OUTPUTALLELES "%4f", $freq;
	}
}

else
{
	print OUTPUTALLELES "NA";
}
}

# Not used
sub mean {
    return sum(@_)/@_;
}

