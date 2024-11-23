#!/usr/bin/perl

use List::Util qw(sum);

# This program takes a vcf file and makes a file in Michigan format
# It scores samples as N that have less than a certain coverage

# open input and output files
open (INPUTFILE, $ARGV[0]);
open (OUTPUTALLELES, "> $ARGV[1]");
open (FILTERFILE, "> $ARGV[2]");

# Pass number of strains as the third argument after removing any line breaks
$numstrains=$ARGV[3];
$numstrains =~ s/\s+$//;

if (@ARGV < 3) {
    die "Improper number of arguments. Usage: perl convert_vcf_MICH_format_filter_depth.pl inputvcffile outMICH_format_file numstrains\n";
}

print "numstrains=$numstrains\n";

# Cols is the number of columns before the real data that are variable between vcf files
$cols = 9;

# Removed mindepth and just used lowdepth

# maxdepth is the maximum fold over average depth for each strain to remove.
$maxdepth = 4;
# lowdepth is the same as max except on the lower end
$lowdepth = 4;
# GQmin is the value if Genotype is lower than this quality is filtered out
$GQmin = 50;

# This is the string for SNPs that are filtered out, currently not used
$filter='hard_filter_1';

# Rolling sum of depths

my $depths;
@depths = (0) x ($numstrains+$cols);

# Rolling sum of positions with depths

my $sumpos;
@sumpos = (0) x ($numstrains+$cols);

# First loop through to get the average depth

while (<INPUTFILE>)
{
	unless ($_ =~ /##/)
	{
# This is the stuff that prints the strain names
		if ($_ =~ /CHROM/)
		{
			@genos = split(/\t/, $_);
# This following line removes the whitespace at the end of the line.
			$genos[$numstrains+$cols-1] =~ s/\s+$//;
			for ($a=$cols;$a<($numstrains+$cols);$a++)
			{	
				$names[$a] = $genos[$a]; 
			}
		}	
		else
		{
#			print "\nthis is the main area for getting the depths
			@genos = split(/\t/, $_);
# This following line removes the whitespace at the end of the line.
			$genos[$numstrains+$cols] =~ s/\s+$//;
# Look to see SNPs that were filtered by GATK VariantFiltration
# Or if the second allele has more than 1 variant in which case there is a ','
			unless (($genos[6] =~ /$filter/)|($genos[4] =~ /\,/))
			{
				for ($b=$cols;$b<($numstrains+$cols);$b++)
				{
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

# Then loop through again and print

open (INPUTFILE, $ARGV[0]);

print OUTPUTALLELES "Contig\tPosition\t ";

while (<INPUTFILE>)
{
	unless ($_ =~ /##/)
	{
		if ($_ =~ /CHROM/)
		{
			@genos = split(/\t/, $_);
#	print "@genos\n";
# This following line removes the whitespace at the end of the line.
			$genos[$numstrains+$cols] =~ s/\s+$//;
			for ($a=$cols;$a<($numstrains+$cols);$a++)
			{	
				print OUTPUTALLELES "$genos[$a]\t"
			}
#		print "\n";
		print OUTPUTALLELES "\n";
		}	
		else
		{
#			print "\nthis is the main printing area";
			@genos = split(/\t/, $_);
# This following line removes the whitespace at the end of the line.
			$genos[$numstrains+$cols] =~ s/\s+$//;
# This prints the Contig number and location if the SNP passed
#			$genos[0] =~ s/bden_JEL423_supercont1.//g;
			$genos[0] =~ s/Supercontig_1.//g;
#			$genos[0] =~ s/bden_JEL423_MT/mt68/g;
			if ($genos[4] =~ /\,/) {print FILTERFILE "Multi-allele\n";}
# Look to see SNPs that were filtered by GATK VariantFiltration
# Or if the second allele has more than 1 variant in which case there is a ','
			unless (($genos[6] =~ /$filter/)|($genos[4] =~ /\,/))
			{
				print OUTPUTALLELES "$genos[0]\t$genos[1]\t";
				$refallele=$genos[3];
				print OUTPUTALLELES "$refallele\t";
				$altallele=$genos[4];
				for ($b=$cols;$b<($numstrains+$cols);$b++)
				{
					@indgenos = split(/:/, $genos[$b]);
#					print "b=$b";
					if ($indgenos[0] =~ /\./)
					{
# Genotype uncertain
						print OUTPUTALLELES "N";
					}
					elsif (($indgenos[2] > ($maxdepth * $meandepth[$b])) | ($indgenos[2] < ($meandepth[$b] / $lowdepth)) | ($indgenos[3] < $GQmin))
					{
						if ($indgenos[2] > ($maxdepth * $meandepth[$b])) {print FILTERFILE "Maxdepth\n";}
						if ($indgenos[2] < ($meandepth[$b] / $lowdepth)) {print FILTERFILE "Lowdepth\n";}
						if ($indgenos[3] < $GQmin) {print FILTERFILE "GQmin\n";}
						$temp1 = $maxdepth * $meandepth[$b];
						$temp2 = $meandepth[$b]/$lowdepth;
						print "depth=$indgenos[2]\t GQ=$indgenos[3] maxdepth*meandepth[$names[$b]]=$temp1 meandepth[$names[$b]]/lowdepth=$temp2 GQmin=$GQmin\n";
# Failed QC
						print OUTPUTALLELES "N";
					}
					elsif ($indgenos[0] =~ /2/)
# Multiallelic, print NA (note shouldn't get to this point as should be filtered earlier)
					{
						print OUTPUTALLELES "N";
					}				
					else {printgeno($refallele,$altallele,$indgenos[0]);}
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
sub printgeno {
my ($base1, $base2, $genotype) = @_;

if (($genotype =~ /0\/0/) | ($genotype eq '0|0'))
{
	print OUTPUTALLELES "$base1";
}
elsif (($genotype =~ /1\/1/) | ($genotype eq '1|1'))
{
	print OUTPUTALLELES "$base2";
}
else 
{
	if (($base1 eq 'A')&&($base2 eq 'C')) {print OUTPUTALLELES "M";}
	if (($base1 eq 'A')&&($base2 eq 'G')) {print OUTPUTALLELES "R";}
	if (($base1 eq 'A')&&($base2 eq 'T')) {print OUTPUTALLELES "W";}
	if (($base1 eq 'C')&&($base2 eq 'A')) {print OUTPUTALLELES "M";}
	if (($base1 eq 'C')&&($base2 eq 'G')) {print OUTPUTALLELES "S";}
	if (($base1 eq 'C')&&($base2 eq 'T')) {print OUTPUTALLELES "Y";}
	if (($base1 eq 'G')&&($base2 eq 'A')) {print OUTPUTALLELES "R";}
	if (($base1 eq 'G')&&($base2 eq 'C')) {print OUTPUTALLELES "S";}
	if (($base1 eq 'G')&&($base2 eq 'T')) {print OUTPUTALLELES "K";}
	if (($base1 eq 'T')&&($base2 eq 'A')) {print OUTPUTALLELES "W";}
	if (($base1 eq 'T')&&($base2 eq 'C')) {print OUTPUTALLELES "Y";}
	if (($base1 eq 'T')&&($base2 eq 'G')) {print OUTPUTALLELES "K";}
}
}

sub mean {
    return sum(@_)/@_;
}

