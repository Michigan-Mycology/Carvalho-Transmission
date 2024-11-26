#!/usr/bin/perl

# This program takes a simple file in array format and makes long list file that can be read in R

# open input and output files
open (INPUTFILE, $ARGV[0]);
open (OUTPUTFILE, "> $ARGV[1]");

if (@ARGV < 3) {
    die "Improper number of arguments. Usage: perl convert_simple_to_long_Bd.pl infile outfile numstrains\n";
}

$numstrains=$ARGV[2] ;

print OUTPUTFILE "Position\tStrain\tGenotype\n";
$a=0;

$firstline = <INPUTFILE>;
print $firstline;
@names = split(/\t/, $firstline);
$names[$numstrains+1] =~ s/\s+$//;
for ($z=1; $z<$numstrains+1; $z++)
{
	$strainname[$z]=$names[$z];
	print "$strainname[$z]\t";
} 

$strainname[$numstrains] =~ s/\s+$//;
while (<INPUTFILE>)
{
	@info = split(/\t/, $_);
	$info[$numstrains] =~ s/\s+$//;
	$data[$a][0] = $info[0];
	for ($b=1; $b<$numstrains+1; $b++)
	{
		$data[$a][$b]=$info[$b];
	}
	$a++;
}

$numbases=$a;
print "numbase=$numbases";

for ($c=1; $c<$numstrains+1; $c++)
{
	for ($d=0; $d<$numbases; $d++)
	{
		print OUTPUTFILE "$data[$d][0]\t$strainname[$c]\t$data[$d][$c]\n";
	}
}
