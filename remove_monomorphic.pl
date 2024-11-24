#!/usr/bin/perl

# This program removes non-polymorphic loci

open (DATAFILE, "$ARGV[0]"); 
open (POLYSONLY, ">$ARGV[1]");

# First lines is the strain names
$strainline = <DATAFILE>;
print POLYSONLY "$strainline";

$numstrains=27;
$minpres=$numstrains/10;
#$minpres=0;

while (<DATAFILE>)
{
	@genos = split(/\t/, $_);
	$genos[$numstrains+2] =~ s/\s+$//;
	$ref = $genos[2];
	$mono='TRUE';
	$sumpres = 0;
	for ($a=0;$a<$numstrains;$a++)
	{
		if (($ref ne $genos[$a+2]) && ($genos[$a+2] ne 'N'))
		{
			$mono='FALSE';
			$sumpres++;
		}
	}
	if (($mono eq 'FALSE') && ($sumpres > $minpres))
	{
		print POLYSONLY $_;
	}
}

exit;