#!/usr/bin/perl

open (DATAFILE, "BdGR.2022.snps.only.qualfiltered.pruned.no_missing.MICH.v4.txt"); 
open (OUTFILE, ">BdGR.2022.snps.only.qualfiltered.pruned.no_missing.no_ref.MICH.v4.txt"); 

$strainline = <DATAFILE>;
print OUTFILE $strainline;

chomp ($strainline);
$numstrains= $strainline;

$SNPline = <DATAFILE>;
print OUTFILE $SNPline;
chomp ($SNPline);
$numsnps = $SNPline;

$contigline = <DATAFILE>;
print OUTFILE $contigline;
$listofnames = <DATAFILE>;
print OUTFILE $listofnames;

while ($snp < $numsnps)
{
	$nextline = <DATAFILE>;
	chomp ($nextline);
	@genos = split(/\t/, $nextline);
# This following line removes the whitespace at the end of the line.
	for ($a=0;$a<$numstrains+3;$a++)
	{
		unless ($a == 2) {if ($a==$numstrains+2) {print OUTFILE "$genos[$a]\n";} else {print OUTFILE "$genos[$a]\t";}}
	}
	$snp++;
}



