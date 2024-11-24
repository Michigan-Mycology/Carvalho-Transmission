#!/usr/bin/perl

# This line is needed to include the min and max functions
use List::Util qw(min max);

# This program calculates Ho for a set of strains "results_het_strain.txt"
# Identifies mutations that are unique to one strain "unique_mutations.txt"
# Finds putative unique mutations in regions of the genome that have undergone LOH "recent_mutations.txt"
# Prints the average heterozygosity in a sliding window for each strain and each contig "sliding_het.txt"
# Identifies sites that violate the infinite sites model "results_violators.txt"
# Calculates Inbreeding coefficients across the genome "results_het_genome.txt"
# Computes genetic distances among strains using the hetequal model "results_distance.txt"
# And tests for linkage disequilibrium within and across contigs "results_LD.txt"

# This code will need to be modified if there is missing data

open (DATAFILE, "BdGR.2022.snps.only.qualfiltered.pruned.no_missing.no_ref.MICH.v4.txt"); 
#open (UNIQUEMUTS, ">Bd_exp_evol.10X.unique_mutations.txt"); 
#open (RECENTMUTS, ">Bd_exp_evol.10X.recent_mutations.txt"); 
open (SLIDEHET, ">BdGR.2022.snps.only.qualfiltered.pruned.no_missing.no_ref.MICH.het_cumul.txt"); 
#open (OUTPUTFILEVIOL, ">Bd_28.ACGT.10X.results_violators.txt");
open (OUTPUTFILEHETSTRAIN, ">BdGR.2022.snps.only.qualfiltered.pruned.no_missing.no_ref.MICH.hetstrain.txt");
open (DISTANCES, ">BdGR.2022.snps.only.qualfiltered.pruned.no_missing.no_ref.MICH.distances.txt");
open (OUTPUTFILEHETGENOME, ">BdGR.2022.snps.only.qualfiltered.pruned.no_missing.no_ref.MICH.het_genome.txt");
#open (OUTPUTFILEHETSITES, ">Bd_exp_evol.10X.het_sites.txt");
#open (OUTPUTFILELD, ">Bd_28.con1.rtarm.10X.results_LD.txt");

# First and second and third ines are the number of strains followed by number of SNPs followed by number of contigs.
$strainline = <DATAFILE>;
chomp ($strainline);
$numstrains= $strainline;

$SNPline = <DATAFILE>;
chomp ($SNPline);
$numsnps = $SNPline;

$contigline = <DATAFILE>;
chomp ($contigline);
$numcontigs = $contigline;

# This value is the maximum reported supercontig cutoff due to size.
$threshold=15;

# These values are the running size of the contigs
# These values are the running size of the contigs
$contig_length[1]=0;
$contig_length[2]=4440149;
$contig_length[3]=$contig_length[2]+2313122;
$contig_length[4]=$contig_length[3]+1829408;
$contig_length[5]=$contig_length[4]+1803316;
$contig_length[6]=$contig_length[5]+1707251;
$contig_length[7]=$contig_length[6]+1545501;
$contig_length[8]=$contig_length[7]+1398854;
$contig_length[9]=$contig_length[8]+1069847;
$contig_length[10]=$contig_length[9]+1057463;
$contig_length[11]=$contig_length[10]+1012305;
$contig_length[12]=$contig_length[11]+979369;
$contig_length[13]=$contig_length[12]+937107;
$contig_length[14]=$contig_length[13]+898261;
$contig_length[15]=$contig_length[14]+857155;
$contig_length[16]=$contig_length[15]+557602;
$contig_length[17]=$contig_length[16]+498254;
$contig_length[18]=$contig_length[17]+243426;

# This value is the window size for the sliding window analyses
$windowsize = 50000;
$stepsize = 10000;
$recentwindow = 50000;

# Read in fourth line to have the strains names
$listofnames = <DATAFILE>;
print OUTPUTFILEHETSITES $listofnames;
chomp ($listofnames);
@names = split(/\t/, $listofnames);

print @names;

# initialize the dataset to all missing data
for ($a=0; $a< $numstrains; $a++)
{
	for ($b=0; $b < $numsnps; $b++)
	{
		$dataset[$a][$b]="?";
	}
}

# Next we will read the data into an array of strains x snps ($dataset). There are two other arrays ($contig & $position) that index the location of the SNPs.

$snp=0;

while ($snp < $numsnps)
{
	$nextline = <DATAFILE>;
	chomp ($nextline);
	@genos = split(/\t/, $nextline);
# This following line removes the whitespace at the end of the line.
	$genos[$numstrains+1] =~ s/\s+$//;
	$contig[$snp] = $genos[0];
	$position[$snp] = $genos[1];
	for ($a=0;$a<$numstrains;$a++)
	{
		$dataset[$a][$snp] = $genos[$a+2];
	}
	$snp++;
}

# First goal is to characterize the heterozygosity of the SNPs by position, and by strain.

$snp = 0;
$HoStrain= ();
@violators = ((0) x $numsnps);

print "contig\tlocation\tH(obs)\tH(exp)\tFis\n";
print OUTPUTFILEHETGENOME "contig\tlocation\tH(obs)\tH(exp)\tFis\n";

while ($snp < $numsnps)
{
	$hos=0;
	$as=0;
	$cs=0;
	$gs=0;
	$ts=0;
	$ns=0;
	$tagged = 'no';
	print "\n$contig[$snp]\t$position[$snp]\t";
	for ($a=0; $a < $numstrains; $a++)
	{
		print "$dataset[$a][$snp]";
		if ($dataset[$a][$snp] =~ /K|M|R|S|W|Y/)
		{
			$hos++;
			$HoStrain[$a]++;
		}
		if ($dataset[$a][$snp] =~ /A|M|R|W/)
		{
			if ($dataset[$a][$snp] =~ /A/)
			{$as+=2;}
			else
			{$as++;}
		}
		if ($dataset[$a][$snp] =~ /C|M|S|Y/)
		{
			if ($dataset[$a][$snp] =~ /C/)
			{$cs+=2;}
			else
			{$cs++;}
		}
		if ($dataset[$a][$snp] =~ /G|R|S|K/)
		{
			if ($dataset[$a][$snp] =~ /G/)
			{$gs+=2;}
			else
			{$gs++;}
		}
		if ($dataset[$a][$snp] =~ /T|W|Y|K/)
		{
			if ($dataset[$a][$snp] =~ /T/)
			{$ts+=2;}
			else
			{$ts++;}
		}
		if ($dataset[$a][$snp] =~ /N/)
		{
			{$ns++;}
		}
	}
	if ((2*$numstrains-$ns*2) == 0 ) {print"\nhey snp=$snp";}
	$pA[$snp] = $as/(2*$numstrains-$ns*2);
	$pC[$snp] = $cs/(2*$numstrains-$ns*2);
	$pG[$snp] = $gs/(2*$numstrains-$ns*2);
	$pT[$snp] = $ts/(2*$numstrains-$ns*2);
	$numN[$snp] = $ns;
	if (($pA[$snp] != 0) && ($pA[$snp] == $pC[$snp])) {$maxallele[$snp] = 'M'; print " A and C";}
	elsif (($pA[$snp] != 0) && ($pA[$snp] == $pG[$snp])) {$maxallele[$snp] = 'R';print " A and G";}
	elsif (($pA[$snp] != 0) && ($pA[$snp] == $pT[$snp])) {$maxallele[$snp] = 'W';print " A and T";} 
	elsif (($pC[$snp] != 0) && ($pC[$snp] == $pG[$snp])) {$maxallele[$snp] = 'S';print " G and C";} 
	elsif (($pC[$snp] != 0) && ($pC[$snp] == $pT[$snp])) {$maxallele[$snp] = 'Y';print " T and C";}
	elsif (($pG[$snp] != 0) && ($pG[$snp] == $pT[$snp])) {$maxallele[$snp] = 'K';print " G and T";} 
	elsif (max ($pA[$snp], $pC[$snp], $pG[$snp], $pT[$snp]) eq $pA[$snp]) {$maxallele[$snp] = 'A';print " A dom";}
	elsif (max ($pA[$snp], $pC[$snp], $pG[$snp], $pT[$snp]) eq $pC[$snp]) {$maxallele[$snp] = 'C';print " C dom";}
	elsif (max ($pA[$snp], $pC[$snp], $pG[$snp], $pT[$snp]) eq $pG[$snp]) {$maxallele[$snp] = 'G';print " G dom";}
	elsif (max ($pA[$snp], $pC[$snp], $pG[$snp], $pT[$snp]) eq $pT[$snp]) {$maxallele[$snp] = 'T';print " T dom";}

# this bit tests if the site violates the infinite sites model and if so tags it and notes it in array violators
	if ((($pA[$snp] && $pC[$snp] && $pG[$snp]) != 0) | (($pA[$snp] && $pC[$snp] && $pT[$snp]) != 0) | (($pA[$snp] && $pG[$snp] && $pT[$snp]) != 0) | (($pC[$snp] && $pG[$snp] && $pT[$snp]) != 0))
	{
		print OUTPUTFILEVIOL "\nlocus $contig[$snp]\t$position[$snp] violates the infinite sites model";
		print "locus $contig[$snp]\t$position[$snp] violates the infinite sites model"; 
		$tagged='yes';
# here we make an array that has a list of all 'tagged' snps that violate infinite sites
		$violators[$snp] = 1;
	
	}
	$He = 1 - $pA[$snp]**2 - $pC[$snp]**2 - $pG[$snp]**2 - $pT[$snp]**2;	
	unless (($He == 0) | ($tagged eq 'yes') | ($numN[$snp] == ($numstrains-1)))
	{
		$Ho = $hos/($numstrains-$ns);
#traditional biased estimator not used
		$Fisbias = 1 - ($Ho/$He);
# equation from pg. 80 in GDA
		$Fis = 1 - ((($numstrains-$ns-1)*$Ho)/((($numstrains-$ns)*$He)-($hos/(2*($numstrains-$ns)))));
# This following line will print to a file that can be readily analyzed in a statistical program
		unless (($contig[$snp] > $threshold) | ($He ==0))
		{
			printf "\t%.3f\t%.3f\t%.3f\t%.3f", $Ho, $He, $Fis, $Fisbias;
			printf OUTPUTFILEHETGENOME "$contig[$snp]\t$position[$snp]\t%.3f\t%.3f\t%.3f\n", $Ho, $He, $Fis;
			for ($d=0; $d < $numstrains; $d++)
			{
				print OUTPUTFILEHETSITES "$dataset[$d][$snp]\t";
			}
			print OUTPUTFILEHETSITES "\n";
		}
	}
	$hos=0;
	$as=0;
	$cs=0;
	$gs=0;
	$ts=0;
	$ns=0;
	$snp++;
	$tagged = 'no';
}
	
# this code prints out the Ho of each strain

print OUTPUTFILEHETSTRAIN "Observed heterozygosity\n";
for ($a=0;$a<$numstrains;$a++)
{
	print OUTPUTFILEHETSTRAIN "$names[$a]\t$HoStrain[$a]\n";
}

# this code computes the LD between each pair of SNPs
# for a roughly 40,000 matrix this is 1.6 billion calculations

# This subroutine makes a global array that houses the last SNP of each contig
getmaxlengthcontig ();
#rununiques ();
slidingwindow ();
#findrecentmuts ();
calculatedistance ();
#calculateLD ();

sub calculateLD {
$snp = 0;
$secondsnp = 1;
print OUTPUTFILELD "\nSupercontig\tcomparison\tdistance\tdeltaAB\n";

while ($snp < $numsnps-1)
{
print "\nLD for marker # $snp";
# This line will skip to the next snp if the first snp has missing data for all but one strain
	while ($secondsnp < $numsnps)
	{
# This line will skip to the next snp if the first snp has missing data for all but one strain
#		if ($numN[$secondsnp] == ($numstrains-1)) {next;}
# this line stops the process if either locus violates the infinite sites model
		unless (($violators[$snp] == 1) | ($violators[$secondsnp] == 1) | ($numN[$snp]==($numstrains-1)) | ($numN[$secondsnp]==($numstrains-1)))
		{
			if ($contig[$snp] == $contig[$secondsnp])
			{$comparison='w';} else {$comparison='b';}
# this subroutine call gets all of the values needed for the calculation of the deltaab
			($nind,$fA,$fB,$n1,$n2,$n4,$n5) = calculatens ($snp, $secondsnp);
#			if ($nind==0) {for ($z=0;$z<$numstrains;$z++){print $dataset[$z][$snp];} print " "; for ($z=0;$z<$numstrains;$z++){print $dataset[$z][$secondsnp];} $topo=<STDIN>;}
# This makes a lot of sense. Don't calculate the values when there are less than 4 points of data. This will save a lot of disk space and hopefully smooth out data driven by few points.
			unless ($nind<4)
			{$nab = (2*$n1 + $n2 + $n4 + (0.5*$n5));
			$deltaab = abs( (1/$nind)*$nab - (2*$fA*$fB));
			print OUTPUTFILELD "$contig[$snp]\t$position[$snp]\t$contig[$secondsnp]\t$position[$secondsnp]\t$comparison\t";
#			print "1stcontig=$contig[$snp]\t$position[$snp]\t2ndcontig=$contig[$secondsnp]\t$position[$secondsnp]$comparison\t";
			if ($comparison =~ /w/)
			{
				$distance=$position[$secondsnp]-$position[$snp]; 
				print OUTPUTFILELD "$distance"; 
#				print "$distance";
			} 
			else 
			{
				print OUTPUTFILELD "\t"; 
#				print "\t";
			}
			printf OUTPUTFILELD "\t%.3f\n", $deltaab;
#			print "\t$deltaab";
			}
		}
		$secondsnp++;
	}
	$snp++;
	$secondsnp=0;
}
}

sub calculatens {

my ($base1, $base2) = @_;

my ($sumn1s) = 0;
my ($sumn2s) = 0;
my ($sumn4s) = 0;
my ($sumn5s) = 0;
my ($summiss) = 0;

if (($maxallele[$base1] =~ /A|M|R|W/) && ($maxallele[$base2] =~ /A|M|R|W/)) 	
{
	for ($a=0; $a < $numstrains; $a++)
	{
			if (($dataset[$a][$base1] =~ /A/) && ($dataset[$a][$base2] =~ /A/))
			{$sumn1s++;}
			elsif (($dataset[$a][$base1] =~ /M|R|W/) && ($dataset[$a][$base2] =~ /M|R|W/))
			{$sumn5s++;}
			elsif (($dataset[$a][$base1] =~ /A/) && ($dataset[$a][$base2] =~ /M|R|W/))
			{$sumn2s++;}
			elsif (($dataset[$a][$base1] =~ /M|R|W/) && ($dataset[$a][$base2] =~ /A/))
			{$sumn4s++;}
			elsif (($dataset[$a][$base1] =~ /N/) | ($dataset[$a][$base2] =~ /N/))
			{$summiss++;}

	}		 
}
elsif (($maxallele[$base1] =~ /A|M|R|W/) && ($maxallele[$base2] =~ /C|S|Y/)) 	
{
	for ($a=0; $a < $numstrains; $a++)
	{
			if (($dataset[$a][$base1] =~ /A/) && ($dataset[$a][$base2] =~ /C/))
			{$sumn1s++;}
			elsif (($dataset[$a][$base1] =~ /M|R|W/) && ($dataset[$a][$base2] =~ /M|S|Y/))
			{$sumn5s++;}
			elsif (($dataset[$a][$base1] =~ /A/) && ($dataset[$a][$base2] =~ /M|S|Y/))
			{$sumn2s++;}
			elsif (($dataset[$a][$base1] =~ /M|R|W/) && ($dataset[$a][$base2] =~ /C/))
			{$sumn4s++;}
			elsif (($dataset[$a][$base1] =~ /N/) | ($dataset[$a][$base2] =~ /N/))
			{$summiss++;}
	}		 
}
elsif (($maxallele[$base1] =~ /A|M|R|W/) && ($maxallele[$base2] =~ /G|K/)) 	
{
	for ($a=0; $a < $numstrains; $a++)
	{
			if (($dataset[$a][$base1] =~ /A/) && ($dataset[$a][$base2] =~ /G/))
			{$sumn1s++;}
			elsif (($dataset[$a][$base1] =~ /M|R|W/) && ($dataset[$a][$base2] =~ /R|S|K/))
			{$sumn5s++;}
			elsif (($dataset[$a][$base1] =~ /A/) && ($dataset[$a][$base2] =~ /R|S|K/))
			{$sumn2s++;}
			elsif (($dataset[$a][$base1] =~ /M|R|W/) && ($dataset[$a][$base2] =~ /G/))
			{$sumn4s++;}
			elsif (($dataset[$a][$base1] =~ /N/) | ($dataset[$a][$base2] =~ /N/))
			{$summiss++;}
	}		 
}
elsif (($maxallele[$base1] =~ /A|M|R|W/) && ($maxallele[$base2] =~ /T/)) 	
{
	for ($a=0; $a < $numstrains; $a++)
	{
			if (($dataset[$a][$base1] =~ /A/) && ($dataset[$a][$base2] =~ /T/))
			{$sumn1s++;}
			elsif (($dataset[$a][$base1] =~ /M|R|W/) && ($dataset[$a][$base2] =~ /W|Y|K/))
			{$sumn5s++;}
			elsif (($dataset[$a][$base1] =~ /A/) && ($dataset[$a][$base2] =~ /W|Y|K/))
			{$sumn2s++;}
			elsif (($dataset[$a][$base1] =~ /M|R|W/) && ($dataset[$a][$base2] =~ /G/))
			{$sumn4s++;}
			elsif (($dataset[$a][$base1] =~ /N/) | ($dataset[$a][$base2] =~ /N/))
			{$summiss++;}
	}		 
}
elsif (($maxallele[$base1] =~ /C|S|Y/) && ($maxallele[$base2] =~ /A|M|R|W/)) 	
{
	for ($a=0; $a < $numstrains; $a++)
	{
			if (($dataset[$a][$base1] =~ /C/) && ($dataset[$a][$base2] =~ /A/))
			{$sumn1s++;}
			elsif (($dataset[$a][$base1] =~ /M|S|Y/) && ($dataset[$a][$base2] =~ /M|R|W/))
			{$sumn5s++;}
			elsif (($dataset[$a][$base1] =~ /C/) && ($dataset[$a][$base2] =~ /M|R|W/))
			{$sumn2s++;}
			elsif (($dataset[$a][$base1] =~ /M|S|Y/) && ($dataset[$a][$base2] =~ /A/))
			{$sumn4s++;}
			elsif (($dataset[$a][$base1] =~ /N/) | ($dataset[$a][$base2] =~ /N/))
			{$summiss++;}
	}		 
}
elsif (($maxallele[$base1] =~ /C|S|Y/) && ($maxallele[$base2] =~ /C|S|Y/)) 	
{
	for ($a=0; $a < $numstrains; $a++)
	{
			if (($dataset[$a][$base1] =~ /C/) && ($dataset[$a][$base2] =~ /C/))
			{$sumn1s++;}
			elsif (($dataset[$a][$base1] =~ /M|S|Y/) && ($dataset[$a][$base2] =~ /M|S|Y/))
			{$sumn5s++;}
			elsif (($dataset[$a][$base1] =~ /C/) && ($dataset[$a][$base2] =~ /M|S|Y/))
			{$sumn2s++;}
			elsif (($dataset[$a][$base1] =~ /M|S|Y/) && ($dataset[$a][$base2] =~ /C/))
			{$sumn4s++;}
			elsif (($dataset[$a][$base1] =~ /N/) | ($dataset[$a][$base2] =~ /N/))
			{$summiss++;}
	}		 
}
elsif (($maxallele[$base1] =~ /C|S|Y/) && ($maxallele[$base2] =~ /G|K/)) 	
{
	for ($a=0; $a < $numstrains; $a++)
	{
			if (($dataset[$a][$base1] =~ /C/) && ($dataset[$a][$base2] =~ /G/))
			{$sumn1s++;}
			elsif (($dataset[$a][$base1] =~ /M|S|Y/) && ($dataset[$a][$base2] =~ /R|S|K/))
			{$sumn5s++;}
			elsif (($dataset[$a][$base1] =~ /C/) && ($dataset[$a][$base2] =~ /R|S|K/))
			{$sumn2s++;}
			elsif (($dataset[$a][$base1] =~ /M|S|Y/) && ($dataset[$a][$base2] =~ /G/))
			{$sumn4s++;}
			elsif (($dataset[$a][$base1] =~ /N/) | ($dataset[$a][$base2] =~ /N/))
			{$summiss++;}
	}		 
}
elsif (($maxallele[$base1] =~ /C|S|Y/) && ($maxallele[$base2] =~ /T/)) 	
{
	for ($a=0; $a < $numstrains; $a++)
	{
			if (($dataset[$a][$base1] =~ /C/) && ($dataset[$a][$base2] =~ /T/))
			{$sumn1s++;}
			elsif (($dataset[$a][$base1] =~ /M|S|Y/) && ($dataset[$a][$base2] =~ /W|Y|K/))
			{$sumn5s++;}
			elsif (($dataset[$a][$base1] =~ /C/) && ($dataset[$a][$base2] =~ /W|Y|K/))
			{$sumn2s++;}
			elsif (($dataset[$a][$base1] =~ /M|S|Y/) && ($dataset[$a][$base2] =~ /T/))
			{$sumn4s++;}
			elsif (($dataset[$a][$base1] =~ /N/) | ($dataset[$a][$base2] =~ /N/))
			{$summiss++;}
	}		 
}
elsif (($maxallele[$base1] =~ /G|K/) && ($maxallele[$base2] =~ /A|M|R|W/)) 	
{
	for ($a=0; $a < $numstrains; $a++)
	{
			if (($dataset[$a][$base1] =~ /G/) && ($dataset[$a][$base2] =~ /A/))
			{$sumn1s++;}
			elsif (($dataset[$a][$base1] =~ /R|S|K/) && ($dataset[$a][$base2] =~ /M|R|W/))
			{$sumn5s++;}
			elsif (($dataset[$a][$base1] =~ /G/) && ($dataset[$a][$base2] =~ /M|R|W/))
			{$sumn2s++;}
			elsif (($dataset[$a][$base1] =~ /R|S|K/) && ($dataset[$a][$base2] =~ /A/))
			{$sumn4s++;}
			elsif (($dataset[$a][$base1] =~ /N/) | ($dataset[$a][$base2] =~ /N/))
			{$summiss++;}
	}		 
}
elsif (($maxallele[$base1] =~ /G|K/) && ($maxallele[$base2] =~ /C|S|Y/)) 	
{
	for ($a=0; $a < $numstrains; $a++)
	{
			if (($dataset[$a][$base1] =~ /G/) && ($dataset[$a][$base2] =~ /C/))
			{$sumn1s++;}
			elsif (($dataset[$a][$base1] =~ /R|S|K/) && ($dataset[$a][$base2] =~ /M|S|Y/))
			{$sumn5s++;}
			elsif (($dataset[$a][$base1] =~ /G/) && ($dataset[$a][$base2] =~ /M|S|Y/))
			{$sumn2s++;}
			elsif (($dataset[$a][$base1] =~ /R|S|K/) && ($dataset[$a][$base2] =~ /C/))
			{$sumn4s++;}
			elsif (($dataset[$a][$base1] =~ /N/) | ($dataset[$a][$base2] =~ /N/))
			{$summiss++;}
	}		 
}
elsif (($maxallele[$base1] =~ /G|K/) && ($maxallele[$base2] =~ /G|K/)) 	
{
	for ($a=0; $a < $numstrains; $a++)
	{
			if (($dataset[$a][$base1] =~ /G/) && ($dataset[$a][$base2] =~ /G/))
			{$sumn1s++;}
			elsif (($dataset[$a][$base1] =~ /R|S|K/) && ($dataset[$a][$base2] =~ /R|S|K/))
			{$sumn5s++;}
			elsif (($dataset[$a][$base1] =~ /G/) && ($dataset[$a][$base2] =~ /R|S|K/))
			{$sumn2s++;}
			elsif (($dataset[$a][$base1] =~ /R|S|K/) && ($dataset[$a][$base2] =~ /G/))
			{$sumn4s++;}
			elsif (($dataset[$a][$base1] =~ /N/) | ($dataset[$a][$base2] =~ /N/))
			{$summiss++;}
	}		 
}
elsif (($maxallele[$base1] =~ /G|K/) && ($maxallele[$base2] =~ /T/)) 	
{
	for ($a=0; $a < $numstrains; $a++)
	{
			if (($dataset[$a][$base1] =~ /G/) && ($dataset[$a][$base2] =~ /T/))
			{$sumn1s++;}
			elsif (($dataset[$a][$base1] =~ /R|S|K/) && ($dataset[$a][$base2] =~ /W|Y|K/))
			{$sumn5s++;}
			elsif (($dataset[$a][$base1] =~ /G/) && ($dataset[$a][$base2] =~ /W|Y|K/))
			{$sumn2s++;}
			elsif (($dataset[$a][$base1] =~ /R|S|K/) && ($dataset[$a][$base2] =~ /T/))
			{$sumn4s++;}
			elsif (($dataset[$a][$base1] =~ /N/) | ($dataset[$a][$base2] =~ /N/))
			{$summiss++;}
	}		 
}
elsif (($maxallele[$base1] =~ /T/) && ($maxallele[$base2] =~ /A|M|R|W/)) 	
{
	for ($a=0; $a < $numstrains; $a++)
	{
			if (($dataset[$a][$base1] =~ /T/) && ($dataset[$a][$base2] =~ /A/))
			{$sumn1s++;}
			elsif (($dataset[$a][$base1] =~ /W|Y|K/) && ($dataset[$a][$base2] =~ /M|R|W/))
			{$sumn5s++;}
			elsif (($dataset[$a][$base1] =~ /T/) && ($dataset[$a][$base2] =~ /M|R|W/))
			{$sumn2s++;}
			elsif (($dataset[$a][$base1] =~ /W|Y|K/) && ($dataset[$a][$base2] =~ /A/))
			{$sumn4s++;}
			elsif (($dataset[$a][$base1] =~ /N/) | ($dataset[$a][$base2] =~ /N/))
			{$summiss++;}
	}		 
}
elsif (($maxallele[$base1] =~ /T/) && ($maxallele[$base2] =~ /C|S|Y/)) 	
{
	for ($a=0; $a < $numstrains; $a++)
	{
			if (($dataset[$a][$base1] =~ /T/) && ($dataset[$a][$base2] =~ /C/))
			{$sumn1s++;}
			elsif (($dataset[$a][$base1] =~ /W|Y|K/) && ($dataset[$a][$base2] =~ /M|S|Y/))
			{$sumn5s++;}
			elsif (($dataset[$a][$base1] =~ /T/) && ($dataset[$a][$base2] =~ /M|S|Y/))
			{$sumn2s++;}
			elsif (($dataset[$a][$base1] =~ /W|Y|K/) && ($dataset[$a][$base2] =~ /C/))
			{$sumn4s++;}
			elsif (($dataset[$a][$base1] =~ /N/) | ($dataset[$a][$base2] =~ /N/))
			{$summiss++;}
	}		 
}
elsif (($maxallele[$base1] =~ /T/) && ($maxallele[$base2] =~ /G|K/)) 	
{
	for ($a=0; $a < $numstrains; $a++)
	{
			if (($dataset[$a][$base1] =~ /T/) && ($dataset[$a][$base2] =~ /G/))
			{$sumn1s++;}
			elsif (($dataset[$a][$base1] =~ /W|Y|K/) && ($dataset[$a][$base2] =~ /R|S|K/))
			{$sumn5s++;}
			elsif (($dataset[$a][$base1] =~ /T/) && ($dataset[$a][$base2] =~ /R|S|K/))
			{$sumn2s++;}
			elsif (($dataset[$a][$base1] =~ /W|Y|K/) && ($dataset[$a][$base2] =~ /G/))
			{$sumn4s++;}
			elsif (($dataset[$a][$base1] =~ /N/) | ($dataset[$a][$base2] =~ /N/))
			{$summiss++;}
	}		 
}
elsif (($maxallele[$base1] =~ /T/) && ($maxallele[$base2] =~ /T/)) 	
{
	for ($a=0; $a < $numstrains; $a++)
	{
			if (($dataset[$a][$base1] =~ /T/) && ($dataset[$a][$base2] =~ /T/))
			{$sumn1s++;}
			elsif (($dataset[$a][$base1] =~ /W|Y|K/) && ($dataset[$a][$base2] =~ /W|Y|K/))
			{$sumn5s++;}
			elsif (($dataset[$a][$base1] =~ /T/) && ($dataset[$a][$base2] =~ /W|Y|K/))
			{$sumn2s++;}
			elsif (($dataset[$a][$base1] =~ /W|Y|K/) && ($dataset[$a][$base2] =~ /T/))
			{$sumn4s++;}
			elsif (($dataset[$a][$base1] =~ /N/) | ($dataset[$a][$base2] =~ /N/))
			{$summiss++;}
	}		 
}

if ($maxallele[$base1] =~ /K|M|R|S|W|Y/) {$freqpA=0.5;}
elsif ($maxallele[$base1] =~ /A/) {$freqpA=$pA[$base1];}
elsif ($maxallele[$base1] =~ /C/) {$freqpA=$pC[$base1];}
elsif ($maxallele[$base1] =~ /G/) {$freqpA=$pG[$base1];}
elsif ($maxallele[$base1] =~ /T/) {$freqpA=$pT[$base1];}
if ($maxallele[$base2] =~ /K|M|R|S|W|Y/) {$freqpB=0.5;}
elsif ($maxallele[$base2] =~ /A/) {$freqpB=$pA[$base2];}
elsif ($maxallele[$base2] =~ /C/) {$freqpB=$pC[$base2];}
elsif ($maxallele[$base2] =~ /G/) {$freqpB=$pG[$base2];}
elsif ($maxallele[$base2] =~ /T/) {$freqpB=$pT[$base2];}

$numberinds = $numstrains -$summiss;

return ($numberinds, $freqpA, $freqpB, $sumn1s, $sumn2s, $sumn4s, $sumn5s);
}

# This subroutine calculates the heterozygous density in a window of size $windowsize 

sub slidingwindow {
my $firstbp = 0;
my $lastbp = $windowsize;
my $workcontig = $contig[$firstbp];
my $middlebp = int($lastbp/2);
my @sumhets = ((0) x $numstrains);
my $totalbp = $middlebp+$contig_length[$workcontig];

print SLIDEHET "Contig\tBase\tCumBase\t";
for ($a=0;$a<$numstrains;$a++)
{
	print SLIDEHET "$names[$a]\t";
}
print SLIDEHET "\n";

while ($workcontig < $numcontigs)
{
	print "\ncontigsize[$workcontig] = $contigsize[$workcontig]";
	if ($contigsize[$workcontig] < $middlebp) {$workcontig++; next;}
	while ($middlebp<$contigsize[$workcontig]) 
	{
		for ($a=0;$a<$numstrains;$a++)
		{
#			for ($b=0;$b<$numsnps;$b++)
			for ($b=$firstcontigSNP{$workcontig};$b<$lastcontigSNP{$workcontig};$b++)
			{
				if (($position[$b]>= $firstbp) && ($position[$b] < $lastbp) && ($contig[$b] == $workcontig))
				{
#					print "\nposition[$b] = $position[$b] contig[$b]=$contig[$b] dataset[$a][$b] = $dataset[$a][$b]";
					if ($dataset[$a][$b] =~ /K|M|R|S|W|Y/)
					{$sumhets[$a]++;}
#					print "\nsumhets[$a]=$sumhets[$a]";
				}
			}
		}
		for ($b=0; $b<$numstrains; $b++)
		{
			$hetaverage[$b] = $sumhets[$b]/$windowsize;
		}
		print SLIDEHET "$workcontig\t$middlebp\t$totalbp\t";
		for ($a=0;$a<$numstrains;$a++)
		{
			printf SLIDEHET "%.6f\t",$hetaverage[$a];
		}
		print SLIDEHET "\n";
		$firstbp += $stepsize;
		$lastbp += $stepsize;
		$middlebp += $stepsize;
		$totalbp += $stepsize;
		@sumhets = ((0) x $numstrains);
		if ($middlebp > $contigsize[$workcontig])
		{
			$firstbp = 0;
			$lastbp = $windowsize;
			$middlebp = int($lastbp/2);
			$workcontig++; 
			$totalbp = $contig_length[$workcontig];
			print "\nworkcontiginslide=$workcontig";
		}
	}
}
}

# This subroutine determines the last bp of each contig which is used for sliding window analysis
# Also it calculates the range in the dataset that each contig covers. This speeds up the sliding window analysis.
sub getmaxlengthcontig {

@contigsize = ((0) x $numcontigs);
my $workcontig=0;
my $worksnp =0;

$firstcontigSNP{$contig[$worksnp]} = $worksnp;

while ($worksnp <= $numsnps)
{
	if ($contig[$worksnp] != $contig[$worksnp+1])
	{
		$contigsize[$contig[$worksnp]] = $position[$worksnp];
		$lastcontigSNP{$contig[$worksnp]} = $worksnp;
		$firstcontigSNP{$contig[$worksnp+1]} = $worksnp+1;		
		print "\ncontigsize $contig[$worksnp] = $contigsize[$contig[$worksnp]]";
	}
	$worksnp++;
}
}

# This subroutine calculates each site where all strains are of one homozygous type except one strain that is uniquely heterozygous
sub rununiques {

print UNIQUEMUTS "Strain\tContig\tPosition\tDataset\tGenotype\n";
for ($a=0;$a<$numsnps;$a++)
{
	if ($numN[$a] == 0)
	{
		if (($pA[$a] == ((2*$numstrains-1)/(2*$numstrains)))|($pC[$a] == ((2*$numstrains-1)/(2*$numstrains)))|($pG[$a] == ((2*$numstrains-1)/(2*$numstrains))) | ($pT[$a] == ((2*$numstrains-1)/(2*$numstrains))))
		{
			for ($b=0;$b<$numstrains;$b++)
			{
				if ($dataset[$b][$a] =~ (/K|M|R|S|W|Y/))
				{
					print UNIQUEMUTS"\n$names[$b]\t$contig[$a]\t$position[$a]\t";
					for ($c=0;$c<$numstrains;$c++)
					{print UNIQUEMUTS "$dataset[$c][$a]";}
					print UNIQUEMUTS "\t$dataset[$b][$a]";
				}
			}
		}
	}
}
}
# This subroutine identifies recent mutations as those with single heterozygous positions that lie within a window of size $recentwindow
sub findrecentmuts {
my $firstbp = 0;
my $lastbp = $recentwindow;
my $workcontig = $contig[$firstbp];
my @sumhets = ((0) x $numstrains);

print RECENTMUTS "Contig\tBase\tGenotypes\tStrain\n";
while ($workcontig < $numcontigs)
{
	print "\nrecentmutcontigsize[$workcontig] = $contigsize[$workcontig]";
	if ($contigsize[$workcontig] < $lastbp) {$workcontig++; next;}
	while ($lastbp<$contigsize[$workcontig]) 
	{
		for ($a=0;$a<$numstrains;$a++)
		{
#			for ($b=0;$b<$numsnps;$b++)
			for ($b=$firstcontigSNP{$workcontig};$b<$lastcontigSNP{$workcontig};$b++)
			{
				if (($position[$b]>= $firstbp) && ($position[$b] < $lastbp) && ($contig[$b] == $workcontig))
				{
#					print "\nposition[$b] = $position[$b] contig[$b]=$contig[$b] dataset[$a][$b] = $dataset[$a][$b]";
					if ($dataset[$a][$b] =~ /K|M|R|S|W|Y/)
					{$sumhets[$a]++;}
#					print "\nsumhets[$a]=$sumhets[$a]";
				}
			}
		}
		for ($b=0; $b<$numstrains; $b++)
		{
			if ($sumhets[$b] == 1)
			{findandprintrecent ($workcontig,$firstbp,$lastbp,$b);}
		}
		$firstbp += $recentwindow;
		$lastbp += $recentwindow;
		@sumhets = ((0) x $numstrains);
		if ($lastbp > $contigsize[$workcontig])
		{
			$firstbp = 0;
			$lastbp = $recentwindow;
			$workcontig++;
		}
	}
}
}

sub findandprintrecent {

my ($cont, $firstnt, $lastnt, $ind) = @_;

# This loop uses a hash to look up where in the list of SNPs to start looking for a particular contig.
# This is because the range along a supercontig can not be specifically pinpointed to a range in the array $dataset

for ($a=$firstcontigSNP{$cont};$a<$lastcontigSNP{$cont};$a++)
{
		if (($position[$a]>= $firstnt) && ($position[$a] < $lastnt) && ($contig[$a] == $cont))
		{
				if ($dataset[$ind][$a] =~ /K|M|R|S|W|Y/)
				{
					$sumhetsites=0;
					for ($c=$firstcontigSNP{$cont};$c<$lastcontigSNP{$cont};$c++)
					{
						if (($position[$c]>= $newfirstnt) && ($position[$c] < $newlastnt) && ($contig[$c] == $cont))
						{		
							if ($dataset[$ind][$c] =~ /K|M|R|S|W|Y/)
							{
								$sumhetsites++;
							}
						}
					}
					unless (($sumhetsites > 1) | ($numN[$a] != 0))
					{			
						print RECENTMUTS "$cont\t$position[$a]\t";	
						print  "$cont\t$position[$a]\t";	
						for ($b=0; $b<$numstrains; $b++)
						{
							print RECENTMUTS "$dataset[$b][$a]";
							print "$dataset[$b][$a]";
						}
						print RECENTMUTS "\t$names[$ind]=$dataset[$ind][$a]\n";
						print "\t$names[$ind]=$dataset[$ind][$a]\n";
					}
				}
		}
}
}

sub calculatedistance {

for ($a=0;$a<$numstrains-1;$a++)
{
	for ($b=0;$b<$numstrains;$b++)
	{
		$sumdiffs[$a][$b] = 0; $summissing[$a][$b]=0; $sumdiffsnoNs[$a][$b] = 0; $summissingnoNs[$a][$b]=0;
	}
}
for ($a=0;$a<$numstrains-1;$a++)
{
	for ($b=$a+1;$b<$numstrains;$b++)
	{
		for ($c=0;$c<$numsnps;$c++)
		{
			if (($dataset[$a][$c] =~ /N/) | ($dataset[$b][$c] =~ /N/))
			{$summissing[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /A/) && ($dataset[$b][$c] =~ /C|G|T/))
			{$sumdiffs[$a][$b]+=2;}
			elsif (($dataset[$a][$c] =~ /A/) && ($dataset[$b][$c] =~ /M|R|W/))
			{$sumdiffs[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /C/) && ($dataset[$b][$c] =~ /A|G|T/))
			{$sumdiffs[$a][$b]+=2;}
			elsif (($dataset[$a][$c] =~ /C/) && ($dataset[$b][$c] =~ /M|S|Y/))
			{$sumdiffs[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /G/) && ($dataset[$b][$c] =~ /A|C|T/))
			{$sumdiffs[$a][$b]+=2;}
			elsif (($dataset[$a][$c] =~ /G/) && ($dataset[$b][$c] =~ /K|R|S/))
			{$sumdiffs[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /T/) && ($dataset[$b][$c] =~ /A|C|G/))
			{$sumdiffs[$a][$b]+=2;}
			elsif (($dataset[$a][$c] =~ /T/) && ($dataset[$b][$c] =~ /K|W|Y/))
			{$sumdiffs[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /K/) && ($dataset[$b][$c] =~ /G|T/))
			{$sumdiffs[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /M/) && ($dataset[$b][$c] =~ /A|C/))
			{$sumdiffs[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /R/) && ($dataset[$b][$c] =~ /A|G/))
			{$sumdiffs[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /S/) && ($dataset[$b][$c] =~ /C|G/))
			{$sumdiffs[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /W/) && ($dataset[$b][$c] =~ /A|T/))
			{$sumdiffs[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /Y/) && ($dataset[$b][$c] =~ /C|T/))
			{$sumdiffs[$a][$b]++;}
		}
	}
}

printf DISTANCES "%-10s\n",$names[$0];
printf  "\n%-10s\n",$names[$0];
for ($a=1;$a<$numstrains;$a++)
{
	printf DISTANCES "%-10s",$names[$a];
	printf  "%-10s",$names[$a];
	for ($b=0;$b<$a;$b++)
	{
		if (($numsnps-$summissing[$b][$a]) == 0) {$distance = 0; print "\nbadguys=$names[$b] vs. $names[$a]";}
		else {$distance = $sumdiffs[$b][$a] / (($numsnps-$summissing[$b][$a])*2);}
		printf DISTANCES "%6f\t", $distance;
		printf  "%6f\t", $distance;
	}
	print DISTANCES "\n";
	print "\n";
}

# This part prints a matrix that does not include any SNPs where any of the taxa are missing data
for ($a=0;$a<$numstrains-1;$a++)
{
	for ($b=$a+1;$b<$numstrains;$b++)
	{
		for ($c=0;$c<$numsnps;$c++)
		{
			if ($numN[$c] != 0)
			{$summissingnoNs[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /A/) && ($dataset[$b][$c] =~ /C|G|T/))
			{$sumdiffsnoNs[$a][$b]+=2;}
			elsif (($dataset[$a][$c] =~ /A/) && ($dataset[$b][$c] =~ /M|R|W/))
			{$sumdiffsnoNs[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /C/) && ($dataset[$b][$c] =~ /A|G|T/))
			{$sumdiffsnoNs[$a][$b]+=2;}
			elsif (($dataset[$a][$c] =~ /C/) && ($dataset[$b][$c] =~ /M|S|Y/))
			{$sumdiffsnoNs[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /G/) && ($dataset[$b][$c] =~ /A|C|T/))
			{$sumdiffsnoNs[$a][$b]+=2;}
			elsif (($dataset[$a][$c] =~ /G/) && ($dataset[$b][$c] =~ /K|R|S/))
			{$sumdiffsnoNs[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /T/) && ($dataset[$b][$c] =~ /A|C|G/))
			{$sumdiffsnoNs[$a][$b]+=2;}
			elsif (($dataset[$a][$c] =~ /T/) && ($dataset[$b][$c] =~ /K|W|Y/))
			{$sumdiffsnoNs[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /K/) && ($dataset[$b][$c] =~ /G|T/))
			{$sumdiffsnoNs[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /M/) && ($dataset[$b][$c] =~ /A|C/))
			{$sumdiffsnoNs[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /R/) && ($dataset[$b][$c] =~ /A|G/))
			{$sumdiffsnoNs[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /S/) && ($dataset[$b][$c] =~ /C|G/))
			{$sumdiffsnoNs[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /W/) && ($dataset[$b][$c] =~ /A|T/))
			{$sumdiffsnoNs[$a][$b]++;}
			elsif (($dataset[$a][$c] =~ /Y/) && ($dataset[$b][$c] =~ /C|T/))
			{$sumdiffsnoNs[$a][$b]++;}
		}
	}
}

printf DISTANCES "\n\nMatrix with no missing data\n\n %-10s\n",$names[$0];
printf  "\n\nMatrix with no missing data\n\n %-10s\n",$names[$0];
for ($a=1;$a<$numstrains;$a++)
{
	printf DISTANCES "%-10s",$names[$a];
	printf  "%-10s",$names[$a];
	for ($b=0;$b<$a;$b++)
	{
		$distance = $sumdiffsnoNs[$b][$a] / (($numsnps-$summissingnoNs[$b][$a])*2);
		printf DISTANCES "%6f\t", $distance;
		printf  "%6f\t", $distance;
	}
	print DISTANCES "\n";
	print "\n";
}

}
			

exit;
