#!/usr/bin/perl
use strict;
use Proviral qw(checkSource readSource);

my $numArg = scalar @ARGV;
if (($numArg!=1)){
	print "Error: The number of parameters is only allowed to be 1.\n";
	&usage();
}
my $type = lc($ARGV[0]);
unless ($type eq "host" || $type eq "virus"){
	print "Error: The allowed values for the type are 'host' or 'virus'\n";
	&usage();
}

my $ref = &readSource();
my %sources = %{$ref};
&checkSource($ref);

my %totalProtein;
my %totalPeptide;
my %totalVirus;
my %totalInfection;
#my %totalPeptideExistance;
#my %totalProduct;

foreach my $source (keys %sources){
	my $version = $sources{$source}{version};
	my $dataFolder = "e:\\\\proviral\\\\datafiles\\\\$type\\\\$source\\\\$version";
	$dataFolder =~s/ /_/g;
	unless (-d $dataFolder){
		print "No data files found for data source $source with version $version\n";
		print "Please make sure that all files for that source should be kept under folder $dataFolder\n";
		next;
	}
	open IN, "$dataFolder\\\\populate_protein_${type}_${source}_${version}.tsv";
	my %proteins;
	while (my $line = <IN>){
		chomp $line;
		my ($ac) = split("\t",$line);
		$totalProtein{$ac}++;
		$proteins{$ac}++;
	}

	open IN, "$dataFolder\\\\populate_peptide_${type}_${source}_${version}.tsv";
	my %peptides;
	while (my $line = <IN>){
		chomp $line;
		my ($seq) = split("\t",$line);
		$totalPeptide{$seq}++;
		$peptides{$seq}++;
	}

	open IN, "$dataFolder\\\\populate_peptide_existance_${type}_${source}_${version}.tsv";
	my %peptideExistance;
	while (my $line = <IN>){
		chomp $line;
		my ($pep, $pro) = split("\t",$line);
		my $key = "$pep-$pro";
#		$totalPeptideExistance{$key}++;
		$peptideExistance{$key}++;
	}

#	open IN, "$dataFolder\\\\populate_product_${type}_${source}_${version}.tsv";
#	my %products;
#	while (my $line = <IN>){
#		chomp $line;
#		my ($pep, $ionType, $serial) = split("\t",$line);
#		my $key = "$pep-$ionType-$serial";
#		$totalProduct{$key}++;
#		$products{$key}++;
#	}

	print "For $type in the source $source $version\n";
	my $count = scalar keys %proteins;
	print "The protein number is $count\n";
	$count = scalar keys %peptides;
	print "The peptide number is $count\n";
	$count = scalar keys %peptideExistance;
	print "The peptide existance number is $count\n";
	#$count = scalar keys %products;
	#print "The product number is $count\n";

	if($type eq "virus"){
		my %virus;
		my %infection;
		open IN, "$dataFolder\\\\populate_virus_${type}_${source}_${version}.tsv";
		while (my $line = <IN>){
			chomp $line;
			my ($virus_taxo) = split("\t",$line);
			$totalVirus{$virus_taxo}++;
			$virus{$virus_taxo}++;
		}
		open IN, "$dataFolder\\\\populate_infection_${type}_${source}_${version}.tsv";
		while (my $line = <IN>){
			chomp $line;
			my ($host_taxo, $virus_taxo) = split("\t",$line);
			$totalInfection{"$host_taxo-$virus_taxo"}++;
			$infection{"$host_taxo-$virus_taxo"}++;
		} 
		$count = scalar keys %virus;
		print "The virus number is $count\n";
		$count = scalar keys %infection;
		print "The infection number is $count\n";
	}
	print "\n";
}

my $count = scalar keys %totalProtein;
print "The total protein number is $count\n";
$count = scalar keys %totalPeptide;
print "The total peptide number is $count\n";
#$count = scalar keys %totalPeptideExistance;
#print "The total peptide existance number is $count\n";
#$count = scalar keys %totalProduct;
#print "The total product number is $count\n";
if ($type eq "virus"){
	$count = scalar keys %totalVirus;
	print "The total virus number is $count\n";
	$count = scalar keys %totalInfection;
	print "The total infection number is $count\n";
}


sub usage(){
	print "Usage: perl getNumber.pl <host|virus>\n";
	print "This script returns the number of entries of proteins, peptides and peptide existances for all sources for the given type and virus and infection number for virus only\n";
	print "Note 1: Due to the limitation of memory, the product ion part has coded and commented out\n";
	print 'Note 2: the path ($dataFolder) indicating where the data files will be stored is hard coded at line 30, please change accordingly';
	exit 1;
}


