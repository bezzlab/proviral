#!/usr/bin/perl
use strict;

use Proteomic::Misc qw(getHydrophobicity ionCalculation calculateMZ);
use Proviral qw(checkSource checkHost);

#hardcoded for sources
#every time need to change the $sourceType and $version variable, which will become one parameter in the command line
my $sourceType = "Swissprot";
$sourceType = "Trembl";
my $version = "Mar 2015";
my $reviewed = 1;

my %newSource;
$newSource{$sourceType}{version} = $version;
$newSource{$sourceType}{reviewed} = $reviewed;
&checkSource(\%newSource);

#generate the data file folder for the given source and version
my $dataFolder = "e:\\proviral\\datafiles\\host\\$sourceType\\$version\\";
$dataFolder =~s/ /_/g;
unless (-d $dataFolder){
	system ("mkdir $dataFolder");
	print "The folder <$dataFolder> has been created\n\n";
}
print "The data files will be saved in the folder <$dataFolder>\n";

my $fastaFile = "test.fasta";
#$fastaFile = "E:\\uniprot data\\$version\\uniprot_sprot.fasta";
#$fastaFile = "C:\\Users\\Jun Fan\\Dropbox\\proviral\\hostSequences\\combinedSwissprotHost.fasta";
$fastaFile = "combinedHostMar2015Trembl.fasta";
#$fastaFile = "C:\\Users\\Jun\\Dropbox\\proviral\\hostSequences\\combinedTremblHost.fasta";

my %species;
my %filenames;
#the same file used to populate Host table
open IN, "populate_host.tsv" or die "Please provide the list of species as host in the populate_host.tsv file in the order of taxo_id name mnemonic";
while (my $line=<IN>){
	next if ($line=~/^#/);
	chomp $line;
	my ($taxo,$name,$mnemonic) = split("\t",$line);
	my $filename ="${dataFolder}host_seq_${name}_${sourceType}_$version.fasta";
	$filename=~s/ /_/g;
	if(-e $filename){
		print "file $filename already exists, now deleting\n";
		system ("del $filename");
	}
	$filenames{$taxo}=$filename;
	$species{$mnemonic} = $taxo;
}

open PRO, ">${dataFolder}populate_protein_host_${sourceType}_$version.tsv" or die "Could not open the output file ${dataFolder}populate_protein_host_${sourceType}_$version.tsv";
open EXISTANCE, ">${dataFolder}populate_peptide_existance_host_${sourceType}_$version.tsv";
open PEP, ">${dataFolder}populate_peptide_host_${sourceType}_$version.tsv";
open ION, ">${dataFolder}populate_product_host_${sourceType}_$version.tsv";

#Host with isoforms
#Inputs: 1) uniprot all species sequence 2)list of host #for swissprot it is feasible, trembl fasta is too big
#read in the host list
#read in the fasta file
#	check species in id
#	if in the host list do protein digest
my %peptidesDone;
my $header;
my $seq = "";
open IN, "$fastaFile" or die "Cannot find the specified fasta file $fastaFile";
while(my $line=<IN>){
	chomp($line);
	if ($line=~/^>/){
		if (length $seq > 0){
			&processProtein($header,$seq);
		}
		$header = $line;
		$seq = "";
	}else{
		$seq.=$line;
	}
}
&processProtein($header,$seq);

#process one protein
#Inputs: 1) protein id 2) protein sequence
#Steps:
#protein digest
#foreach peptide <= 25AA
#	save to peptide_existance table
#	calculate hydro, mz1, mz2, mz3, mz4
#	save to peptide table
sub processProtein(){
	my ($header, $seq) = @_;
	my $protein_id = "";
	if ($sourceType eq "Swissprot") {
		$protein_id = &processUniprotHeader ($header);
	}elsif ($sourceType eq "Trembl") {
		$protein_id = &processUniprotHeader ($header);
	}else{
		die "\nThe source type <$sourceType> has not been coded for.\n";
	}
	return if (length $protein_id == 0); 
	print PRO "$seq\n";

	my ($ac,$id,$taxo_id) = split("#",$protein_id);
	my $out = $filenames{$taxo_id};
	open OUT,">>$out";
	print OUT ">$ac|$id\n$seq\n";
	close OUT;
	
	my @peptides = split (/(?<=[KR])(?=[^P])/,$seq); #alternative RE /(?!P)(?<=[RK])/
	my %peptidesCount;
	foreach my $peptide (@peptides) {
		$peptidesCount{$peptide}++ if (length $peptide <= 30 && length $peptide >= 3);
	}
	foreach my $peptide (sort {$a cmp $b} keys %peptidesCount){
		next if ($peptide=~/[UOBZJX]/);
		print EXISTANCE "$peptide\t$ac\t$peptidesCount{$peptide}\n";
		next if (exists $peptidesDone{$peptide});
		$peptidesDone{$peptide}=1;
		my $len = length $peptide;
		my $hydro = &getHydrophobicity($peptide);
		my $mz = &calculateMZ($peptide);
		print PEP "$peptide\t$len\t$hydro\t$mz\tnull\n";
		my @results = @{&ionCalculation($peptide)};
		foreach my $line(@results){
			print ION "$line\n";
		}
	}
}

sub processUniprotHeader(){
	my $header = $_[0];
	if ($header=~/^>((sp|tr)\S+)\s/){
		my ($type,$ac,$id) = split(/\|/,$1);
#		print "$type,$ac,$id\n";
		my (undef, $species) = split("_",$id);
		my $taxo_id = 0;
		if (exists $species{$species}){ #check whether in the host list
			$taxo_id = $species{$species};
			print PRO "$ac\t$id\t$taxo_id\t\\N\t$sourceType\t";
		}else{
			return "";
		}
		return "$ac#$id#$taxo_id";
	}
}

#sort firstly by length, secondly by alphabetical order
#@uniquePep = sort {
#	if((length $a)<(length $b)) {return -1;}
#	elsif((length $a)>(length $b)) {return 1;}
#	else {return $a cmp $b}
#}  @uniquePep;
