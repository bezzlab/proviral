#!/usr/bin/perl
package Proviral;

use strict;
use warnings;
use Exporter;
use DBD::mysql;

our @ISA= qw( Exporter );

# these CAN be exported.
our @EXPORT_OK = qw (checkHost checkSource readSource getPeptidesBySQL);
my $dsn = "dbi:mysql:proviral:localhost:3306";
my $connection = DBI->connect($dsn,"proviral","proviral") or die "Can't connect to the DB\n";#often named as dbh
my $sourceFile = "populate_source.tsv";
my $hostFile = "populate_host.tsv";

my $pepBySQLforHostAll = $connection->prepare("select pe.peptide_seq from protein p,peptide_existance pe where p.host_taxo_id = ? and p.accession = pe.protein_accession");
my $pepBySQLforHostReviewed = $connection->prepare("select pe.peptide_seq from source s, protein p, peptide_existance pe where s.reviewed > 0 and s.source = p.source and p.host_taxo_id = ? and p.accession = pe.protein_accession");
my $pepBySQLforVirusAll = $connection->prepare("select pe.peptide_seq from protein p,peptide_existance pe where p.virus_taxo_id = ? and p.accession = pe.protein_accession");
my $pepBySQLforVirusReviewed = $connection->prepare("select pe.peptide_seq from source s, protein p, peptide_existance pe where s.reviewed > 0 and s.source = p.source and p.virus_taxo_id = ? and p.accession = pe.protein_accession");

#reviewed parameter has value either 0 (all, i.e. reviewed+unreviewed) or 1 (reviewed)
#species parameter is for the taxonomy id
#type parameter has value either host or virus
sub getPeptidesBySQL(){
	my $reviewed = $_[0];
	my $species = $_[1];
	my $type = $_[2];
	my $handle;

	if($reviewed == 0){#all data
		if($type eq "host"){
			$handle = $pepBySQLforHostAll;
		}else{
			$handle = $pepBySQLforVirusAll;
		}
	}else{#reviewed only
		if($type eq "host"){
			$handle = $pepBySQLforHostReviewed;
		}else{
			$handle = $pepBySQLforVirusReviewed;
		}
	}
	my %hash;
	$handle->execute($species);
	while (my ($pep) = $handle->fetchrow_array()) {
		$hash{$pep}++;
	}
	return \%hash;
}

sub checkHost(){
	print "Checking hosts now.\n\n";
	my $ref = &readHost();
	my %hostParam = %{$ref};
	my %hostDB;
	my $hostHandle = $connection->prepare("select taxo_id, name, mnemonic from host");
	$hostHandle->execute();
	while (my ($taxo,$name,$mnemonic)=$hostHandle->fetchrow_array()){
		$hostDB{$taxo}{name}=$name;
		$hostDB{$taxo}{mnemonic}=$mnemonic;
	}
	foreach my $hostParam(keys %hostParam){
		if (exists $hostDB{$hostParam}){
			unless ($hostDB{$hostParam}{name} eq $hostParam{$hostParam}{name} 
				&& $hostDB{$hostParam}{mnemonic} eq $hostParam{$hostParam}{mnemonic}){
				print "The information for host with taxo id $hostParam in the database is different from the host file $hostFile\n";
				print "In the database: name <$hostDB{$hostParam}{name}> mnemonic <$hostDB{$hostParam}{mnemonic}>\n";
				print "From the file: name <$hostParam{$hostParam}{name}> mnemonic <$hostParam{$hostParam}{mnemonic}>\n\n";
			}
		}else{
			print "New host with taxonomy $hostParam in the file, which will be inserted into table\n";
			my $insertHostHandle = $connection -> prepare("insert into host values ($hostParam,'$hostParam{$hostParam}{name}','$hostParam{$hostParam}{mnemonic}')");
			$insertHostHandle -> execute();
		}
	}
}

sub checkSource(){
	print "Checking sources now.\n\n";
	my %sourceParam = %{$_[0]};
	my %sourceDB;
	my $sourceHandle = $connection->prepare("select source,version,reviewed from source");
	$sourceHandle->execute();
	print "Existing sources in the database:\n";
	while (my ($source,$version,$reviewed)=$sourceHandle->fetchrow_array()){
		$sourceDB{$source}{version} = $version;
		$sourceDB{$source}{reviewed} = $reviewed;
		print "<$source>\t<$version>\t<$reviewed>\n";
	}
	foreach my $sourceParam (keys %sourceParam){
		if (exists $sourceDB{$sourceParam}){
			unless ($sourceParam{$sourceParam}{version} eq $sourceDB{$sourceParam}{version}
				&& $sourceParam{$sourceParam}{reviewed} eq $sourceDB{$sourceParam}{reviewed}) {
				print "The information for source with name $sourceParam in the database is different from the given source\n";
				print "In the database: version <$sourceDB{$sourceParam}{version}> reviewed <$sourceDB{$sourceParam}{reviewed}>\n";
				print "From the given value: version <$sourceParam{$sourceParam}{version}> reviewed <$sourceParam{$sourceParam}{reviewed}>\n\n";
			}
		}else{
			print "There is a new source $sourceParam given which has been inserted into table automatically\n";
			my $insertSourceHandle = $connection -> prepare("insert into source values ('$sourceParam','$sourceParam{$sourceParam}{version}',$sourceParam{$sourceParam}{reviewed})");
			$insertSourceHandle->execute();
		}
	}

}

sub readSource(){
	my %result;
	open IN, "$sourceFile";
	while (my $line = <IN>){
		next if (substr($line,0,1) eq "#");
		chomp $line;
		my ($source, $version, $reviewed) = split("\t",$line);
		$result{$source}{version} = $version;
		$result{$source}{reviewed} = $reviewed;
	}
	return \%result;
}

sub readHost(){
	my %result;
	open IN, "$hostFile";
	while (my $line = <IN>){
		next if (substr($line,0,1) eq "#");
		chomp $line;
		my ($taxo, $name, $mnemonic) = split("\t",$line);
		$result{$taxo}{name} = $name;
		$result{$taxo}{mnemonic} = $mnemonic;
	}
	return \%result;
}

1;