#this is for the simplest situation: static one, not dynamic, i.e. based on user selection
#for single host to calculate for all related viruses, again static (all), not dynamic (some virus definitely in/not in the sample by user)
#again while calculating strain signature, only consider that virus, not any other virus(es) found by virus signature

#select a host by taxo id

#from infection table, get all viruses
#for each virus, get all proteins group by virus column in virus table
#for each virus group, get all tryptic peptides
#find unique peptide within virus groups
use strict;
use DBD::mysql;
use POSIX qw(strftime);
use Proviral qw(getPeptidesBySQL);

unless (scalar @ARGV==1) {
	print "Error: Wrong number of parameters\n\n";
	&usage();
}
my $reviewed = $ARGV[0];
unless ($reviewed=~/^\d+$/){
	print "The parameter should be integer, ideally 0 (all) or 1 (reviewed only)\n";
	&usage();
}

my $type = "reviewed";
$type = "all" if ($reviewed == 0);

my $sigFolder = "signaturesDesktop\\";
$sigFolder =~s/ /_/g;
unless (-d $sigFolder){
	system ("mkdir $sigFolder");
	print "The folder <$sigFolder> has been created\n\n";
}
print "The signature files will be stored in the folder $sigFolder\n\n";

my $totalStart = time();
my $dsn = "dbi:mysql:proviral:localhost:3306";
my $connection = DBI->connect($dsn,"proviral","proviral") or die "Can't connect to the DB\n";#often named as dbh
#my $proteinHostHandle = $connection->prepare("select accession from protein where host_taxo_id = ?");
#my $proteinVirusHandle = $connection->prepare("select accession from protein where virus_taxo_id = ?");
my $existanceHandle = $connection->prepare("select peptide_seq from peptide_existance where protein_accession = ?");
my $infectionHandle = $connection->prepare("select virus_taxo_id from infection where host_taxo_id = ?");

my $virusHandle = $connection->prepare("select taxo_id, name from virus");
$virusHandle->execute();
my %viralNames;
while (my ($virus_taxo,$name)=$virusHandle->fetchrow_array()){
	$name=~s/ /_/g;
	$viralNames{$virus_taxo} = $name;
}
#my $pepInOneGo = $connection->prepare("select pe.peptide_seq, i.virus_taxo_id from infection i, protein p, peptide_existance pe where i.host_taxo_id = ? and i.virus_taxo_id = p.virus_taxo_id and p.accession = pe.protein_accession");

#my %abc = %{&getPeptidesBySQL(10247,$pepBySQLforVirus)};
#my %abc = %{&getPeptidesBySQL(9606,$pepBySQLforHost)};
#my %abc = %{&getPeptidesForSpecies(9606,$proteinHostHandle)};
#foreach my $pep(sort keys %abc){
#	print "$pep\t$abc{$pep}\n";
#}
#exit;

#&dealOneHost(9606);
#exit;

#open HOST_INTER,">populate_host_interference.tsv";
#open HOST_INTER_DETAIL,">populate_host_interference_detail.tsv";
#open SIG, ">populate_virus_level_signature.tsv";

my $hostHandle = $connection->prepare("select taxo_id from host");
$hostHandle->execute();
while(my ($host)=$hostHandle->fetchrow_array()){
	&dealOneHost($host);
}

print "Total run time for all hosts: ", time-$totalStart, " seconds\n";

sub dealOneHost(){
	my $host_taxo = $_[0];
	my $startTime = time();
	my %virusSigSeqs;
#	print "Processing host $host_taxo started at ", strftime("%a %b %e %H:%M:%S %Y", localtime($startTime)), "\n";
#	my $now = strftime "%a %b %e %H:%M:%S %Y", localtime;
#	print "haha: $now\n";

	#get peptides for host
	my %hostPeptides = %{&getPeptidesBySQL($reviewed, $host_taxo,"host")};
	my @peps = keys %hostPeptides;
#	print "Peptides in host $host_taxo: ", (scalar @peps),"\n";

	#get all peptides for viruses
	my %virusPeptides; #keys are peptides, values are arrays of viruses the peptide exists
	$infectionHandle->execute($host_taxo);
	my $totalVirus = 0;
	while (my ($virus_taxo) = $infectionHandle->fetchrow_array()){
#%peptides contains all peptides for the given species, keys are peptides, values are counts
#		my %peptides = %{&getPeptidesForSpecies($virus_taxo,$proteinVirusHandle)};
		my %peptides = %{&getPeptidesBySQL($reviewed, $virus_taxo,"virus")};
		my @peptides = keys %peptides;
#		print "Peptides in virus $virus_taxo: ", (scalar @peptides),"\n";
		foreach my $pep(@peptides){
			push (@{$virusPeptides{$pep}},$virus_taxo);
		}
		$totalVirus++;
	}
	#get all peptides from all viruses infecting the same host by one SQL
	#the result showed that not quicker than byVirus method, kept here only for reference
#	my %tmp; #one peptide may appear in multiple proteins in the same species, here we only care about species, not particular proteins
#	$pepInOneGo->execute($host_taxo);
#	while (my ($pep,$virus_taxo)=$pepInOneGo->fetchrow_array()){
#		$tmp{$pep}{$virus_taxo}++;
#	}
#	foreach my $pep (keys %tmp){
#		@{$virusPeptides{$pep}}= sort {$a <=> $b} keys %{$tmp{$pep}};
#	}

	foreach my $pep(sort keys %virusPeptides){
		my @viruses = @{$virusPeptides{$pep}};
		if (exists $hostPeptides{$pep}){ #the peptide exist in the host, not signature
#			print HOST_INTER "$pep\t$host_taxo\t$hostPeptides{$pep}\n";
#			if ($hostPeptides{$pep} == 1){
#				print "proteotypic peptide $pep for host $host_taxo can also be found in virus @viruses\n";
#			}else{
#				print "peptide $pep in both host $host_taxo with $hostPeptides{$pep} proteins and virus @viruses\n";
#			}
#			foreach my $virus_taxo(@viruses){
#				print HOST_INTER_DETAIL "$pep\t$host_taxo\t$virus_taxo\n";
#			}
		}else{#not exist in host
			my $count = scalar @viruses;
			if ($count == 1){
#				print SIG "$pep\t$viruses[0]\n";
				push (@{$virusSigSeqs{$viruses[0]}},$pep);
#				print "signature $pep in $viruses[0]\n";
#			}else{
#				print "not signature $pep as multiple existances in @viruses\n";
			}
		}
	}

	my $endTime = time();
	print "Run time for processing host $host_taxo : ", $endTime-$startTime, " seconds\n";

	return if (scalar keys %virusSigSeqs == 0);

	my $total = 0;
	my $count = 0;
	my $min = 99999999;
	my $max = 0;

	open SEQ, ">${sigFolder}species_in_host_${host_taxo}_${type}_sig.fasta";
	foreach my $virus_taxo(keys %virusSigSeqs){
		my @sigPeptides = @{$virusSigSeqs{$virus_taxo}};
		my $seq = join("",@sigPeptides);
		my $num = scalar @sigPeptides;
		print SEQ ">species_sig_$viralNames{$virus_taxo}_${virus_taxo}_in_${host_taxo}_total_$num\n$seq\n";

		$total += $num;
		$count++;
		$max = $num if ($num>$max);
		$min = $num if ($num<$min);
	}
	print "For host $host_taxo there are total $total signature peptides in $count viruses out of all $totalVirus related with max $max and min $min\n\n";
}

sub getPeptidesForSpecies(){
	my $species = $_[0];
	my $handle = $_[1];
	my %hash;
	$handle->execute($species);
	while (my ($acc) = $handle->fetchrow_array()) {
		$existanceHandle->execute($acc);
		while(my ($pep)=$existanceHandle->fetchrow_array()){
			$hash{$pep}++;
		}
	}
	return \%hash;
}

sub usage(){
	print "Usage: signature_generation_species.pl <reviewed>\n";
	print "This script generates the signatures at the strain level, i.e. for each taxonomy\n";
	print "The reviewed parameter indicate whether only include (value as 1) sequences from reviewed sources (e.g. swissprot) or not (value as 0)\n";
	exit;
}