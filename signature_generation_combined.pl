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

$| = 1;

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

my $sigFolder = "signaturesFeb\\";
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
my $infectionAllHandle = $connection->prepare("select virus_taxo_id from infection where host_taxo_id = ?");
my $infectionReviewedHandle = $connection->prepare("select distinct(i.virus_taxo_id) from infection i, source s, protein p where i.host_taxo_id = ? and p.virus_taxo_id = i.virus_taxo_id and s.source = p.source and s.reviewed > 0");

my $virusHandle = $connection->prepare("select taxo_id, name, virus from virus");
$virusHandle->execute();
my %viralNames;
my %viralTaxos;
while (my ($virus_taxo,$name, $virus)=$virusHandle->fetchrow_array()){
	$name=~s/ /_/g;
	$viralNames{$virus_taxo} = $name;
	$virus=~s/ /_/g;
	$viralTaxos{$virus_taxo} = $virus;
}

&dealOneHost(9606);
exit;

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
	my %viralProteins;
#	print "Processing host $host_taxo started at ", strftime("%a %b %e %H:%M:%S %Y", localtime($startTime)), "\n";
#	my $now = strftime "%a %b %e %H:%M:%S %Y", localtime;
#	print "haha: $now\n";

	#get peptides for host
	my %hostPeptides = %{&getPeptidesBySQL($reviewed, $host_taxo,"host")};
#	my @peps = keys %hostPeptides;
#	print "Peptides in host $host_taxo: ", (scalar @peps),"\n";

	#get all peptides for viruses
	my %virusPeptideSpecies; #keys are peptides, values are arrays of viruse speciess the peptide exists
	my $infectionHandle = $infectionReviewedHandle;
	$infectionHandle = $infectionAllHandle if ($reviewed == 0);
	$infectionHandle->execute($host_taxo);
	my $totalVirus = 0;
	while (my ($virus_taxo) = $infectionHandle->fetchrow_array()){
#%peptides contains all peptides for the given species, keys are peptides, values are counts
#		my %peptides = %{&getPeptidesForSpecies($virus_taxo,$proteinVirusHandle)};
		my %peptides = %{&getPeptidesBySQL($reviewed, $virus_taxo,"virus")};
		my @peptides = keys %peptides;
#		print "Peptides in virus $virus_taxo: ", (scalar @peptides),"\n";
		my $virus = $viralTaxos{$virus_taxo};
		foreach my $pep(@peptides){
			push (@{$virusPeptideSpecies{$pep}},$virus_taxo);
			$viralProteins{$virus}{$pep} += $peptides{$pep};
		}
		$totalVirus++;
	}

	my %virusPeptideViruses; #keys are peptides, values are arrays of viruses the peptide exists
	foreach my $virus (sort {$a cmp $b} keys %viralProteins){
		my %tmp = %{$viralProteins{$virus}};
		foreach my $pep (keys %tmp){
			push (@{$virusPeptideViruses{$pep}},$virus);
		}		
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

	foreach my $pep(sort keys %virusPeptideSpecies){
		my @viruses = @{$virusPeptideSpecies{$pep}};
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
				push (@{$virusSigSeqs{'species'}{$viruses[0]}},$pep);
				delete $virusPeptideViruses{$pep}; #already be a signature for species, must be signature for virus, but redundant
#			}else{
#				print "not signature $pep as multiple existances in @viruses\n";
			}
		}
	}

	foreach my $pep(sort keys %virusPeptideViruses){
		my @viruses = @{$virusPeptideViruses{$pep}};
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
				push (@{$virusSigSeqs{'virus'}{$viruses[0]}},$pep);
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

	open SEQ, ">${sigFolder}combined_in_host_${host_taxo}_${type}_sig.fasta";
	if(exists $virusSigSeqs{'species'}){
		my %virusSigSeqSpecies = %{$virusSigSeqs{'species'}};
		foreach my $virus_taxo(keys %virusSigSeqSpecies){
			my @sigPeptides = @{$virusSigSeqSpecies{$virus_taxo}};
			my $seq = join("",@sigPeptides);
			my $num = scalar @sigPeptides;
			print SEQ ">species_sig_$viralNames{$virus_taxo}_${virus_taxo}_in_${host_taxo}_total_$num\n$seq\n";

			$total += $num;
			$count++;
			$max = $num if ($num>$max);
			$min = $num if ($num<$min);
		}
	}
	if(exists $virusSigSeqs{'virus'}){
		my %virusSigSeqVirus = %{$virusSigSeqs{'virus'}};
		foreach my $virus(keys %virusSigSeqVirus){
			my @sigPeptides = @{$virusSigSeqVirus{$virus}};
			my $seq = join("",@sigPeptides);
			my $num = scalar @sigPeptides;
			print SEQ ">viral_sig_${virus}_in_${host_taxo}_total_$num\n$seq\n";

			$total += $num;
			$count++;
			$max = $num if ($num>$max);
			$min = $num if ($num<$min);
		}
	}
	close SEQ;
	print "For host $host_taxo there are total $total signature peptides in $count entries out of all $totalVirus related with max $max and min $min\n\n";
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