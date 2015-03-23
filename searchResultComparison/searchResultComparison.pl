use strict;
my $numArg = scalar @ARGV;
if (($numArg%2==1)){
	print "Error: The number of parameters is odd.\n";
	&usage();
}

my %inputs;
my @inputs;
for (my $i=0;$i<$numArg;$i+=2){
#	push (@inputs,$ARGV[$i+1]);
#	$inputs{$ARGV[$i+1]} = $ARGV[$i];
#	print "<$ARGV[$i+1]> ==> <$ARGV[$i]>\n";
}

push (@inputs, "human");
#push (@inputs, "species_signature");
push (@inputs, "virus_signature");
$inputs{"human"}= "Peptide_result_list_reviewed_human.csv";
$inputs{"human"}= "Peptide_result_list_unreviewed_human.csv";
$inputs{"human"}= "Peptide_result_list_all_species_swissprot.csv";

#$inputs{"species_signature"}="Galaxy9-[peptide_list_for_signature_Human].csv";
$inputs{"virus_signature"}="Peptide_result_list_reviewed_human_species_signature.csv";
#$inputs{"virus_signature"}="Peptide_result_list_reviewed_human_virus_signature.csv";
#$inputs{"virus_signature"}="Peptide_result_list_reviewed_human_conserved_signature.csv";
#$inputs{"virus_signature"}="Peptide_result_list_unreviewed_human_species_signature.csv";
#$inputs{"virus_signature"}="Peptide_result_list_unreviewed_human_virus_signature.csv";
#$inputs{"virus_signature"}="Peptide_result_list_unreviewed_human_conserved_signature.csv";

my %data; #first key is the peptide sequence, second key is the dataset id, the value is the protein list in that dataset
foreach my $id(@inputs){
	my $file = $inputs{$id};
	&readMzidlibResultFile($id, $file);
}

my %result;
my %proteinResult;
foreach my $peptide (keys %data){
	my %dataset = %{$data{$peptide}};
	my @datasets = sort {$a cmp $b} keys %dataset;
	my $size = scalar @datasets;
	my $index  = join ("-",@datasets);
	push(@{$result{$size}{$index}},$peptide);
	foreach my $datasetID (@datasets){
		my %acc = %{$dataset{$datasetID}};
		foreach my $tmp(keys %acc){
			push(@{$proteinResult{$tmp}{$datasetID}},$peptide);
		}
	}
}

#$"="\n";
#print "Unique peptides\n";
#my %uniquePeptides = %{$result{"1"}};
#foreach my $id(keys %uniquePeptides){
#	print "For dataset $id\n";
#	my @peptides = @{$uniquePeptides{$id}};
#	print "@peptides\n\n";
#}


$"=",";

#print "Unique proteins with more than one peptide support\n";
foreach my $pro(keys %proteinResult){
#	next;
	my %tmp = %{$proteinResult{$pro}};
	my @datasets = keys %tmp;
	my $num = scalar @datasets;
#	print "\nHAHA $num\n$pro\n\n" if ($pro=~/species_sig_Human_immunodeficiency_virus_1_11676_in_9606/);
	if ($num == 1){
		my @supportPeptides = @{$tmp{$datasets[0]}};
		my $numOfSupportPeptides = scalar @supportPeptides;
		my $totalPeptides;
		if ($pro=~/_total_(\d+)/){
			$totalPeptides = $1;
			my $coverage = $numOfSupportPeptides/$totalPeptides;
			print "In dataset $datasets[0] protein $pro has identified peptides @supportPeptides\nwith coverage: $coverage ($numOfSupportPeptides/$totalPeptides)\n" if ($coverage > 0.1);
		}else{
			print "In dataset $datasets[0] protein $pro has identified peptides @supportPeptides\n" if ($numOfSupportPeptides > 1);
		}
	}
}

sub readPeptideExtractFile(){
	my ($id, $file) = @_;
	die "Could not find the file $file\n" unless (-e $file);
}

sub readMzidlibResultFile(){
	my ($id, $file) = @_;
	die "Could not find the file $file\n" unless (-e $file);
	#the parse is based on r541 07 Oct 2014 https://code.google.com/p/mzidentml-lib/source/browse/trunk/src/main/java/uk/ac/liv/mzidlib/converters/MzIdentMLToCSV.java
	#line 209 if(exportOption.equals("exportPSMs"))
	#line 211 header consists of spectrumHeader psmHeader scoreHeader and endPsmHeader
	#line 48 offset 4 spectrumHeader = "Raw data location" + sep + "Spectrum ID" + sep + "Spectrum Title" + sep + "Retention Time (s)" + sep;
	#line 49 the needed sequence and passThreshold are 7th (total 11) and 3rd (total 7) psmHeader = "PSM_ID" + sep + "rank" + sep + "Pass Threshold" + sep + "Calc m/z" + sep + "Exp m/z" + sep + "Charge" + sep + "Sequence" + sep + "Modifications"
	#protein in the endPsmHeader line 54 sep + "proteinacc_start_stop_pre_post_;" + sep + "Is decoy";
	open IN, "$file";
	<IN>;#remove the header
	while (my $line=<IN>){
		chomp $line;
		my @elmts = split(",",$line);
		my $passThreshold =  $elmts[6];
		next unless ($passThreshold eq "true");	#must pass threshold
		my $decoy = lc($elmts[-1]);
		next unless ($decoy eq "false");#not deal with decoy

		my $peptide = $elmts[10]; #11th element # line 439, wrapped by "\""
		$peptide = substr($peptide, 1, (length $peptide)-2); # substr(str, start, length) so last parameter is -2 removing two \"
		
		#line 490-506, wrapped with "\"", separated by ";" if multiple proteins. Each protein is accession_start_end_pre_post, be aware _ may exist in protein accession
		my $protein = $elmts[-2]; #last second
		$protein = substr($protein, 1, (length $protein)-2); # substr(str, start, length) so last parameter is -2 removing two \"
		my @tmp = split(";",$protein);
		foreach my $tmp(@tmp){
			my @aa = split("_",$tmp);
			splice (@aa,-4,4); #remove the last four elements from the array
			my $acc = join("_",@aa);
			$data{$peptide}{$id}{$acc} = 1;#the peptide list is organized by spectrum, so one peptides can be detected several times, so do the proteins => hash
		}
	}	
}

sub usage(){
	print "Usage: perl msSearchResultComparison.pl <peptide list 1 file> <peptide list 1 name> <peptide list 2 file> <peptide list 2 name> [... <peptide list N file> <peptide list N name>]\n";
	print "The script is to compare the peptide lists generated by convertion routine in mzid-lib software. The peptide lists are from the same spectral file searched against different fasta databases\n";
	exit;
}