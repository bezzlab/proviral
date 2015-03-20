#!/usr/bin/perl
use strict;
use Proteomic::Misc qw(getHydrophobicity ionCalculation calculateMZ);
use Proviral qw(checkSource checkHost);

#Virus
#Inputs: swissprot virus dat
#read in the file
#	parse ID, AC, OH, OS, OX, SQ lines
#	insert infection table (host, virus)
#	do protein digest
my $source = "Trembl";
$source = "Swissprot";
my $version = "Mar 2015";

my $dataFolder = "e:\\proviral\\datafiles\\virus\\$source\\$version\\";
$dataFolder =~s/ /_/g;
unless (-d $dataFolder){
	system ("mkdir $dataFolder");
	print "The folder <$dataFolder> has been created\n\n";
}
print "The data files will be saved in the folder <$dataFolder>\n";

my $datFile = "virus_sample.dat";
$datFile = "e:\\uniprot data\\$version\\uniprot_sprot_viruses.dat";
#$datFile = "e:\\uniprot data\\$version\\uniprot_trembl_viruses.dat";

my %species;
my %hosts;
#the same file used to populate Host table which is used to only analyze viruses which infect the required hosts
open IN, "populate_host.tsv" or die "Please provide the list of species as host in the populate_host.tsv file in the order of taxo_id name mnemonic";
while (my $line=<IN>){
	chomp $line;
	my ($taxo,$name) = split("\t",$line);
	my $filename = "${dataFolder}viral_seqs_for_host_${name}_${source}_$version.fasta";
	$filename=~s/\s/_/g;
	if(-e $filename){
		print "file $filename already exists, now deleting\n";
		system ("del $filename");
	}
	$hosts{$taxo} = $filename;
}

open INF, ">${dataFolder}populate_infection_virus_${source}_$version.tsv";
open PRO, ">${dataFolder}populate_protein_virus_${source}_$version.tsv";
open EXISTANCE, ">${dataFolder}populate_peptide_existance_virus_${source}_$version.tsv";
open PEP, ">${dataFolder}populate_peptide_virus_${source}_$version.tsv";
open ION, ">${dataFolder}populate_product_virus_${source}_$version.tsv";
open VIRUS, ">${dataFolder}populate_virus_virus_${source}_$version.tsv";

my %doneVirus; # keys are taxo id, value is another hash with predefined keys including name, virus, subtype, strain
#my $count = 0;
open IN, "$datFile" or die "Can not find the specified file $datFile";
my @lines=();
while(my $line=<IN>){
	chomp $line;
	if($line eq "\/\/"){
 		&dealOneProtein(\@lines);
#		$count++;
		@lines=();
	}else{
		push(@lines,$line);
	}
}
#print "Count: $count\n";

#foreach my $species(sort {$a cmp $b} keys %species){
#	print "$species\t$species{$species}\n";
#}
#one protein section starts with ID line and finishes with //
sub dealOneProtein(){
	my @lines = @{$_[0]};
	my $id;
	my $ac;
	my $taxo;
	my $species;
	my $seq = "";
	my @osLines;
	my @ohLines;
	#collect all needed information for the current protein entry
	for(my $i=0;$i<scalar @lines;$i++){
		my $line = $lines[$i];
		if($line =~/^ID\s+(\S+)/){
			$id = $1;
			#the following 4 lines of codes for manual manipulation if stopped unexpected 
#			return if ($id ne "14310_ARATH");
#			return if ($id ne "CARM1_RAT");
#			return if ($id ne "PLEC_MOUSE");
#			print "@lines\n";
#			exit;
		}elsif($line=~/^AC\s+/){
			#Researchers who wish to cite entries in their publications should always cite the first accession number. This is commonly referred to as the 'primary accession number'.
			#if > 0, already assigned
			($ac)=split(";",$') unless ((length $ac) > 0);
		#The sequence data line has a line code consisting of two blanks rather than the two-letter codes used until now. The sequence counts 60 amino acids per line, in groups of 10 amino acids, beginning in position 6 of the line.
		}elsif($line=~/^\s+/){
			$seq.=$line;
		#The OX (Organism taxonomy cross-reference) line is used to indicate the identifier of a specific organism in a taxonomic database. The format of the OX line is:
		#OX   Taxonomy_database_Qualifier=Taxonomic code;
		#Currently the cross-references are made to the taxonomy database of NCBI, which is associated with the qualifier 'TaxID' and a taxonomic code.		
		}elsif($line=~/^OX\s+NCBI_TaxID=(\d+)/){
			$taxo = $1;
		#The OS (Organism Species) line specifies the organism which was the source of the stored sequence. In the rare case where all the species information will not fit on a single line, more than one OS line is used. The last OS line is terminated by a period.
		#The species designation consists, in most cases, of the Latin genus and species designation followed by the English name (in parentheses). For viruses, only the common English name is given.
		#The names (official name, common name, synonym) concerning one species are cut across lines when they do not fit into a single line:		
		#OS   Rous sarcoma virus (strain Schmidt-Ruppin A) (RSV-SRA) (Avian leukosis
		#OS   virus-RSA).
		}elsif($line=~/^OS\s+/){
			push (@osLines,$line);
			while(1){
				$i++;
				my $newLine = $lines[$i];
				last unless ($newLine =~/^OS\s+/);
				push (@osLines, $newLine);
			}
		#The OH (Organism Host) line is optional and appears only in viral entries. It indicates the host organism(s) that are susceptible to be infected by a virus.
		#A virus being an inert particle outside its hosts, the virion has neither metabolism, nor any replication capability, nor autonomous evolution. Identifying the host organism(s) is therefore essential, because features like virus-cell interactions and posttranslational modifications depend mostly on the host.
		#The format of the OH line is:
		#OH   NCBI_TaxID=TaxID; HostName.

		}elsif($line=~/^OH\s+/){
			push (@ohLines,$line);
			while(1){
				$i++;
				my $newLine = $lines[$i];
				last unless ($newLine =~/^OH\s+/);
				push (@ohLines, $newLine);
			}
		}
	}
	#parse the collected information
#	print "$ac\t$id\t$taxo\n";
	$species = join (" ",@osLines);
	$species=~s/OS\s+//g; #multiple lines, so g is needed to remove all OS tags
	$species = substr($species,0,(length $species) - 1); #remove the . at the end
	#print "$species\n";
	my $len = length $species;
	print "WARNING: $species length $len > 255\n" if ($len > 254); 

	my @hosts;
	foreach my $line(@ohLines){
		if ($line=~/NCBI_TaxID=(\d+)/){
			push (@hosts, $1) if (exists $hosts{$1});
		}
	}
#	$species{$species} = $taxo if (scalar @hosts > 0);
#	$species{$species} = $taxo;
	if (scalar @hosts > 0){
		$seq =~s/\s//g;
		unless(exists $doneVirus{$taxo}){#populated in the parseOS subroutine
			&parseOS($species, $taxo);
			#based on the assumption that all possible hosts are given for the virus, i.e. every time the host list is identical for the same virus
			foreach my $host(@hosts){
				my $filename = $hosts{$host};
				open OUT, ">>$filename";
				print OUT ">$ac|$id\n$seq\n";
				close OUT;
				print INF "$host\t$taxo\n";
			}
		}
		print PRO "$ac\t$id\t\\N\t$taxo\t$source\t$seq\n";
		#the codes below are the same codes as in host_tables_prep.pl
		#the reason not to extract into subroutine is due to the need to output into different result files
		my @peptides = split (/(?<=[KR])(?=[^P])/,$seq); #alternative RE /(?!P)(?<=[RK])/
		my %peptidesCount;
		foreach my $peptide (@peptides) {
			$peptidesCount{$peptide}++ if (length $peptide <= 30 && length $peptide >= 3);
		}
		foreach my $peptide (sort {$a cmp $b} keys %peptidesCount){
			next if ($peptide=~/[UOBZJX]/);
			print EXISTANCE "$peptide\t$ac\t$peptidesCount{$peptide}\n";
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
#	print "Sequence: $seq\n";
}

foreach my $virusTaxo (keys %doneVirus){
	my %hash = %{$doneVirus{$virusTaxo}};
	print VIRUS "$virusTaxo\t$hash{'name'}\t$hash{'virus'}\t$hash{'subtype'}\t$hash{'strain'}\n";
}

#http://web.expasy.org/docs/userman.html#OS_line
#The species designation consists, in most cases, of the Latin genus and species designation followed by the English name (in parentheses). For viruses, only the common English name is given.
sub parseOS(){
	my $desc = $_[0];
	my $taxo = $_[1];
	my $len = length $desc;
	my @segments;
	my $level = 0;
	my @tmp = ();
	my $name = "";
	#only consider one level of (), i.e. any additional () in the () will be kept intact
	for (my $i = 0; $i < $len; $i++){
		my $curr = substr ($desc, $i, 1);
		if ($curr eq "("){
			if ($level == 0){
				my $tmp  = join ("", @tmp);
				$tmp = trim($tmp);
				$name = $tmp if($name eq "");
				push (@segments, $tmp) if (length $tmp > 0);
				@tmp=();
			}else{
				push (@tmp, $curr);
			}
			$level++;
		}elsif($curr eq ")"){
			$level --;
			if ($level == 0){
				my $tmp  = join ("", @tmp);
				$tmp = trim($tmp);
				push (@segments, $tmp);
				@tmp=();
			}else{
				push (@tmp, $curr);
			}
		}else{
			push (@tmp, $curr);
		}
	}
	my $tmp  = join ("", @tmp);
	$tmp = trim($tmp);
	$name = $tmp if (length $name == 0); #if no (), then at this point $name is empty
	push (@segments, $tmp) if (length $tmp > 0);
	
	my @shortNameCandidates;
	my @strainCandidates;
	foreach my $segment(@segments){
		if($segment=~/(strain|isolate)/i){
			push (@strainCandidates, $segment) ;
#		if ($segment=~/^\S+$/){#normally no space in the short name/abbreviation
		}elsif (length $segment<9){ #as suggest by the name "short name", the length normally is short, 9 is picked by testing the different values also due to the length of strings "strain" and "isolate"
			push (@shortNameCandidates, $segment); 
#		}elsif($segment=~/(subtype|genotype|serotype|strain|isolate)/i){
		}
	}
#species type distribution codes	
#	my $a = scalar @shortNameCandidates;
#	my $b = scalar @strainCandidates;
#	my $c = scalar @remaining;
#	$count{"$a-$b"}++;
#	push(@{$entries{"$a-$b"}},"original: <$desc>\nname: <$name>\n");

	my $virus = "";
	my $subtype = "";
	my $strain = "";
	if ($name=~/(subtype|genotype|serotype)/i){
		$virus = $`;
		$virus = trim($virus);
		$subtype = "$&$'";
	}
	#give preference to the segment which contain only strain info
	for (my $i=0;$i<scalar @strainCandidates;$i++){
		my $abc = $strainCandidates[$i];
		if ($abc=~/^(strain|isolate)/i){
			$strain = $abc;
			last;
		}
	}
	#if no such segment found, any strain segment will do
	$strain = $strainCandidates[0] if (length $strain == 0 && scalar @strainCandidates > 0);
	
	#AFV-1, AFV-2 treat as AFV subtype 1 and 2
	#Acidianus filamentous virus 1 (isolate United States/Yellowstone) (AFV-1)
	if(length $subtype == 0){
		my $newName;
		if ($name=~/(virus)\s+(\d+)$/){
			my $serial = $2;
			$newName = "$`$1";
			foreach my $tmp(@shortNameCandidates){
				if ($tmp=~/\d+$/){
					if ($serial == $&){#if the same number appears in both the full name and short name, it should be the subtype
						$subtype = $serial;
						$virus = $newName;
						last;
					}
				}
			}
		}
	}
	
	$virus = $1 if ($name =~/^(.+\s+virus)(\s+\d+)?$/i && length $virus == 0);
	$virus = $name if (length $virus == 0);
#	print "original: <$desc>\nname: <$name>\nvirus: <$virus>\nsubtype: <$subtype>\nstrain <$strain>\n\n";
#	print "original: <$desc>\nvirus: <$virus>\nsubtype: <$subtype>\nstrain <$strain>\n\n";
	$doneVirus{$taxo}{'name'} = $desc;
	$doneVirus{$taxo}{'virus'} = $virus;
	$doneVirus{$taxo}{'subtype'} = $subtype;
	$doneVirus{$taxo}{'strain'} = $strain;
}

sub trim(){
	my $s = shift; 
	$s =~ s/^\s+|\s+$//g; 
	return $s 
}