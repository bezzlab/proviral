#!/usr/bin/perl
use strict;
use constant HOST_PRE => "host_seq_";
use constant HOST_SUF => ".fasta";
use constant VIRAL_PRE => "viral_seqs_for_host_";
use constant VIRAL_SUF => ".fasta";
use constant SIG_VIRUS_PRE => "virus_in_host_";
use constant SIG_SUF => "_sig.fasta";
use constant SIG_SPECIES_PRE => "species_in_host_";

#to combine host sequences and either signatue or all related virus sequences
open IN, "populate_host.tsv" or die "Please provide the list of species as host in the populate_host.tsv file in the order of taxo_id name mnemonic";
my %species;
while (my $line=<IN>){
	chomp $line;
	my ($taxo,$name) = split("\t",$line);
	$name=~s/ /_/g;
	$species{$name} = $taxo;
}


my @todo;
foreach my $name(keys %species){
#	print "$name\t$species{$name}\n";
	my $hostFile = HOST_PRE."$name".HOST_SUF;
	my $viralFile = VIRAL_PRE."$name".VIRAL_SUF;
	my $sigSpeciesFile = SIG_SPECIES_PRE."$species{$name}".SIG_SUF;
	my $sigVirusFile = SIG_VIRUS_PRE."$species{$name}".SIG_SUF;
#	print "$hostFile\n$viralFile\n$sigFile\n\n";
	if(-e $hostFile){
		if (-e $viralFile && -e $sigSpeciesFile && -e $sigVirusFile){
			print "Dealing with host $name\n";
			&output("host_with_sig_virus_$name.fasta","host_with_sig_species_$name.fasta","host_with_all_viral_seqs_$name.fasta",$hostFile,$viralFile,$sigVirusFile,$sigSpeciesFile);
		}
	}else{
		print "No sequence file found for host $name\n";
	}
}

sub output(){
	my ($resultSigVirus,$resultSigSpecies,$resultViral,$hostFile,$viralFile,$sigVirusFile,$sigSpeciesFile) = @_;
#	print "<$hostFile>\t<$viralFile>\t<$sigFile>\n";
	open IN, "$hostFile";
	open VIRAL, ">$resultViral";
	open SIG_VIRUS, ">$resultSigVirus";
	open SIG_SPECIES, ">$resultSigSpecies";
	while (my $line = <IN>){
		print VIRAL "$line";
		print SIG_VIRUS "$line";
		print SIG_SPECIES "$line";
	}
	open IN ,"$viralFile";
	while (my $line = <IN>){
		print VIRAL "$line";
	}
	open IN ,"$sigVirusFile";
	while (my $line = <IN>){
		print SIG_VIRUS "$line";
	}
	open IN, "$sigSpeciesFile";
	while (my $line = <IN>){
		print SIG_SPECIES "$line";
	}
	close IN;
	close VIRAL;
	close SIG_VIRUS;
	close SIG_SPECIES;
}
