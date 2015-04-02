#!/usr/bin/perl
use strict;

my %data;

my $suffix = "_in_host_9606_all_sig.fasta";
#my $suffix = "_in_host_9606_reviewed_sig.fasta";
&readin("combined");
&readin("species");
print "\n";

my %species = %{$data{"species"}};
my %combined = %{$data{"combined"}};

foreach my $id (keys %species){
	if (exists $combined{$id}){
		my $seqSpecies = $species{$id};
		my $seqCombined = $combined{$id};
		if ($seqSpecies ne $seqCombined){
			print "Difference for $id\n";
			print "combined: <$seqCombined>\n";
			print "species : <$seqSpecies>\n";
		}
		delete $combined{$id};
	}else{
		print "$id not found in combined";
	}
}

print "Only in combined\n";
foreach my $id (keys %combined){
	print "$id\n$combined{$id}\n\n";
}

sub readin(){
	my $type = $_[0];
	my $filename="$type$suffix";
	print "reading $filename\n";
	open IN, "$filename";
	my $header = "";
	my $seq = "";
	while(my $line=<IN>){
		chomp($line);
		if ($line=~/^>/){
			if (length $seq > 0){
				$data{"$type"}{$header} = $seq;
			}
			$header = $line;
			$seq = "";
		}else{
			$seq.=$line;
		}
	}
	$data{"$type"}{$header} = $seq;
}
