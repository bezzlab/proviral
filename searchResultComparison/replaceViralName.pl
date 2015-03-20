use strict;

my %hash;
open IN, "..\\species_in_host_9606_reviewed_sig.fasta" or die "Could not find the file";
#open IN, "..\\species_in_host_9606_all_sig.fasta" or die "Could not find the file";
while(my $line = <IN>){
	if ($line=~/(species_sig_.+)_total_\d+/){
		$hash{$1} = $&;
#		print "$1\n\n$&\n";
#		exit;
	}
}
#my @tmp = keys %hash;
#my $num = scalar @tmp;
#print "find total $num entries\n";

open CSV, "Peptide_result_list_reviewed_human_species_signature_beforeReplacing.csv";
#open CSV, "Peptide_result_list_unreviewed_human_species_signature.csv";
my $count = 0;
#for detailed comment about protein location please refer to searchResultComparison.pl
my $line = <CSV>;
print "$line";
while ($line = <CSV>){
	chomp $line;
	my @elmts = split(",",$line);
	my $protein = $elmts[-2]; #last second
	$protein = substr($protein, 1, (length $protein)-2); # substr(str, start, length) so last parameter is -2 removing two \"
	my @tmp = split(";",$protein);
	for (my $i=0;$i < scalar @tmp;$i++){
		my $tmp = $tmp[$i];
		if ($tmp =~/(species_sig_.+_\d+_in_\d+)/){
			$tmp[$i] = "$`".$hash{$1}."$'";
		}
	}
	$protein = join (";",@tmp);
	$elmts[-2] = "\"$protein\"";
	$line = join(",",@elmts);
	print "$line\n";
}
#print "$count\n";