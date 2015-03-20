package Proteomic::Misc;

use strict;
use warnings;
use Exporter;

our @ISA= qw( Exporter );

# these CAN be exported.
our @EXPORT_OK = qw (mw getHydrophobicity ionCalculation calculateMZ);

my $water = 18.01056;
my $H = 1.007825;
my $O = 15.99491;
#http://www.matrixscience.com/help/aa_help.html
#for Pyrrolysine (O) http://www.jbc.org/content/280/44/36962.full.pdf
my %mw = ("A"=>71.03711,"C"=>103.00919,"D"=>115.02694,"E"=>129.04259,"F"=>147.06841,
	  "G"=>57.02146,"H"=>137.05891,"I"=>113.08406,"K"=>128.09496,"L"=>113.08406,
	  "M"=>131.04049,"N"=>114.04293,"P"=>97.05276,"Q"=>128.05858,"R"=>156.10111,
	  "S"=>87.03203,"T"=>101.04768,"V"=>99.06841,"W"=>186.07931,"Y"=>163.06333,
	  "U"=>150.95363,"O"=>237.1456);
my %retentionCoefficient = ("A" => 0.8, "R" => -1.3, "N" => -1.2, "D" => -0.5, "C" => -0.8, "G" => -0.9, "E" => 0.0, "Q" => -0.9, "H" => -1.3, "I" => 8.4, "L" => 9.6, "K" => -1.9, "M" => 5.8, "F" => 10.5, "P" => 0.2, "S" => -0.8, "T" => 0.4, "W" => 11.0, "Y" => 4.0, "V" => 5.0);
my %retentionCoefficientNterm = ("A" => -1.5, "R" => 8.0, "N" => 5.0, "D" => 9.0, "C" => 4.0, "G" => 5.0, "E" => 7.0, "Q" => 1.0, "H" => 4.0, "I" => -8.0, "L" => -9.0, "K" => 4.6, "M" => -5.5, "F" => -7.0, "P" => 4.0, "S" => 5.0, "T" => 5.0, "W" => -4.0, "Y" => -3.0, "V" => -5.5);

sub mw(){
	my $seq = $_[0];
	my $mw=$water;
	my @aa = split("",$seq);
	foreach(@aa){
		$mw+=$mw{$_};
	}
	return $mw;
}

sub getHydrophobicity(){
	my @aa = split("",uc($_[0]));
        my $sumRc = 0;
        foreach(@aa){
            if (exists $retentionCoefficient{$_}){
				$sumRc += $retentionCoefficient{$_};
			}else{
				print "AA $_ not found in the hash\nWhole sequence is $_[0]\n";
			}
        }
        my $correction;
        my $len = scalar @aa;
        if ($len<10){
        	$correction = 1-0.027*(10-$len);
        }elsif($len>20){
        	$correction = 1-0.014*($len-20);
        }else{
        	$correction = 1;
        }
        my $hydrophobicity = $correction*($sumRc+0.42*$retentionCoefficientNterm{$aa[0]}+
                0.22*$retentionCoefficientNterm{$aa[1]}+
                0.05*$retentionCoefficientNterm{$aa[2]});
        $hydrophobicity = $hydrophobicity - 0.3*($hydrophobicity-38) if ($hydrophobicity>38);
        return $hydrophobicity;
}

#calculate theoritical b/y ions m/z values
sub ionCalculation(){
	my $pep = $_[0];

	my @results;
	my $mw=0;
	my @aa = split("",$pep);
	my $count = 0;
	my $ion;
	my $len = length $pep;
	for (my $count = 0; $count < $len; $count++){
		my $aa = $aa[$count];
		$mw+=$mw{$aa}+$water;
		$ion.="$aa";
		$mw-=$water if ($count>0);
		my $value = $mw - $H - $O; #mw = $mw-$H-$O, mz = (mw + $H)/z where z = 1, so mz = $mw - $O
		push (@results,"$pep\tb\t".($count+1)."\t$value");
	}
	#print "$mw\n";
	my @reverseAA = reverse @aa;
	my $tmp_y=$O+3*$H;
	
	for (my $i = 0; $i < $len; $i++){
		my $aa = $reverseAA[$i];
		$tmp_y+=$mw{$aa};
		#my $value = $tmp_y + $H;
		#print "$pep\ty\t".($i+1)."\t$value\n";
		push (@results, "$pep\ty\t".($i+1)."\t$tmp_y");
	}
	return \@results;
}

sub calculateMZ(){
	my $seq = $_[0];
	my $mw = &mw($seq);
	my @mzValues;
	for (my $i=1;$i<5;$i++){
		my $mz = ($mw + $H * $i)/$i;
		push (@mzValues, $mz);
	}
	my $mzValue = join("\t", @mzValues);
	return $mzValue;
}

1;