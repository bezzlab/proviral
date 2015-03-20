#!/usr/bin/perl
use strict;
opendir DIR, ".";
my @files = readdir DIR;
closedir DIR;

open SP,">combinedSwissprotHost.fasta";
open TR,">combinedTremblHost.fasta";
foreach my $file(@files){
	if ($file=~/^swissprot\s\d+.+\.fasta$/i){
		print "sw <$file>\n";
		&printToFile($file,"sp");
	}elsif($file=~/^trembl\s\d+.+\.fasta$/i){
		print "tr <$file>\n";
		&printToFile($file,"tr");
	}
}

sub printToFile(){
	my ($inFile,$type) = @_;
	open IN,"$inFile";
	while(my $line=<IN>){
		if ($type eq "sp"){
			print SP "$line";
		}else{
			print TR "$line";
		}
	}
}