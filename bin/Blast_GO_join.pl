#!/usr/bin/perl
use warnings;
use strict;

die "Please specify (1) LOC to Prot hash ,  (2) Blast2Go result table \n" unless(@ARGV==2);
print "Starting script\n\nComparing\n\n";
#specify table in ARGV
my $Table = $ARGV[0];
my $Blast2GO=$ARGV[1];
my $outfile="Lanio_GO_hash.tab";

open(my $IN, "<", $Table)   or die "Could not open $Table \n";
open(my $B2G, "<", $Blast2GO)   or die "Could not open $Blast2GO \n";
open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";

my %blast2go_res;
my $head=<$B2G>;

while (my $line2=<$B2G>){
	chomp $line2;
	my @linesplit= split("\t", $line2);
	if ($linesplit[1] && $linesplit[7] && $linesplit[8]){
		$blast2go_res{$linesplit[1]}{"GOs"}="$linesplit[7]";
		$blast2go_res{$linesplit[1]}{"Names"}="$linesplit[8]";
	}
}


my %prot_to_LOC;

while (my $line=<$IN>){
	chomp $line;
	my @linesplit= split("\t", $line);
	#not necessary
	my $LOC=shift(@linesplit);
	$LOC=~ s/gene-//g;
	foreach my $bit (@linesplit){
		$prot_to_LOC{$bit}=$LOC;
	}
	#Main code.
	foreach my $bit (@linesplit){
		if ($blast2go_res{$bit}{"GOs"}){
			my @split_GO=split ("\; ", $blast2go_res{$bit}{"GOs"});
			my @split_Name=split ("\; ", $blast2go_res{$bit}{"Names"});
			my $n=0;
			foreach my $GO (@split_GO){
				my $new=substr $GO, 2;
				my $new_name=substr $split_Name[$n], 2;
				print $outhandle "$LOC\t$new\t$new_name\n";

				$n++;
			}
			
		}
	}
}

`sort Lanio_GO_hash.tab | uniq > Lanio_GO_hash.uniq.tab`;

