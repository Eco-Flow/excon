#!/usr/bin/perl
use warnings;
use strict;


die "Please specify (1)fasta file\n" unless(@ARGV==1);

my $fastafile = $ARGV[0];
my $outfile="$fastafile\.largestIsoform.fa";

open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";


use Bio::SeqIO;
my $seqio = Bio::SeqIO->new(-file => "$fastafile", '-format' => 'Fasta');
my %fastadictionary=();
my @headersplit=();
while (my $seq = $seqio->next_seq){ ## selects one sequence at a time
   	## set variables for THIS sequence
    my $id = $seq->desc;
    #print "here $id\n";
	my $string = $seq->seq;
	my @split=split(/\;/, $id);
	foreach my $name (@split){
		my @obs= split("\=", $name);
		if ($obs[0] eq "gene"){
			$id=$obs[1];
		}
	}
	
	#my @spl2=split("\_", $split[0]);
	#my $waste=pop @spl2;
	my $new_id=$id;
	my $len=length($string);
	#print "length = $len\n";
	if ($fastadictionary{$new_id}){
		my $len_old=length($fastadictionary{$new_id});
		if ($len >= $len_old){
			$fastadictionary{$new_id}=$string;
		}
	}
	else{
		$fastadictionary{$new_id}=$string;
	}
	
}

print "Now print new fasta with one main protein per gene.\n";

foreach my $key ( sort keys %fastadictionary){
	if ($fastadictionary{$key} eq "Sequenceunavailable"){

	}
	else{
		print $outhandle ">$key\n$fastadictionary{$key}\n";
	}
	
}

print "Finished:  input lines, output lines\n";

