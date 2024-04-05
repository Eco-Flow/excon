#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;


die "Please specify (1)fasta file\n" unless(@ARGV==1);

my $fastafile = $ARGV[0];
my $outfile="$fastafile\.noNs";

open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";


use Bio::SeqIO;
my $seqio = Bio::SeqIO->new(-file => "$fastafile", '-format' => 'Fasta');
my %fastadictionary=();
my @headersplit=();


while (my $seq = $seqio->next_seq){ ## selects one sequence at a time
   	## set variables for THIS sequence
    my $id = $seq->display_id;
	my $string = $seq->seq;
	my @split=split(" ", $id);
	print $outhandle ">$id\n$string\n";
}

`mv $fastafile\.noNs $fastafile`;


