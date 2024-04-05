#!/usr/bin/perl
use warnings;
use strict;


die "Please specify (1)fasta file\n" unless(@ARGV==1);

my $fastafile = $ARGV[0];
my $outfile="$fastafile\.largestIsoform.fa";

open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";

my $outfile2="$fastafile\.summary.txt";
open(my $outhandle2, ">", $outfile2)   or die "Could not open $outfile2 \n";
my $uniq_gene_count=0;
my $total_gene_count=0;

use Bio::SeqIO;
my $seqio = Bio::SeqIO->new(-file => "$fastafile", '-format' => 'Fasta');
my %fastadictionary=();
my @headersplit=();



while (my $seq = $seqio->next_seq){ ## selects one sequence at a time
   	## set variables for THIS sequence
    my $id = $seq->desc;
    chomp $id;
    #print "here $id\n";
	my $string = $seq->seq;
    chomp $string;
    my $new_id;
    my $len=length($string);



	#If entry looks like it has multiple entry info, sep by ;
	if ($id=~ m/\;/g){
		my @split=split(/\;/, $id);
        #Check it has a gene name:
		foreach my $name (@split){
            #print "NAME = $name\n";
			my @obs= split("\=", $name);
			if ($obs[0] eq "gene" || $obs[0] eq "gene_id"){
				$new_id=$obs[1];
                #print "K $new_id $obs[1]\n";
			}
		}
        #Check it has a ID name if not:
        if ($new_id){
            #Do nothing, a name was found.
        }
        else{
            my @split=split(/\;/, $id);
            #Check it has a gene name:
            foreach my $name (@split){
                my @obs= split("\=", $name);
                if ($obs[0] eq "ID"){
                    $new_id=$obs[1];
                }
            }
        }

        #Check if Name is the key, if not found in last two checks:
        if ($new_id){
            #Do nothing, a name was found.
        }
        else{
            my @split=split(/\;/, $id);
            #Check it has a gene name:
            foreach my $name (@split){
                my @obs= split("\=", $name);
                if ($obs[0] eq "Name"){
                    $new_id=$obs[1];
                }
            }
        }

        #Check we got a name, if not die
        if ($new_id){
            #All good
        }
        else{
            die "This gene did not have an ID or a gene entry: \n$id\n\nWe cannot run with GFF files that don't have a ID or name.\n";
        }
	}
    else{
        #IF header is not sep by ;.... probably will be by full stop, e.g. simple AUGUSTUS names.
        my $id2 = $seq->id;
        my @split=split(/\./, $id2);
        $new_id=$split[0];
    }
	
	
	if ($fastadictionary{$new_id}){
        #print "it exists $new_id\n";
		my $len_old=length($fastadictionary{$new_id});
		if ($len >= $len_old){
			$fastadictionary{$new_id}=$string;
		}
	}
	else{
        #print "it didn't $new_id\n";
		$fastadictionary{$new_id}=$string;
        $uniq_gene_count++;
	}
    $total_gene_count++;
	
}

#print off changes in gene to protein number.

print $outhandle2 "$total_gene_count\t$uniq_gene_count\n";

print "For file : $fastafile\n";
print "Total genes\tUnique proteins\n";
print "$total_gene_count\t$uniq_gene_count\n";

#print "Now print new fasta with one main protein per gene.\n";

foreach my $key ( sort keys %fastadictionary){
	if ($fastadictionary{$key} eq "Sequenceunavailable"){

	}
	else{
		print $outhandle ">$key\n$fastadictionary{$key}\n";
	}
	
}

#print "Finished:  input lines, output lines\n";
