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

my %LOC_to_transcriptid;
my %transcriptid_to_LOC;
my %transcriptid_to_seq;
my %transcriptid_to_shortid;

while (my $seq = $seqio->next_seq){ ## selects one sequence at a time
    ## set variables for THIS sequence
    my $genename= $seq->id;
    my $id = $seq->desc;
    chomp $id;
    my $string = $seq->seq;
    chomp $string;
    my $new_id;
    my $len=length($string);
    $total_gene_count++;
    #If entry looks like it has multiple entry ENSEMBL info, sep by ;
    if ($id=~ m/\;/g){
        my @split=split(/\;/, $id);

        #Check it has a gene name:
        foreach my $name (@split){
            my @obs= split("\=", $name);
            if ($obs[0] eq "gene"){
                if ($LOC_to_transcriptid{$obs[1]}){
                    my $old=$LOC_to_transcriptid{$obs[1]};
                    $LOC_to_transcriptid{$obs[1]}="$old\t$genename";
		    #print "$obs[1] $LOC_to_transcriptid{$obs[1]}\n";
                }
                else{
                    $LOC_to_transcriptid{$obs[1]}=$genename;
                }
                #should be unique
                if ($transcriptid_to_LOC{$genename}){
                    die "This should not happen, means a transcript ID was found twice.\n";
                }
                else{
		    #print "$genename\t$string\n";
                    $transcriptid_to_seq{$genename}=$string;
                    $transcriptid_to_LOC{$genename}=$obs[1];
                    my @tmp_split=split(/\./, $genename);
                    $transcriptid_to_shortid{$genename}=$tmp_split[0];
                }
                
            }
        }
    }
}


foreach my $key ( sort keys %LOC_to_transcriptid){
    my @vals=split("\t", $LOC_to_transcriptid{$key});
    my $longest=0;
    my $best;
    my $bestseq;
    $uniq_gene_count++;
    #print "$key $LOC_to_transcriptid{$key}\n";
    foreach my $trans (@vals){
        my $length=length($transcriptid_to_seq{$trans});
        if ($length >= $longest){
            $best=$trans;
            $longest=length($transcriptid_to_seq{$trans});
            $bestseq=$transcriptid_to_seq{$trans};
        }
        else{
            #do nothing;
        }
    }
    print $outhandle ">$transcriptid_to_LOC{$best}\n$bestseq\n";
	#$transcriptid_to_shortid{$best}
}


print "For NCBI file : $fastafile\n";
print "Total genes\tUnique proteins\n";
print "$total_gene_count\t$uniq_gene_count\n";


