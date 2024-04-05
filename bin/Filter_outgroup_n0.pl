#!/usr/bin/perl
use warnings;
use strict;


die "Please specify (1) N0.tsv (2) list of excluded species\n" unless(@ARGV==2);

my $n0 = $ARGV[0];
my $excluded = $ARGV[1];

my $outfile="$n0\.filt";
open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";

#store neg abmes in hash
open(my $exhandle, "<", $excluded)   or die "Could not open $excluded \n";
my $line_ex=<$exhandle>;
chomp $line_ex;
my @split_ex=split(/\ /, $line_ex);
my %negatives;
foreach my $neg (@split_ex){
    print "\<$neg\>";
    $negatives{$neg}="hit";
}

#open tab file and target cols with negatives.
open(my $inhandle, "<", $n0)   or die "Could not open $n0 \n";

my $head=<$inhandle>;
chomp $head;
my @h_split=split("\t", $head);

my $n=0;
my %cols_to_remove;
my @store;
foreach my $col (@h_split){
    
    if ($negatives{$col}){
        $cols_to_remove{$n}="TBD";
        print "$col $n \n";
    }
    else{
        
        push (@store, $col);
    }
    $n++;
    
    
}
my $joint=join("\t", @store);
print $outhandle "$joint\n";


while (my $line2=<$inhandle>){
    chomp $line2;
    my @split=split("\t", $line2);
    my $n=0;
    my @stor2;
    foreach my $col (@split){
        if ($cols_to_remove{$n}){
            #do nothing
        }
        else{
            push (@stor2, $col);
            #print $outhandle "\t$col";
        }
        $n++;
    }
    my $join_here=join("\t", @stor2);
    print $outhandle "$join_here\n";
}


