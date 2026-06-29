#!/usr/bin/perl
use warnings;
use strict;


print "Please be in folder with all the Species Go Summarys\n";

# Which p-value column of the *_TopGo_results_ALL.tab tables to summarise on.
# The correct choice depends on the topGO algorithm, not on user preference:
#   classic_fisher  -> each GO term is tested independently, so a multiple-testing
#                      correction is appropriate. Use column 9 (bonferroni).
#   weight01_*/elim/weight -> the algorithm already accounts for the GO hierarchy,
#                      so adjusting on top is over-conservative. Use column 13
#                      (raw "none" p-value).
# Column indices (0-based): 5=topGOresult(raw) 9=bonferroni 10=BH 12=fdr 13=none
my $go_algo = $ARGV[0] // "classic_fisher";
my $pcol    = ($go_algo eq "classic_fisher") ? 9 : 13;
print "Summarising on column $pcol for go_algo='$go_algo' "
    . ($pcol == 9 ? "(bonferroni-adjusted)\n" : "(raw p-value)\n");

my %go_key;
my %go_names;
my %species_list;
my @gos=`ls *results_ALL.tab`;
my $species;
my $subset;
foreach my $sp (@gos){
    chomp $sp;
    my @split=split(/\./, $sp);
    #my @sp_folder=split("\/", $split[0]);
    $species=$split[0];
    $subset=$split[1];
    #print "$species $subset\n";
    $species_list{$species}="Done";
    open(my $filein, "<", $sp)   or die "Could not open $sp\n";
    my $header=<$filein>;
    while (my $line = <$filein>){
        chomp $line;
        #print "$line\n";
        my @split2=split("\t", $line);
        # $pcol is chosen by algorithm above (raw for hierarchy-aware algorithms,
        # bonferroni for classic_fisher). This value also feeds Count_significant.
        $go_key{$split2[0]}{$subset}{$species}=$split2[$pcol];
        $go_names{$split2[0]}=$split2[1];
        #print "$split2[0] $subset $species $split2[$pcol]\n";
    }
}



#Create a species array:
my @species_array;
foreach my $sp (keys %species_list) {
    push (@species_array, $sp);
}

#Create job type array:
my @job_type_array=("pos","neg");

#Create output files
my $outname="Go_summary_posneg_merged.tsv";
open(my $outhandle, ">", $outname)   or die "Could not open $outname merged\n";
my $outname2="Go_summary_pos.tsv";
open(my $outhandle2, ">", $outname2)   or die "Could not open $outname2 pos\n";
my $outname3="Go_summary_neg.tsv";
open(my $outhandle3, ">", $outname3)   or die "Could not open $outname3 neg\n";


print $outhandle "GO_ID\tGO_term";
print $outhandle2 "GO_ID\tGO_term";
print $outhandle3 "GO_ID\tGO_term";



foreach my $types ( @job_type_array ) {
    foreach my $species ( @species_array ) {
        print $outhandle "\t$species\_$types";
    }
}

foreach my $species ( @species_array ) {
	print $outhandle2 "\t$species";
    print $outhandle3 "\t$species";
}

print $outhandle "\tCount_posSynteny\tCount_negSynteny\n";
print $outhandle2 "\tCount_significant\n";
print $outhandle3 "\tCount_significant\n";



foreach my $goterms (keys %go_key) {
    print $outhandle "$goterms\t$go_names{$goterms}";
    print $outhandle2 "$goterms\t$go_names{$goterms}";
    print $outhandle3 "$goterms\t$go_names{$goterms}";
    
    #All job_type_array
    #My tests_signif count for all sampled:
    my %count_pval;

    #For each test (e.g. topSynteny","botSynteny)
    foreach my $types ( @job_type_array ) {

        foreach my $species ( @species_array ) {
            if (exists $go_key{$goterms}{$types}{$species} ){
                print $outhandle "\t$go_key{$goterms}{$types}{$species}";

                #if p value is less than 0.05, add to counts for GO term:
                if ($go_key{$goterms}{$types}{$species} <= 0.05){
                    if ($count_pval{$goterms}{$types}){
                        $count_pval{$goterms}{$types}++;
                    }
                    else{
                        $count_pval{$goterms}{$types}=1;
                    }

                }

                if($types eq "pos"){
                	print $outhandle2 "\t$go_key{$goterms}{$types}{$species}";
                }
                if($types eq "neg"){
                	print $outhandle3 "\t$go_key{$goterms}{$types}{$species}";
                }    
                
            }
            else{
            	#ALL types:
                print $outhandle "\tNA";

                if($types eq "pos"){
                	print $outhandle2 "\tNA";
                }
                if($types eq "neg"){
                	print $outhandle3 "\tNA";
                }
                
            }
        }

	    #Now print off the significant counts per row:
	    foreach my $types ( @job_type_array ) {
	        if ($count_pval{$goterms}{$types}){
	            print $outhandle "\t$count_pval{$goterms}{$types}";
	        }
	        else{
	            print $outhandle "\t0";
	        }

	    }

	    #Now do for each output table:
	    if($types eq "pos"){
	    	if ($count_pval{$goterms}{$types}){
	    	print $outhandle2 "\t$count_pval{$goterms}{$types}";
	    	}
	    	else{
	    		print $outhandle2 "\t0";
	    	}
	    }
	    if($types eq "neg"){
	    	if ($count_pval{$goterms}{$types}){
	    	print $outhandle3 "\t$count_pval{$goterms}{$types}";
	    	}
	    	else{
	    		print $outhandle3 "\t0";
	    	}
	    }
    }

    print $outhandle "\n";
    print $outhandle2 "\n";
    print $outhandle3 "\n";

}




