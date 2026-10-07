#!/usr/bin/perl
use warnings;
use strict;

print "Please be in folder with focal gff3 file and GO hashes\n\n";

# GFF3 column 9 as key => value, so ID/Parent are found whatever order they come in.
sub parse_attributes {
    my ($col9) = @_;
    my %attr;
    foreach my $bits (split("\;", $col9)){
        my ($key, $value) = split("\=", $bits, 2);
        next if !defined $value;
        $key   =~ s/^\s+|\s+$//g;
        $value =~ s/^\s+|\s+$//g;
        $attr{$key} = $value;
    }
    return %attr;
}

# Orthogroup of a gene: try its gene ID, the gene ID without NCBI's "gene-" prefix,
# then its transcript IDs, since Orthogroups.tsv can name proteins either way.
# Used for both the chromosome and the background tables, so they always agree.
sub find_og {
    my ($og_of, $gene, @trans) = @_;
    (my $stripped = $gene // "") =~ s/^gene-//;
    foreach my $id ($gene, $stripped, @trans){
        return $og_of->{$id} if defined $id && $og_of->{$id};
    }
    return;
}

my $go_algo = $ARGV[0] // "classic_fisher";

my @goes=`ls *.go.txt`;
my @in_gfffile=`ls *.longest.gff`;
my $ortho="Orthogroups.tsv";


#Store orthogroup names hash/
# OrthoFinder writes Orthogroups.tsv with Windows (CRLF) line endings, and chomp only
# removes the \n: the \r left behind was glued to the last gene of the last species
# column on every row, so that gene matched nothing. Strip it from every input read here.
my %orthogroup_hash;
open(my $orthin, "<", $ortho)   or die "Could not open $ortho\n";
my $header=<$orthin>;
$header =~ s/\r?\n\z//;
my @colHeadsplit=split("\t", $header);
while (my $lineOrtho=<$orthin>){
    $lineOrtho =~ s/\r?\n\z//;
    my $n=1;
    my @colsplit=split("\t", $lineOrtho);
    my $OG= shift(@colsplit);
    foreach my $row (@colsplit){
        my @isoform_sp=split(", ", $row);
        my @big_name=split(/\./, $colHeadsplit[$n]);
        my $tidyname=$big_name[0];
        foreach my $iso (@isoform_sp){
            $orthogroup_hash{$tidyname}{$iso}=$OG;
        }
        $n++;
    }
}



my @jobs;   #To store pairs of gff and go files to run

#Check we have the right files for each species:
foreach my $gofile (@goes){
    chomp $gofile;
    my @sp_gp =split(/\./, $gofile);
    my $match=0;
    foreach my $gfffile (@in_gfffile){
        chomp $gfffile;
        my @sp_gff=split(/\./, $gfffile);
        if ($sp_gp[0] eq $sp_gff[0]){
            $match=1;
            push (@jobs, "$gofile $gfffile");
        }
    }
    if ($match){
        #print "$sp_gp[0] is matched\n";
    }
    else{
        print "ERROR: $sp_gp[0] does not have an equivalent gff file\n"
    }
}


#Now run through the jobs and prepare the input files.
foreach my $species (@jobs){
    my @sp=split(/\ /, $species);
    my $go=$sp[0];
    my $gff=$sp[1];
    chomp $gff;
    
    my @go_split =split(/\./, $go);
    my $species_name=$go_split[0];

    my $out="$species_name\.go_r_file.txt";
    open(my $fileout, ">", $out)   or die "Could not open $out\n";

    my $out2="$species_name\.go_r_file.noDuplicates.txt";
    open(my $fileout2, ">", $out2)   or die "Could not open $out2\n";

    # A gene can have more than one transcript line (BRAKER writes GeneMark models
    # as both an mRNA and a transcript), so only write each gene/OG to a scaffold once.
    my %seen_gene_scaffold;
    my %seen_og_scaffold;

    my $og_of = $orthogroup_hash{$species_name} // {};
    my %trans_of_gene;

    open(my $filein, "<", $gff)   or die "Could not open $gff\n";
    while (my $line=<$filein>){
        $line =~ s/\r?\n\z//;
        my @split=split("\t", $line);
        my $gene;
        my $tran;
        my $scaffold=$split[0];
        # Make scaffold name R-safe (R can't handle purely numeric names)
        $scaffold =~ s/^(\d+)$/chr_$1/;
        if ($line =~ /^#/ || scalar(@split) < 9){
            #do nothing
        }
        else{
            # BRAKER/TSEBRA/AGAT annotations write AUGUSTUS models as "transcript"
            # and only GeneMark models as "mRNA", so mRNA lines alone miss most genes.
            if ($split[2] eq "mRNA" || $split[2] eq "transcript"){

                my %attr = parse_attributes($split[8]);
                if ($split[1] eq "AUGUSTUS" || $split[1] eq "maker"){
                    $gene=$attr{"Parent"};
                    $tran=$attr{"ID"};
                }
                else{
                    #Its probably a normal NCBI type:
                    my $fullgene=$attr{"Parent"};
                    my @fullsp=split("\:", $fullgene // "");
                    $gene=$fullsp[-1];
                    # Keep full transcript ID (including rna- prefix) for orthogroup lookup
                    $tran=$attr{"ID"};
                }
                next if !defined $gene;

                # Transcripts of each gene, for the background lookup below
                (my $gene_stripped = $gene) =~ s/^gene-//;
                push @{$trans_of_gene{$gene_stripped}}, $tran if defined $tran;

                # Write to OG duplicates file if found in orthogroups
                if (my $og_id = find_og($og_of, $gene, $tran)){
                    # Sanitise OG ID for R
                    $og_id =~ s/\-/\_/g;
                    $og_id =~ s/\:/\_/g;
                    print $fileout2 "$og_id\t$scaffold\n" unless $seen_og_scaffold{"$og_id\t$scaffold"}++;
                }

                # Sanitise gene ID for R before writing to go_r_file.txt
                my $gene_r = $gene;
                $gene_r =~ s/^gene-//;
                $gene_r =~ s/\-/\_/g if $gene_r;
                $gene_r =~ s/\:/\_/g if $gene_r;
                print $fileout "$gene_r\t$scaffold\n" unless $seen_gene_scaffold{"$gene_r\t$scaffold"}++;
            }
        }
    }

    #Now make a Orthogroup background file:
    open(my $filego, "<", $go)   or die "Could not open $go\n";

    my $out3="$species_name\.go_r_file.noDuplicates_BK.txt";
    open(my $fileout3, ">", $out3)   or die "Could not open $out3\n";

    my %exist_hit;
    while (my $linego=<$filego>){
        $linego =~ s/\r?\n\z//;
        my @splitgo=split("\t", $linego);
        my $go_gene = $splitgo[0];

        # Same lookup as for the chromosome table, with the gene's transcripts from the GFF.
        # The GO file's gene ID alone misses genes that Orthogroups.tsv names by transcript
        # (e.g. rna-XM_...), which the chromosome table still finds.
        (my $go_gene_stripped = $go_gene) =~ s/^gene-//;
        my $found_og = find_og($og_of, $go_gene, @{$trans_of_gene{$go_gene_stripped} // []});

        if ($found_og){
            # Sanitise OG ID for R
            my $og_clean = $found_og;
            $og_clean =~ s/\-/\_/g;
            $og_clean =~ s/\:/\_/g;
            my $hit_key = "$og_clean\t$splitgo[1]";
            if (!$exist_hit{$hit_key}){
                print $fileout3 "$og_clean\t$splitgo[1]\n";
                $exist_hit{$hit_key}="YES";
            }
        }
    }

    close $fileout;
    close $filein;
    print "Now running go chromosome analysis for $species_name\n";

    `sort $species_name\.go_r_file.noDuplicates.txt | uniq > $species_name\.go_r_file.noDuplicates.uniq.txt`;

    my $out4="$species_name\.go_r_file.noDuplicates.noSameScaffold.txt";
    open(my $fileout4, ">", $out4)   or die "Could not open $out4\n";
    open(my $filein4, "<", "$species_name\.go_r_file.noDuplicates.uniq.txt")   or die "Could not open $species_name\.go_r_file.noDuplicates.uniq.txt\n";
    my %hit;
    my %bad;
    while (my $line4=<$filein4>){
        chomp $line4;
        my @split=split("\t", $line4);
        if ($hit{$split[0]}){
            $bad{$split[0]}="duplicate_scaffold";
        }
        $hit{$split[0]}=$split[1];
    }

    foreach my $key ( keys %hit ){
        if ($bad{$key}){
            #don't use it
        }
        else{
            print $fileout4 "$key\t$hit{$key}\n";
        }
    }

    close $fileout4;
    close $filein4;

    my $out5="$species_name\.go_r_file.noDuplicates_BK.noDuplicates.txt";
    open(my $fileout5, ">", $out5)   or die "Could not open $out5\n";
    open(my $filein5, "<", "$species_name\.go_r_file.noDuplicates_BK.txt")   or die "Could not open $species_name\.go_r_file.noDuplicates_BK.txt\n";
    while (my $line5=<$filein5>){
        chomp $line5;
        my @split=split("\t", $line5);
        if ($bad{$split[0]}){
            #do nothing
        }
        else{
            print $fileout5 "$line5\n";
        }
    }
    close $fileout5;
    close $filein5;

    `mkdir Unfiltered_Go_$species_name`;
    `ChopGO_ChromoGoatee.pl -i $species_name\.go_r_file.txt --GO_file $go -sp $species_name --go_algo $go_algo`;
    `mv *_res.tab Unfiltered_Go_$species_name`;
    `mv *_res.tab.pdf Unfiltered_Go_$species_name`;
    `mkdir Filtered_dup_Go_$species_name`;
    `ChopGO_ChromoGoatee.pl -i $species_name\.go_r_file.noDuplicates.noSameScaffold.txt --GO_file $species_name\.go_r_file.noDuplicates_BK.noDuplicates.txt -sp $species_name --go_algo $go_algo`;
    `mv *_res.tab Filtered_dup_Go_$species_name`;
    `mv *_res.tab.pdf Filtered_dup_Go_$species_name`;
}

print "\n\nScript finished\n";
