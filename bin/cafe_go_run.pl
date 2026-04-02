#!/usr/bin/perl
use warnings;
use strict;

my $target      = $ARGV[0];   # e.g. Belonogaster_juncea.pos.txt
my $background  = $ARGV[1];   # e.g. Belonogaster_juncea.BK.txt.uniq
my $go_file     = $ARGV[2];   # OG_GO_format.tsv
my $pval        = $ARGV[3];
my $type        = $ARGV[4];
my $go_max_plot = $ARGV[5];
my $go_algo     = $ARGV[6] // "classic_fisher";

my $lc = `wc -l $target | awk '{print \$1}'`;
chomp $lc;

if ($lc < 10) {
    print "WARNING: Only $lc genes in $target — skipping GO (too few)\n";
    exit 0;
}

if ($lc < 100) {
    print "WARNING: Only $lc genes in $target — GO results may be unreliable\n";
}

print "Running: ChopGO_VTS2.pl -i $target --GO_file $go_file -bg $background -pval $pval -pval_type $type -max_plot $go_max_plot --go_algo $go_algo\n";

my $exit = system("ChopGO_VTS2.pl -i $target --GO_file $go_file -bg $background -pval $pval -pval_type $type -max_plot $go_max_plot --go_algo $go_algo");

exit($exit >> 8);
