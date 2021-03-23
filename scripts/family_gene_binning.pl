# test script to find families and put their data to array
#! /usr/bin/perl

use strict;
use warnings;


open DATA, "</3_6_15_exons_transcripts_genes_per_family.txt";
open OUT, ">binning_genes_per_family.txt";
my @rawdata = <DATA>;
my %families=();
my @fam_data=();
my %binning_fams=();

#print "family\ttotal_exons\ttotal_transcripts\tmean_size_of_exons\tfamily_discription\n";


foreach (@rawdata) #creation of hash with keys the family IDs.
{
    chomp;
    my @array = split("\t", $_);
    $families{$array[0]}=$array[3];
}

for my $family (sort keys(%families))
{
    $binning_fams{$families{$family}}++;
}

print OUT "number of genes\tnumber of families\n";
for my $bin (sort keys(%binning_fams))
{
    print OUT $bin, "\t", $binning_fams{$bin}, "\n";
}
