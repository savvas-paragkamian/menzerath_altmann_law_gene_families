# script gia upologismo exonion, metagrafon, gonidion, kai mesou megethous exonion gia kathe oikogeneia tou anthropou. Dedomena apo ENSEMBL. Gia ton antropo to script thelei peripou 90 lepta gia na trexei (2.4 GHz Intel Core i5 Haswell, 8gb ram).

#! /usr/bin/perl

use strict;
use warnings;

open DATA, "<final_merge_biomart_exon_transcript_6_2015.txt";
open OUT, ">3_6_15_exons_transcripts_genes_per_family2.txt";
my @rawdata = <DATA>;
my %families=();
my @fam_data=();

foreach (@rawdata) #creation of hash with keys the family IDs.
{
    chomp;
    my @array = split("\t", $_);
    if (scalar(@array) ==10)
    {
        $families{$array[8]} = $array[9];
    }

}

my $i=0;
for my $family (sort keys(%families)) # in this loop there is the calculations and plot export. for each family an array is created which contains all the lines that include this family ID.
{
    foreach (@rawdata)
    {
        chomp;
        if ($_=~ /$family/)
        {
            push(@fam_data, $_); # gia kathe grammi tou arxeiou topothetountai sto arxeio mono autes pou exoun to kleidi (family ID). Etsi gemizei to array me dedomena 1 mono oikogeneias. Sto telos ths loupas adeiazetai to array kai xanagemizei me to kainourgio kleidi.Stis grammes autou tou array ginontai oles oi katametriseis tis oikogeneias
        }
    }

    ### pairno ta fam_data gia na ginou ta epomena
    my %sumsize=();
    my $countexon=0;
    my %avg=();
    my %standsum=();
    my %stderror=();
    my %binsums=();
    my %count_transcripts=();
    my %count_genes=();
    my %fam_number_transcripts=();
    my $exon_count_verification;
    my $total_exon_size=0;

    ### thelo na metrithoun ola ta exonia!! kathe transcript 1 fora.
    
    foreach (@fam_data) #mpainoun ta dedomena. meta upologismos mean megethous ton exonion kai arithmos tous gia kathe transcript ID
    {
        chomp;
        my @array = split("\t", $_);
        
        $total_exon_size += ($array[5] - $array[4]); # megethos(value) kathe exoniou(key)
        $countexon++; #  arithmos pou emfanizontai ta exonia
        $count_transcripts{$array[1]}++; # arithmos pou emfanizetai kathe transcripts(value)=exons pou vriskontai sta dedomena pou trexei h loupa
        $count_genes{$array[2]}++; #katametrisi gonidion
    
    }
    for my $transcript (sort keys(%count_transcripts))
    {
        $exon_count_verification += $count_transcripts{$transcript};
    }
  
    my $avg= $total_exon_size/$countexon;
    
    print OUT $family,"\t", $countexon, "\t", scalar(keys(%count_transcripts)),"\t",scalar(keys(%count_genes)),"\t", $avg,"\t", $families{$family},"\n";
    
    my $leftfam = scalar(keys(%families))-$i; #katametrisi poson oikogeneion apomenoun.
    $i++;
    print $leftfam,"\n";
    @fam_data=(); # adeiasma tou array.IMPORTANT!!!!!!!!!!!
}


