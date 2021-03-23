#enosi 2 arxeion keimenou pou exoun koini stili

#! /usr/bin/perl
use warnings;
use strict;


open EXA1, "<biomart Data exons ids april 2015.txt";
open EXA2, "<uniq_biomart_families.txt";
open OUT, ">merge_biomart_exon_family.txt";

my %hash = ();
my @coordinates = <EXA1>;
my @ids = <EXA2>;


foreach (@coordinates)
{
    chomp;
    my ($char, $start, $end, $rank, $tr_count, $exon_id, $transcript_id, $gene_id)=split("\t", $_);
    my $value = "$gene_id\t$char\t$start\t$end\t$rank\t$tr_count";
    my $pair = "$exon_id\t$transcript_id"; #monadikotita grammis.
    $hash{$pair} = $value; # topothetish stoixeion sto hash me kleidia stili 0. An kapoio stoixeio kleidiou iparxei 1< fores to kleidi apodidei perissoteres times.
}
foreach (@ids)
{
    chomp;
    my ($gene_id2, $transcript_id2, $exon_id2, $family_id, $family_discription)=split("\t", $_);
    
    my $pair = "$exon_id2\t$transcript_id2";
    if (exists $hash{$pair})
    {
        print OUT $exon_id2, "\t",$transcript_id2, "\t",$hash{$pair},"\t", $family_id, "\t", $family_discription, "\n";
    }
    else
    {
        print "$gene_id2\n";
    }
    
}


