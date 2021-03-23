#me auton ton kodika ginetai diagramma mesou megethous exonion se sxesi me ton arithmo exonion pou exei kathe metagrafo gia ola ta exons tou anthropou.

#! /usr/bin/perl

use strict;
use warnings;
use Statistics::R;

open DATA, "<biomart Data exons ids april 2015 no HEADER.txt";

my @rawdata = <DATA>;


MEAN_SIZE_STDERROR(@rawdata);

R_CALC_PLOT();


exit;
##########-----------------------------------------END OF PROGRAM------------------------------##############

sub R_CALC_PLOT
{
    ####----------------------------------enter R for plots------------------------------#######
    #$_=$_;
    # Create a communication bridge with R and start R
    
    #my ($family_name, $family_description) = @_;
    my $R = Statistics::R->new();
    
    # Run simple R commands
    #my $output_file = "$_.pdf";
    $R -> run (q'mydata <- read.table("bashfileforR.txt", header = TRUE)'); # use of temporary file in R


    $R->run(q`pdf("whole_data_set_legent2.pdf")`);
    
    $R->run(
    q 'fit <- lm(log10(mydata$mean_size_of_exons)~log10(mydata$exon_count))',
    q 'rsquared <- summary(fit)$r.squared',
    q 'plot(log10(mydata$exon_count),log10(mydata$mean_size_of_exons), pch=1, main="", xlab="", ylab="" ,sub="")',
    qq 'title(main = "All human exons", xlab=bquote(log[10]~"(number of exons)"), ylab=bquote(log[10]~"(mean size of exons)"))',
    q'arrows(log10(mydata$exon_count),log10(mydata$mean_size_of_exons-mydata$standard_error), log10(mydata$exon_count), log10(mydata$mean_size_of_exons+mydata$standard_error), length=0.05, angle=90, code=3)',
    q'abline(coef=coef(fit))',
    q'calcslope <- coef(fit)[2]',
    q'calcspearman <-cor(log10(mydata$exon_count),log10(mydata$mean_size_of_exons), method = "spearman")',
    q'calcpearson <-cor(log10(mydata$exon_count),log10(mydata$mean_size_of_exons), method = "pearson")',
    q'slope <- paste0("Slope = ",format(round((calcslope), 3), nsmall = 3))',
    q'pastersquared <- paste0("Determ. coef. = ", round(rsquared,2), nsmall=2)',
    q'spearman <- paste0("Spearman = ", round(calcspearman,2), nsmall=2)',
    q'pearson <- paste0("Pearson = ", round(calcpearson,2), nsmall=2)',
    qq'exons <- paste0("Exons = 1005980")',
    qq'transcripts <- paste0("Transripts = 164434")',
    q'legend("top", legend=c(slope, spearman,pearson, pastersquared, exons, transcripts),col = "black", cex = .6, bty = "n")');
    
    my $spearman = $R -> get('spearman');
    my $slope = $R -> get('slope');
    $R->run(q`dev.off()`);
    
    $R->stop();

    }

sub MEAN_SIZE_STDERROR
{
    #=----------------EXON MEAN SIZE PER NUMBER OF EXONS calculations from raw data ----------------------------###
    ### pairno ta fam_data gia na ginou ta epomena
    my %sumsize=();
    my %countexon=();
    my %avg=();
    my %standsum=();
    my %stderror=();
    my %binsums=();
    my %count_transcripts=();
    my %families=();
    
    my @data=@_;
    
    foreach (@data) { #mpainoun ta dedomena. meta upologismos mean megethous ton exonion kai arithmos tous gia kathe transcript ID
    
        chomp;
        my @array = split("\t", $_);
        my $exonsize = ($array[2] - $array[1]); # megethos kathe exoniou
        $sumsize{$array[6]} += $exonsize;
        $countexon{$array[6]}++;
    }
    
    for (sort keys (%countexon))
    {
        $binsums{$countexon{$_}} += $sumsize{$_}; #athrisma ton megethon ton exonion olon ton transcripts pou exount ton idio arithmo exonion
        $count_transcripts{$countexon{$_}}++; #sunolo transcripts pou exoun ton idio arithmo exons.
        
    }
    
    for (sort keys %binsums)
    {
        $avg{$_} = $binsums{$_}/($_*$count_transcripts{$_}); #upologismos mean megethous ton exonion gia ola ta transcript ID pou exoun ton idio arutmo exons.
    }
    
    foreach (@data){ # gia kathe grammi tou arxeiou (ola ta exons) upologizetai to arhroisma toy standard deviation
        chomp;
        my @array = split("\t", $_);
        my $exonsize = ($array[2] - $array[1]);
        $standsum{$countexon{$array[6]}} += ($exonsize - $avg{$countexon{$array[6]}})**2; # dhmiourgia tou athrismatos gia ton tipo tou standard deviation. gia kathe timi transcript dinetai i timi tou hash countexon i opoia apotelei kleidi tou hash standsum.
    }
    
    
    open OUT, ">bashfileforR.txt";
    
    print OUT "exon_count\tmean_size_of_exons\tstandard_error\tstandard_deviation\tnumber_of_transcripts\ttotal_exons\n";
    
    for (sort keys(%standsum)) # here the calculations of mean and standard error for the size of exons per exon number are printed in a file. This file is overwritten for each family ID.
    {
        my $stddev;
        if ($_*$count_transcripts{$_} == 1){
            $stddev = ($standsum{$_}/(($_*$count_transcripts{$_})))**0.5; # to evala LEIPEI to -1 !!!!!!!!!!!!!!!!!!!!!! An se mia oikogeneia uparxei 1 mono transcript me 1 exonio tote to sunolo ton exonion gia bin=1 einai 1.
        }
        else {
            $stddev = ($standsum{$_}/(($_*$count_transcripts{$_})-1))**0.5; # to evala LEIPEI to -1 !!!!!!!!!!!!!!!!!!!!!! An se mia oikogeneia uparxei 1 mono transcript me 1 exonio tote to sunolo ton exonion gia bin=1 einai 1.
            
        }
        $stderror{$_} = $stddev/($_*$count_transcripts{$_})**0.5;
        my $exons = $_*$count_transcripts{$_};
        print OUT "$_\t$avg{$_}\t$stderror{$_}\t$stddev\t$count_transcripts{$_}\t$exons\n"; # thelei >>???
        
    }
}