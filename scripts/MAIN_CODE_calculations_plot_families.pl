#script gia ton upologismo gia kathe oikogeneia ton gonidion tou anthropou: mesou megethous exonion ta opoia anoikoun sta metagrafa pou exoun ton idio arithmo exonion. Epeita diagramma stin R me Y axona meso megethos exonion kai X axona arithmo exonion metagrafon. O kodikas trexei gia 1200 oikogeneies se 10 lepta peripou (2.4 GHz Intel Core i5 Haswell, 8gb ram).


#! /usr/bin/perl

use strict;
use warnings;
use Statistics::R;

open DATA, "<final_merge_biomart_exon_transcript_6_2015.txt";
open FAM, "<3_6_15_exons_transcripts_genes_per_family.txt";
my @fam_selection = <FAM>;
my @rawdata = <DATA>;
my %families=();
my %families_description=();
my %family_exons=();
my %family_transcripts=();
my %transcript_number=();
my @fam_data=();
our $family;

foreach (@fam_selection) #creation of hash with keys the family IDs.
{
    chomp;
    my @array = split("\t", $_);
    $families{$array[0]} = $array[1];
    $transcript_number{$array[0]} = $array[2];
    $families_description{$array[0]} = $array[5];
    $family_exons{$array[0]}=$array[1];
    $family_transcripts{$array[0]}=$array[2];
    }

for $family (sort keys(%families)) # in this loop there is the calculations and plot export. for each family an array is created which contains all the lines that include this family ID.
{
    #print $_,"\n";
    if (($families{$family} > 100) && ($transcript_number{$family} >= 10) ) # epilogh oikogeneion me vasi ton arithmo exon pou exoun KAI ton arithmo exonion.
    {
    
        foreach my $line (@rawdata)
        {
            chomp($line);
            if ($line =~ /$family/)
            {
                push(@fam_data, $line); # the array with lines that include the key(family ID). With each iteration the array is empted and new lines are imported for the specific key. This array is the source for the rest of the program for the calculations and plotting
                #print "$_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
                #print "nai\n";
            }
        }
    
        MEAN_SIZE_STDERROR(@fam_data);
        
        R_CALC_PLOT($family, \%families_description, \%family_exons, \%family_transcripts);
    }
    @fam_data=(); # ADEIASMA ARRAY!!!!
}

exit;
##########-----------------------------------------END OF PROGRAM------------------------------##############

sub R_CALC_PLOT # regression analysis and plot in R for every family
{
    ####----------------------------------enter R for plots------------------------------#######
    
    open STATISTIC, ">>family_slope_spearman.txt";
    my ($family_name, $family_description, $family_exons, $family_transcripts) = @_;
    my $R = Statistics::R->new();
    
    $R -> run (q'mydata <- read.table("bashfileforR.txt", header = TRUE)'); # use of temporary file in R
    $R->run(qq`pdf("$family.pdf")`);
    $R->run(
    q 'fit <- lm(log10(mydata$mean_size_of_exons)~log10(mydata$exon_count))', # regression analysis
    q 'rsquared <- summary(fit)$r.squared', # coefficient of determination
    q 'plot(log10(mydata$exon_count),log10(mydata$mean_size_of_exons), pch=1, main="", xlab="", ylab="" ,sub="")', # dimiourgia plot me dipli logarithmiki klimaka ton metavlitvn
    qq 'title(main = "$family_name",sub = "${$family_description}{$family_name}", xlab=bquote(log[10]~"(number of exons)"), ylab=bquote(log[10]~"(mean size of exons)"))', # titloi kai axones
    q'arrows(log10(mydata$exon_count),log10(mydata$mean_size_of_exons-mydata$standard_error), log10(mydata$exon_count), log10(mydata$mean_size_of_exons+mydata$standard_error), length=0.05, angle=90, code=3)', # topothetish error bars. einai ena mikro hack pou xrisimopoiountai arrows.
    q'abline(coef=coef(fit))', # eutheia grammikis pallidromisis
    q'calcslope <- coef(fit)[2]', # klisi eutheias
    q'calcspearman <-cor(log10(mydata$exon_count),log10(mydata$mean_size_of_exons), method = "spearman")', # upologismos Spearman coefficient
    q'calcpearson <-cor(log10(mydata$exon_count),log10(mydata$mean_size_of_exons), method = "pearson")', #upologismos Pearson coefficient
    q'slope <- paste0("Slope = ",format(round((calcslope), 3), nsmall = 3))', # ektiposi
    q'pastersquared <- paste0("Determ. coef. = ", round(rsquared,2), nsmall=2)',
    q'spearman <- paste0("Spearman = ", round(calcspearman,2), nsmall=2)', # ektiposi
    q'pearson <- paste0("Pearson = ", round(calcpearson,2), nsmall=2)',# ektiposi
    qq'exons <- paste0("Exons = ","${$family_exons}{$family_name}")',# ektiposi
    qq'transcripts <- paste0("Transripts = ", "${$family_transcripts}{$family_name}")',# ektiposi
    q'legend("top", legend=c(slope, spearman,pearson, pastersquared, exons, transcripts),col = "black", cex = .6, bty = "n")');# topothetish ton prohgoumenon ektiposeon sto upomnhma tou plot
    
    my $spearman = $R -> get('calcspearman'); #krataei metavliti
    my $pearson = $R -> get('calcpearson');
    my $slope = $R -> get('calcslope');
    $R->run(q`dev.off()`);

    $R->stop();
    
    print STATISTIC $family, "\t", $spearman, "\t", $pearson, "\t", $slope, "\n"; # ektiposi gia kathe oikogeneia ta statistika pou tin xaraktirizoun

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
    my @fam_data=@_;
    
    foreach (@fam_data) #mpainoun ta dedomena. meta upologismos mean megethous ton exonion kai arithmos tous gia kathe transcript ID
    {
        chomp;
        my @array = split("\t", $_);
        my $exonsize = ($array[5] - $array[4]); # megethos kathe exoniou
        $sumsize{$array[1]} += $exonsize; # athrisma megethous exonion
        $countexon{$array[1]}++; # arithmos exonion ana transcript
        #printf $array[0],"\t", $array[1],"\t", $array[2],"\t", $array[3], "\n";
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
    
    foreach (@fam_data){ # gia kathe grammi tou arxeiou (ola ta exons) upologizetai to arhroisma toy standard deviation
        chomp;
        my @array = split("\t", $_);
        my $exonsize = ($array[5] - $array[4]);
        $standsum{$countexon{$array[1]}} += ($exonsize - $avg{$countexon{$array[1]}})**2; # dhmiourgia tou athrismatos gia ton tipo tou standard deviation. gia kathe timi transcript dinetai i timi tou hash countexon i opoia apotelei kleidi tou hash standsum.
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
            $stddev = ($standsum{$_}/(($_*$count_transcripts{$_})-1))**0.5; #
        }
        $stderror{$_} = $stddev/($_*$count_transcripts{$_})**0.5;
        my $exons = $_*$count_transcripts{$_};
        print OUT "$_\t$avg{$_}\t$stderror{$_}\t$stddev\t$count_transcripts{$_}\t$exons\n";
    }
}
