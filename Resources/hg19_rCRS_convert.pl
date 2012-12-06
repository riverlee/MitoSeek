#!/usr/bin/perl
#############################################
#Author: Jiang Li
#email: riverlee2008@gmail.com
#Creat Time: 
#Vanderbilt Center for Quantitative Sciences
#############################################
use strict;
use warnings;

#Aim:
#Position and nucleotide mapping between hg19 and rCRS
#perl hg19_rCRS_convert.pl genome/hg19.fasta genome/rCRS.fasta hg19_to_rCRS_table.txt
my ($hg19fasta,$rCRSfasta,$output) = @ARGV;
die "perl $0 hg19fasta rCRSfasta output\n" unless (@ARGV==3);

my $hg19seq=scanreads($hg19fasta);
my $rCRSseq=scanreads($rCRSfasta);

print "Total reads in hg19:",length($hg19seq),"\n";
print "Total reads in rCRS:",length($rCRSseq),"\n";

#print OUT
=convert
1 to 309	0
310	missing
311 to  316	-1
317	missing
318 to 3108	-2
missing->3107
3109 to 16194	-1
16195	missing
16196 to 16571	-2
=cut

open(OUT,">$output") or die $!;
print OUT join "\t",("hg19pos","hg19ref","rCRSpos","rCRSref","isEqual\n");
for(my $i=1;$i<=309;$i++){
    my $hg19ref=substr($hg19seq,$i-1,1);
    my $rCRSref=substr($rCRSseq,$i-1,1);
    my $isEqual=($hg19ref eq $rCRSref)?1:0;
    print OUT join "\t",($i,$hg19ref,$i,$rCRSref,$isEqual);
    print OUT "\n";
}
print OUT join "\t",(310,substr($hg19seq,309,1),"N","","0\n");
for(my $i=311;$i<=316;$i++){
    my $hg19ref=substr($hg19seq,$i-1,1);
    my $rCRSref=substr($rCRSseq,$i-2,1);
    my $isEqual=($hg19ref eq $rCRSref)?1:0;
    print OUT join "\t",($i,$hg19ref,$i-1,$rCRSref,$isEqual);
    print OUT "\n";
}
print OUT join "\t",(317,substr($hg19seq,316,1),"N","","0\n");
for(my $i=318;$i<=3108;$i++){
    my $hg19ref=substr($hg19seq,$i-1,1);
    my $rCRSref=substr($rCRSseq,$i-3,1);
    my $isEqual=($hg19ref eq $rCRSref)?1:0;
    print OUT join "\t",($i,$hg19ref,$i-2,$rCRSref,$isEqual);
    print OUT "\n";
}
print OUT join "\t",("N","","3107",substr($rCRSseq,3106,1),"0\n");

for(my $i=3109;$i<=16194;$i++){
    my $hg19ref=substr($hg19seq,$i-1,1);
    my $rCRSref=substr($rCRSseq,$i-2,1);
    my $isEqual=($hg19ref eq $rCRSref)?1:0;
    print OUT join "\t",($i,$hg19ref,$i-1,$rCRSref,$isEqual);
    print OUT "\n";
}
print OUT join "\t",(16195,substr($hg19seq,16194,1),"N","","0\n");

for(my $i=16196;$i<=16571;$i++){
    my $hg19ref=substr($hg19seq,$i-1,1);
    my $rCRSref=substr($rCRSseq,$i-3,1);
    my $isEqual=($hg19ref eq $rCRSref)?1:0;
    print OUT join "\t",($i,$hg19ref,$i-2,$rCRSref,$isEqual);
    print OUT "\n";
}





sub scanreads{
    my $str="";
    my $in=shift;
    open(IN,$in) or die $!;
    while(<IN>){
        next if (/^>/ || /^$/);
        s/\r|\n//g;
        $str.=$_;
    }
    close IN;
    return $str;
}



