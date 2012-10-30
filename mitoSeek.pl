#!/usr/bin/perl
#############################################
#Author: Jiang Li
#email: riverlee2008@gmail.com
#Creat Time: Tue 23 Oct 2012 01:37:54 PM CDT 
#Vanderbilt Center for Quantitative Sciences
#############################################
use strict;
use warnings;
use List::Util qw(sum);

our $samtools="samtools";

my $debug=1;
print $ENV{PATH};

_check_samtools();


_get_mitochondrial_bed("1032QC01_sorted.bam");




#
sub _mito_qc_stat{
    
}


#Use samtools view -H read the header, and prepare the mithochondrial bed file.
sub _get_mitochondrial_bed{
    my ($in) = @_;
    open(IN,"samtools view -H $in|") or die $!;
    my $flag=0;
    my $m="";
    my $len="";
    while(<IN>){
        if(/SN:M|SN:chrM/i){
            if(/SN:(.*?)\s+LN:(\d+)/){
                $m=$1;
                $len=$2;
                $flag=$1;
                last;
            }
        }
    }
    unless($flag){
        exit(1);
    }
    open(OUT,">tmp/m.bed");
    print OUT join "\t",($m,1,$len);
    close OUT;
}


#Determine heteroplasmy mutations from a basecall format file with given parameters
#Parameters:
#$inbase  input file of allele count parsed from pileup file
#$hp      (-hp) Heteroplasmy threshold using percent alternative allele observed
#$ha      (-ha) Heteroplasmy threshold using allele observed
#$isall   (-A)  0/1, 1 denotes that the total read count is the total allele count of all allele observed,
#         while 0 indicates the total read count is the sum of major and minor allele counts, default=0
#$sb      (-sb) Remove all sites with strand bias score in the top %, default = 10
sub _determine_heteroplasmy{
    my($inbase,$hp,$ha,$isall,$sb) = @_;
    open(IN,$inbase) or die $!;
    <IN>; #skip the header
    my @result;
    my %sb;
    while(<IN>){
        s/\r|\n//g;
        my($chr,$loc,$ref,$forward_A,$forward_T,$forward_C,$forward_G,$reverse_A,$reverse_T,$reverse_C,$reverse_G) = split "\t";
        # Get the major allele and minor allele.
        my %atcg;
        $atcg{A}=$forward_A+$reverse_A;
        $atcg{T}=$forward_T+$reverse_T;
        $atcg{C}=$forward_C+$reverse_C;
        $atcg{G}=$forward_G+$reverse_G;
        my @atcg = sort {$atcg{$b}<=>$atcg{$a}} keys %atcg;
        my $major_allele=$atcg[0];
        my $minor_allele=$atcg[1];
        
        #Calculate heteroplasmy ratio
        my $heteroplasmy=0;
        if($isall){
            $heteroplasmy = $atcg{$minor_allele}/($atcg{A}+$atcg{C}+$atcg{T}+$atcg{G});
        }else{
            $heteroplasmy = $atcg{$minor_allele}/($atcg{$major_allele}+$atcg{$minor_allele});
        }
        
        #Calculate 95% confidence interval
        my ($lower,$upper) = (0,0);
        if ($heteroplasmy>0){
            my $n=0;
            if($isall){
                $n=$atcg{A}+$atcg{C}+$atcg{T}+$atcg{G};
            }else{
                $n=$atcg{$major_allele}+$atcg{$minor_allele};
            }
            ($lower,$upper) = _ci($heteroplasmy,$n);
        }
        
        #Stat values could be 
        #0 (not a heteroplasmy mutation)
        #1 (show heteroplasmy, but not pass the cutoff)
        #2 (heteroplasmy mutation pass cutoff and not have strong strand bias)
        #3 (heteroplasmy mutation pass cutoff but show strong strand bias)
        my $stat=0;
        if($heteroplasmy>0){
            #determine 1 or 2 at this stage (for value 3, need to first calculate strand bias across all the site, and then modify those with stat=2, only if when $sb is not equal to 0)
            if($heteroplasmy>$hp/100 && $atcg{$minor_allele}>$ha){
                $stat=2;
            }else{
                $stat=1;
            }
        }
        
        #If provided strand bias cutoff, need to re-loop the @result, change state=2 to 3 when necessary.
        push @result,[$chr,$loc,$ref,$forward_A,$forward_T,$forward_C,$forward_G,$reverse_A,$reverse_T,$reverse_C,$reverse_G,$stat,$heteroplasmy,$lower,$upper,$major_allele,$minor_allele,$atcg{$major_allele},$atcg{$minor_allele}];

        #calculate strand bias if $sb>0
        if($sb>0){
            my($forward_primary_allele_count,$forward_non_primary_allele_count,
               $reverse_primary_allele_count,$reverse_non_primary_allele_count)=(0,0,0,0);
            my %tmp;
            $tmp{'A'}->{'forward'}=$forward_A;
            $tmp{'A'}->{'reverse'}=$reverse_A;
            $tmp{'T'}->{'forward'}=$forward_T;
            $tmp{'T'}->{'reverse'}=$reverse_T;
            $tmp{'C'}->{'forward'}=$forward_C;
            $tmp{'C'}->{'reverse'}=$reverse_C;
            $tmp{'G'}->{'forward'}=$forward_G;
            $tmp{'G'}->{'reverse'}=$reverse_G;
            
            #Now whether isall=0 or 1, strand bias is calculated based on major allele and minor allele
            ($forward_primary_allele_count,$forward_non_primary_allele_count,
             $reverse_primary_allele_count,$reverse_non_primary_allele_count)=
            ($tmp{$major_allele}->{'forward'},$tmp{$minor_allele}->{'forward'},
             $tmp{$major_allele}->{'reverse'},$tmp{$minor_allele}->{'reverse'});
             
             $sb{$loc}=_sb($forward_primary_allele_count,$forward_non_primary_allele_count,
             $reverse_primary_allele_count,$reverse_non_primary_allele_count);
        }
    }
    
    #Re-loop the @result and write out
    if($sb>0){
        my @sbvalue = sort {$b<=>$a} values %sb;
        my $index = int (scalar(@sbvalue)*$sb/100)-1;
        $index=0 if ($index<0);
        $index=$#sbvalue if ($index>$#sbvalue);
        my $cutoff = $sbvalue[$index];
        foreach my $arrayref(@result){
            if($sb{$arrayref->[1]}>=$cutoff &&
                $arrayref->[11]==2){
                    $arrayref->[11]=3
            }
        }
    }
    
   #write out
   
    
}

#Calculate 95% confidence interval
#Parameters:
#$p  rate range from 0-1
#$n  total number of allele
sub _ci{
    my ($p,$n) = @_;
    my $tmp = 1.96*sqrt(($p*(1-$p))/$n);
    my $lower=0;
    my $upper=1;
    $lower=$p-$tmp if(($p-$tmp)>0);
    $upper=$p+$tmp if (($p+$tmp)<1);
    return($lower,$upper);
}


#Give a bam/sam file($inbam) and a mithochondrial bed file($mbed),
#fetch the alignment reads in mithochondrial genome
#and output in bam format
#Parameters:
#$inbam  input bam or sam format file
#$isbam  0/1, 0 indicates it's in sam format while 1 indicates a bam format
#$mbed   bed file of mithochondrial genome
#$mmq    minimum map quality
sub _get_mitochondrial_bam{
    my ($inbam,$isbam,$mbed,$mmq,$outbam) = @_;
    my $comm="";
    if($isbam){
        $comm = "$samtools view -b -q $mmq -L $mbed $inbam > $outbam";
    }else{
        $comm = "$samtools view -bS -q $mmq -L $mbed $inbam > $outbam";
    }
    system($comm);
}


#Given two basecall files and then determine the somatic mutations
#The first input basecall file is tumor sample while the second input basecall is normal sample
#Parameters:
#$tumorbase  tumor basecall input file
#$normalbase normal basecall input file
#$sp  (-sp) cutoff used in determine somatic mutation, means % of alternative allele observed in tumor
#$sa  (-sa) cutoff used in determine somatic mutation, means number of alternative allele observed in tumor
#$isall   (-A)  0/1, 1 denotes that the total read count is the total allele count of all allele observed,
#         while 0 indicates the total read count is the sum of major and minor allele counts, default=0
sub _determine_somatic{
    my ($tumorbase,$normalbase,$sp,$sa,$isall) = @_;
    my %tumor = _read_basecall($tumorbase);
    my %normal = _read_basecall($normalbase);
    
    my %result;
    
    #Get the overlapped locations and then loop them to determine whether it is a somatic mutation
    my @commloc = grep {exists($tumor{$_})} sort {$a<=>$b} keys %normal; #numeric sort
    foreach my $loc (@commloc){
        #Only if normal is not heteroplasmy
        my $normalGeno="";
        my $normalDP="";
        my $tumorGeno="";
        my $tumorDP="";
        my $issomatic=0;
        
        if($normal{$loc}->{'minor_allele_count'}==0){
            if($tumor{$loc}->{'minor_allele_count'}==0){
                #If at this site, both normal and tumor are not heteroplasmy, then determine whether their
                #major_allele is the same
                if($normal{$loc}->{'major_allele'} ne $tumor{$loc}->{'major_allele'}){
                    #a somatic mutation
                    $issomatic =1 ;
                }
            }else{
                if($normal{$loc}->{'major_allele'} ne $tumor{$loc}->{'major_allele'}){
                   #No matter the minor_allele_count meet the cutoff of not, it is already a somatic mutation
                   $issomatic=1;
                }else{
                    #Only the minor_allele_count meet the cutoff and it is a somatic mutation
                    my $ratio=0;
                    if($isall){
                        $ratio = $tumor{$loc}->{'minor_allele_count'}/(sum @{$tumor{$loc}->{'alleles_count'}});
                    }else{
                        $ratio = $tumor{$loc}->{'minor_allele_count'}/($tumor{$loc}->{'major_allele_count'}+$tumor{$loc}->{'minor_allele_count'});
                    }
                    
                    if($ratio>$sp/100 && $tumor{$loc}->{'minor_allele_count'}>$sa){
                        $issomatic=1;
                    }
                }
            }
        }
        
        #Assign this loc to %result
        if($issomatic){
            $normalGeno=$normal{$loc}->{'major_allele'};
            $normalDP=$normal{$loc}->{'major_allele_count'};
            
            my @tmp = grep {$_ != 0} @{$tumor{$loc}->{'alleles_count'}};
            $tumorDP=join "|", @tmp;
            $tumorGeno=join "|",@{$tumor{$loc}->{'alleles'}}[0..$#tmp];
            $result{$loc}->{'ref'}=$normal{$loc}->{'reference_allele'};
            $result{$loc}->{'normalGeno'}=$normalGeno;
            $result{$loc}->{'normalDP'}=$normalDP;
            $result{$loc}->{'tumorGeno'}=$tumorGeno;
            $result{$loc}->{'tumorDP'} = $tumorDP;
        }
    }
    return %result;
}


#Given a basecall format file, parse it and return a hash table
#Parameters:
#$in  input basecall file
sub _read_basecall{
    my($in) = @_;
    my %return;
    open(IN,$in) or die $!;
    <IN>;
    while(<IN>){
        s/\r|\n//g;
        my($chr,$loc,$ref,$forward_A,$forward_T,$forward_C,$forward_G,$reverse_A,$reverse_T,$reverse_C,$reverse_G) = split "\t";
        # Get the major allele and minor allele.
        my %atcg;
        $atcg{A}=$forward_A+$reverse_A;
        $atcg{T}=$forward_T+$reverse_T;
        $atcg{C}=$forward_C+$reverse_C;
        $atcg{G}=$forward_G+$reverse_G;
        my @atcg = sort {$atcg{$b}<=>$atcg{$a}} keys %atcg;
        my $major_allele=$atcg[0];
        my $minor_allele=$atcg[1];
        
        $return{$loc}->{'major_allele'}=$atcg[0];
        $return{$loc}->{'major_allele_count'}=$atcg{$atcg[0]};
        $return{$loc}->{'minor_allele'}=$atcg[1];
        $return{$loc}->{'minor_allele_count'}=$atcg{$atcg[1]};
        $return{$loc}->{'reference_allele'}=$ref;
        $return{$loc}->{'reference_allele_count'}=$atcg{$ref};
        $return{$loc}->{'alleles'}=[@atcg];
        $return{$loc}->{'alleles_count'}=[@atcg{@atcg}]; #Another ways to fetch value from a hash
        $return{$loc}->{'basecallline'}=[$chr,$loc,$ref,$forward_A,$forward_T,$forward_C,$forward_G,$reverse_A,$reverse_T,$reverse_C,$reverse_G];
    }
    close IN;
    return %return;
}



#Given a bam/sam file and run samtools mpileup to pile up the file
#Parameters:
#$inbam  input bam or sam format file
#$isbam  0/1, 0 indicates it's in sam format while 1 indicates a bam format
#$mbq    minimum base quality
#$refseq reference sequence to determine the reference allele
#$outpileup  pileup output file
sub _pileup{
    my($inbam,$isbam,$mbq,$refseq,$outpileup) = @_;
    my $comm="";
    if($isbam){
        $comm="$samtools -Q $mbq -f $refseq $inbam > $outpileup";
    }else{
        $comm="$samtools view -Su $inbam | $samtools -Q $mbq -f $refseq - > $outpileup";
    }
    system($comm);
}

#Given a pileup file, parse the pileup file into basecall format and write it out
#Parameters:
#$inpileup  one-sample pileup format file
#$mbq       minimum base quality
#$outbase   output allele count file parsed from pileup file

sub _parse_pileup{
    my($inpileup,$mbq,$outbase) = @_;
    open(IN,$inpileup) or die $!;
    open(OUT,">$outbase") or die $!;
    #Uppeer case is from forward strand, while lower case is from reverse strand
    print OUT "Chr\t"."Loc\t"."Ref\t"."A\t"."T\t"."C\t"."G\t"."a\t"."t\t"."c\t"."g\n";
    while(<IN>){
        s/\r|\n//g;
        my($chr,$loc,$ref,$dp,$bases,$bq) = split /\s+/;
        $ref=uc($ref);
        #do some modificaton on $bases to remove additional characters
        #1,remove the ^. pattern (marks the start of a read segment and the ASCII of the character following `^' minus 33 gives the mapping quality)
        $bases=~s/\^.//g;
        #2,remove the $ pattern (marks the end of a read segment)
        $bases=~s/\$//g;
        #3,remove -[0-9]+[ACGTNacgtn]+ pattern (denotes a deletion of one or more bases)
        my %hash=();
        while($bases=~/-(\d+)/g){
            $hash{$1}=1;
        }
        foreach my $k (keys %hash){
            $bases=~s/-$k[ACGTNacgtn]{$k}//g;
        }
        #4,remove +[0-9]+[ACGTNacgtn]+ pattern (denotes a insertion of one or more bases)
        %hash=();
        while($bases=~/\+(\d+)/g){
            $hash{$1}=1;
        }
        foreach my $k (keys %hash){
            $bases=~s/\+$k[ACGTNacgtn]{$k}//g;
        }
        #Now @base and @bq have the same length (Note that: the < or > in the $bases denote a gap)
        my @base=split (//,$bases);
        my @bq=split(//,$bq);
        my $A=0;
        my $T=0;
        my $C=0;
        my $G=0;
        my $a=0;
        my $t=0;
        my $c=0;
        my $g=0;
        #start the loop
        for(my $i=0;$i<@base;$i++){
            my $ch=$base[$i];
            my $score=ord($bq[$i])-32; #Need to be more robust later
            if($score>=$mbq){
                if($ch eq "A"){
                    $A++;
                }elsif($ch eq "T"){
                    $T++;
                }elsif($ch eq "C"){
                    $C++;
                }elsif($ch eq "G"){
                    $G++;
                }elsif($ch eq "a"){
                    $a++;
                }elsif($ch eq "t"){
                    $t++;
                }elsif($ch eq "c"){
                    $c++;
                }elsif($ch eq "g"){
                    $g++;
                }elsif($ch eq "."){
                    if($ref eq "A"){
                        $A++;
                    }elsif($ref eq "T"){
                        $T++;
                    }elsif($ref eq "C"){
                        $C++;
                    }elsif($ref eq "G"){
                        $G++;
                    }
                }elsif($ch eq ","){
                    if($ref eq "A"){
                        $a++;
                    }elsif($ref eq "T"){
                        $t++
                    }elsif($ref eq "C"){
                        $c++;
                    }elsif($ref eq "G"){
                        $g++;
                    }
                }
            }#end the condition  $score>=$mbq
        }#end the loop
        
        if($A+$T+$C+$G+$a+$t+$c+$g>0){
            print OUT "$chr\t$loc"."\t".$ref."\t".$A."\t".$T."\t".$C."\t".$G."\t".$a."\t".$t."\t".$c."\t".$g."\n";
        }
    }
    close IN;
    close OUT;
}

#Calculate Strandbias
#Parameters:
#$a  forward primary allele count
#$b  forward non-primary allele count
#$c  reverse primary allele count
#$d  reverse non-primary allele count
#Fomular:
#   _b_ _ _d_
#   a+b   c+d
#~~~~~~~~~~~~~~~~~
#   ___b+d____
#    a+b+c+d
#
sub _sb{
     my ($a,$b,$c,$d)=@_;
    if($a+$b==0 || $c+$d==0 || $b+$d==0){
        return 0;
    }else{
        return abs($b/($a+$b)-$d/($c+$d))/(($b+$d)/($a+$b+$c+$d));
    }
}


#Check whether samtools exist in your $PATH, If not, the program will exit
sub _check_samtools{
    my $r=`which samtools`;
    if($r){
        chomp($r);
        info("Checking: samtools in $r");
    }else{
       info("Error: samtools not installed yet");
       exit(1);
    }
}


#Customed print with time included, if providing $flag then there will be no "\n"
sub info{
    my ($s,$flag)=@_;
    print "[",scalar(localtime),"] $s";
    print "\n" unless($flag);
}

