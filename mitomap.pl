#!/usr/bin/env perl
###################################
# Author: Jiang (River) Li
# Email:  riverlee2008@gmail.com
# Date:   Wed Feb 13 10:00:36 2013
###################################
use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use FindBin;

# After you extract mitochondria from exome sequencing, you will have a mitochondria only Bam file.  Do a bwa alignment on non mitochondria human reference.  If anything mapped to human reference remove those reads.

# Global variables
my $samtools    = "$FindBin::Bin/Resources/samtools/samtools";   #Where is the samtools file
my $bwa         = "$FindBin::Bin/Resources/bwa/bwa";

# Default parameters
my $mmq         = 20;   #minimal mapping quanlity
my $inbam;
my $mitobed;
my $bwaindex;
my $outbam;
my $help;

my $paired=0;  
# Assign values to variables
unless(
    GetOptions(
        "i=s"   =>\$inbam,
        "b=s"   =>\$mitobed,
        "r=s"   =>\$bwaindex,
        "o=s"   =>\$outbam,
        "mmq=i" =>\$mmq,
        "samtools=s"=>\$samtools,
        "bwa=s" =>\$bwa,
        "h"     =>\$help)
){
    print $!,"\n";
    usage(1);
}

#
info("Checking parameters ...",1);
check();
$paired = determine_paired($inbam);

my %tempfiles=("mitobam1"=>"mito_step1.bam",          #extracted mapped reads in mito
               "single_sai"=>"mito_single.sai",       
               "pair_sai1"=>"mito_pair_1.sai",
               "pair_sai2"=>"mito_pair_2.sai",
               "mitosam2" =>"mito_step2.sam",
               "mitobam2"=>"mito_step2.bam");        #remap mitoreads to human genome
# Step 1, Get mapped read in mito
info("Extracting mapped reads on mitochondria ...",1);
get_mito_bam($inbam,$mitobed,$mmq,$tempfiles{'mitobam1'});

# Step 2, Remap mitoreads to non-mitochondrial human reference
info("Remapping those reads to non-mitochondria human genome ...",1);
bwa_map_bam($tempfiles{'mitobam1'},$bwaindex,$paired,$tempfiles{'single_sai'},$tempfiles{'pair_sai1'},
            $tempfiles{'pair_sai2'},$tempfiles{'mitosam2'},$tempfiles{'mitobam2'});

# Step 3, Remove those could be remapped to non-mitochondria human genome
info("Removing those remapped reads and get the finial mitochondrial read ...",1);
get_final_mito_bam($tempfiles{'mitobam1'},$tempfiles{'mitobam2'},$outbam);


# remvoe the tempfiles
END{
    foreach my $k (%tempfiles){
        if(-e $tempfiles{$k}){
            unlink $tempfiles{$k};
        }
    }
}


sub get_final_mito_bam{
    my($initalbam,$remappedbam,$out) = @_;
    # Read mapped read name in remappedbam file and store it in %hash
    my %hash;
    open(IN,"$samtools view $remappedbam|") or die $!;
    my $flag_read_unmapped=0x0004;
    while(<IN>){
        s/\r|\n//g;
        my (
            $qname, $flag,  $rname, $pos, $mapq, $cigar,
            $rnext, $pnext, $tlen,  $seq, $qual, $others
        ) = split "\t";
        next if($flag & $flag_read_unmapped); #read not mapped
        $hash{$qname}++;
    }
    close IN;
    
    # Read initialbam and remove those remapped
    open(IN,"$samtools view -h $initalbam|") or die $!;
    my $tmpsam="${initalbam}.tmp.sam";
    open(TMP,">$tmpsam") or die $!;
    while(<IN>){
        if(/^@/){ #print header
            print TMP $_;
            next;
        }
        my (
            $qname, $flag,  $rname, $pos, $mapq, $cigar,
            $rnext, $pnext, $tlen,  $seq, $qual, $others
        ) = split "\t"; 

        # output read if not being remappd
        if(!exists($hash{$qname})){
            print TMP $_;
        }
    }
    close TMP;
    close IN;

    #Conver to bam
    my $comm = "$samtools view -bS $tmpsam > $out";
    run($comm);
    unlink $tmpsam;
}

sub get_final_mito_bam2{
    my($initalbam,$remappedbam,$out) = @_;
    my $comm="$intersectBed -abam $initalbam -bbam $remappedbam -v > $out";
    run($comm);
}

sub bwa_map_bam{
    my($in,$index,$paired,$singlesai,$pairsai1,$pairsai2,$outsam,$out) = @_;
    if($paired){
        my $comm="bwa aln -b1 $index $in > $pairsai1";
        run($comm);
        $comm="bwa aln -b2 $index $in >$pairsai2";
        run($comm);
        $comm="bwa sampe $index $pairsai1 $pairsai2 $in $in >$outsam";
        run($comm);
        $comm="$samtools view -bS $outsam>$out";
        run($comm);

    }else{
        my $comm="bwa aln -b $index $in>$singlesai";
        run($comm);
        $comm="bwa samese $index $singlesai $in >$outsam";
        run($comm);
        $comm="$samtools view -bS $outsam>$out";
        run($comm);
    }
}


# Give a bed file, extract reads in this region
sub get_mito_bam{
    my ($in,$mbed,$mmq,$out) = @_;
    my $comm="$samtools view -b -q $mmq -L $mbed $in>$out";
    run($comm);
}


sub run{
    my ($comm) = @_;
    if($comm !~/>/){  #suppress the output
        $comm.=">/dev/null";
    }
    system($comm) == 0 or die $!;
}

# Read the first 100 reads from a bam and determine whether single-end or pair-end
sub determine_paired{
    my ($in) = @_;
    open(IN,"$samtools view $in|") or die $!;

    my $count=0;
    my $flag_paired=0x0001;
    my $isparied=0;

    while(<IN>){
        s/\r|\n//g;
        my (
            $qname, $flag,  $rname, $pos, $mapq, $cigar,
            $rnext, $pnext, $tlen,  $seq, $qual, $others
        ) = split "\t";
        if ( $flag & $flag_paired ) {
            $isparied = 1;
            last;
        }
        if($count>100){
            last;   
        }                            
    }
    close IN;
    return ($isparied);
}

# Check variables
sub check{
    my $msg ="";

    if(!defined($inbam) || $inbam eq ""){
        $msg="Input bam (-i) hasn't been defined\n";
    }elsif(! -e $inbam){
        $msg.="Input bam (-i) '$inbam' doesn't exists\n";
    }

    if(!defined($mitobed) || $mitobed eq ""){
        $msg.="Mitochondrial bed file (-b) hasn't been defined\n";
    }elsif(! -e $mitobed){
        $msg.="Mitochondrial bed file (-b) '$mitobed' doesn't exists\n";
    }

    if(!defined($bwaindex) || $bwaindex eq ""){
        $msg.="bwa index (-r) hasn't been defined\n";
    }elsif(! -e $bwaindex){
        $msg.="bwa index (-r) '$bwaindex' doesn't exists\n";
    }

    if(!defined($outbam) || $outbam eq ""){
        $msg.="output bam (-o) hasn't been defined\n";
    }

    #check samtools
    my $r = `which $samtools`;
    if(!-e $samtools && !$r){
        $msg.="samtools (-samtools) '$samtools' doesn't exists\n"
    }

    #check bwa
    $r = `which $bwa`;
    if(!-e $bwa && !$r){
        $msg.="bwa (-bwa) '$bwa' doesn't exists\n"
    }

    if($msg){
        print $msg,"\n";
        usage(1);
    }
}
sub usage{
    my ($flag) = @_; # if $flag defined, will exit the program
    print <<USAGE;
Usage: perl mitoSeek.pl -i <inbam> -b <mito.bed> -r <bwa_index> -o <output.bam> -mmq [20] -samtools [samtools] -bwa [bwa]
-i <bam>                Input bam file
-b <bed>                Input bed file region of mitochondria
-r <index>              bwa index file of non-mitochondria human genome
-o <bam>                Output bam file only containning mapped reads in mitochondrial
-mmq <int>              Minimum mapping quality, default=20
-samtools [samtools]    Tell where is the samtools program, default is 
                        your mitoseek directory/Resources/samtools/samtools
-bwa [bwa]              Tell where is the bwa program, default is
                        your mitoseek directory/Resources/bwa/bwa
-h                      print out this page

USAGE
    if($flag){
        exit(1);
    }
}


sub info{
    my ($str,$flag) = @_;
    print "[",scalar(localtime),"] $str";
    print "\n" if ($flag);
}
