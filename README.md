#Table of Content
================
* [Overview](#overview)
* [Usage](#usage)
* [Change log & Download](#change)
   * [Release version 1.3 on May 18, 2015](#v1.3)
   * [Release version 1.2 on Jan 15, 2014](#v1.2)
   * [Release version 1.1 on Feb 15, 2013](#v1.1)
   * [Release version 1.0 on Dec 14, 2012](#v1.0)
* [Statistical framework for heteroplasmy detection](#statistics)
   * [Fisher's exact test](#fisher)
   * [Empirical Bayesian for Binomial proportion with conjugate Beta prior](#bayes)
* [Prerequisites](#Prerequisites)
  * [Step1: Intall perl packages required by circos](#step1)
  * [Step2: Intall perl packages required by MitoSeek](#step2)
  * [Step3: Build samtools](#step3)
* [Others](#Others)
  * [Configure your perl environment](#perlsetup)
  * [Mitochondrial genome information (hg19/rCRS)](#mitogenome)
     * [Mitochondrial genome reference](#mitoreference)
     * [Mitochondrial genome annotation](#mitoanno)
     * [Known pathogenic mutations](#pathogenic)
  * [Whole genome exon coordinates](#exon)
    

<a name="overview"/>

# Overview
MitoSeek is an open-source software tool to reliably and easily extract mitochondrial genome information from exome sequencing data. MitoSeek evaluates **mitochondrial genome alignment quality, estimates relative mitochondrial copy numbers, and detects heteroplasmy, somatic mutation, and structural variance of the mitochondrial genome**.

<a name="usage"/>

## Usage
Example code for running MitoSeek with given toy dataset could be

```shell
perl mitoSeek.pl -i Examples/brca_tumor.bam -j Examples/brca_normal.bam -t 4 -sb 0 -hp 1 -d 5 -str 4 -sp 1 -sa 0
```

This example report could be accessed at [Here]( http://htmlpreview.github.com/?https://github.com/riverlee/MitoSeek/blob/release/brca_tumor/mitoSeek.html)


Details for each parameters are here
```shell
Usage: perl mitoSeek.pl -i inbam 
-i [bam]                Input bam file
-j [bam]                Input bam file2, if this file is provided, it will conduct somatic mutation mining, and it will be 
                        taken as normal tissue.
-t [input type]         Type of the bam files, the possible choices are 1=exome, 2=whole genome, 3= RNAseq, 4 = mitochondria only,default = 1
-d [int]                The minimum recommended depth requirement for detecting heteroplasmy. Lower depth will severely damage the 
                        confidence of heteroplasmy calling, default=50
-ch                     Produce circos plot input files and circos plot figure for heteroplasmic mutation,
                        (-noch to turn off and -ch to turn on), default = on
-hp [int]               Heteroplasmy threshold using [int] percent alternative allele observed, default = 5
-ha [int]               Heteroplasmy threshold using [int] allele observed, default = 0
-alpha [float]          Shape1 parameter of Beta prior distribution, default is 3.87 which is estimated from 600 BRCA samples
-beta  [float]          Shape2 parameter of Beta prior distribution, default is 174.28, which is estimated from 600 BRCA samples
-A                      If - A is used, the total read count is the total allele count of all allele observed. 
                        Otherwise, the total read count is the sum of major and minor allele counts. Default = off
-mmq [int]              Minimum map quality, default =20
-mbq [int]              Minimum base quality, default =20
-sb [int]               Remove all sites with strand bias score in the top [int] %, default = 10 
-cn                     Estimate relative copy number of input bam(s), does not work with mitochondria targeted sequencing bam files,
                        (-noch to turn off and -ch to turn on) default = off.
-sp [int]               Somatic mutation detection threshold,int = percent of alternative allele observed in tumor, default int=5
-sa [int]               Somatic mutation detection threshold,int = number of alternative allele observed in tumor, default int=3
-cs                     Produce circos plot input files and circos plot figure for somatic mutation, 
                        (-nocs to turn off and -cs to turn on), default = off
-r [ref]                The reference used in the bam file, the possible choices are hg19 and rCRS, default=hg19
-R [ref]                The reference used in the output files, the possible choices are hg19 and rCRS, default=hg19
-str [int]              Structure variants cutoff for those discordant mapping mates, 
                        int = number of spanning reads supporting this structure variants, default = 2
-strf [int]             Structure variants cutoff for those large deletions,
                        int = template size in bp, default=500
-QC                     Produce QC result, (--noQC to turn off and -QC to turn on), default=on
-samtools[samtools]     Tell where is the samtools program, default is your mitoseek directory/Resources/samtools/samtools
-bwa [bwa]              Tell where is the bwa program, default value is 'bwa' which is in your $PATH
-bwaindex [bwaindex]    Tell where is the bwa index of non-mitochondrial human genome, no default value
-advance                Will get mitochondrial reads in an advanced way, generally followed by 1) Initially extract mitochrodrial reads from 
                        a bam file, then 2) remove those could be remapped to non-mitochondrial human genome by bwa. Advanced extraction needs 
                        -bwaindex option. Default extraction without removing step.

```

<a name="change"/>

## Change log & Download

<a name="v1.3">
  
### Release version 1.3 on May 18, 2015
Fix the bug when parsing pileup file with depth=0


* **[Download Zip](https://github.com/riverlee/MitoSeek/archive/v1.3.zip)**
* **[Browse Code ](https://github.com/riverlee/MitoSeek/tree/v1.3)**


<a name="v1.2">
  
### Release version 1.2 on Jan 15, 2014
Version 1.2 has been improved by reading options from a configure file instead of reading from command line

* **[Download Zip](https://github.com/riverlee/MitoSeek/archive/v1.2.zip)**
* **[Browse Code](https://github.com/riverlee/MitoSeek/tree/v1.2)**

Changes are here:
  * Reading options from a configure instead of reading from commmand line. The new program is called **mitoSeek_new.pl**, an example of configure file is called **'para.txt'**. Meanwhile, the original program which reads options from command line is also kept, called **mitoSeek.pl**.

<a name="v1.1">
  
### Release version 1.1 on Feb 15, 2013
Version 1.1 has been improved according to reviewers' advise.

* **[Download Zip](https://github.com/riverlee/MitoSeek/archive/v1.1.zip)**
* **[Browse Code ](https://github.com/riverlee/MitoSeek/tree/v1.1)**

Changes are here:
  * Add option **-advance** to improve mitochondrial reads extraction, which is done by 1) Initially extract mitochrodrial reads from a bam file, then 2) remove those could be remapped to non-mitochondrial human genome by bwa. The version **v1.0** only implementes the step 1
  * Add Statistical framework for heteroplasmy detection, which are [fisher test](#fisher) and [empirical bayesian](#bayes). 
  * Add option **-samtools** for people who would like to use their specified samtools
  * Add options **-bwa** and **--bwaindex** which are needed if use **-advance** option
  * In the heteroplasmy output, it will contain columns of **fisher.pvalue, fisher.adjust.pvalue,fisher.phred.score,empirical.probability,** and **empirical.phred.score**
  * Thanks to Peter's code for beta calculation in perl at <http://home.online.no/~pjacklam/perl/modules/>


<a name="v1.0"/>

### Release version 1.0 on Dec 14, 2012
Initial version for the paper

* **[Download Zip](https://github.com/riverlee/MitoSeek/archive/964cd2e61735a60283f0280020cadbb53be3617e.zip)**
* **[Browse Code ](https://github.com/riverlee/MitoSeek/tree/964cd2e61735a60283f0280020cadbb53be3617e)**

<a name="statistics"/>

## Statistical framework for heteroplasmy detection
We have implemented statistical framework in addition to the empirical filters. We added a function to perform a one-tail Fisher’s exact test to determine if the rate of heteroplasmy is significantly greater than user defined cutoff (-hp). Moreover, MitoSeek also provides Phred quality scores based on the p-value.


<a name="fisher"/>
### Fisher's exact test

```
           major   minor    
observed    n11     n12 | n1p 
expected    n21     n22 | n2p   
           -----------------
            np1     np2   npp
where n11 and n12 are observed number of major and minor alleles,
n21 = (n11+n12)*(1-hp/100) in which hp is defined by -hp
n22 = (n11+n12)*hp/100  in which hp is defined by -hp
```

The phred score of heteroplasmy for Fisher's exact test calucated as

```
phred.score.fisher = - log10 * log10(fisher.pvalue)
```

<a name="bayes"/>
### Empirical Bayesian for Binomial proportion with conjugate Beta prior

```
                                      hp          
                                    _ ---         
                                   /  100         
likelihood.of.heteroplasmy =  1 -  |      f(x) dx 
                                  _/  0           


```

where

```
                        1             a  + n12 - 1         b  + n11 - 1  
f(x) =  -------------------------- * x            * (1 - x)
        beta(n12 +  a ,n11 +  b )                             
                                                                     

```

**_a_** and **_b_** are share parameters of Beta prior distribution, which are  estimated from 600 BRCA samples.


The phred score of heteroplasmy for Empirical Bayesian approach is calculated as

```
phred.score.empirical = - log10 * log10(1 - likelihood.of.heteroplasmy)


```




<a name="Prerequisites"/>

# Prerequisites

MitoSeek runs on 32-bit or 64-bit GNU/Linux and request perl packages like **GD::Graph::boxplot**, etc. Here are the steps you need to do to install packages needed for the **MitoSeek** program.


<a name="step1"/>

## Step1: Install perl packages required by [circos](http://circos.ca/)
MitoSeek utilizes [circos](http://circos.ca/) to plot heteroplasmy and somatic mutation, thus, perl packages required by [circos](http://circos.ca/) needed to been installed first.

```
#1) go to the circos package folder which is included as part of the MitoSeek 
cd Resources/circos-0.56/bin

#2) Check whether all the packages needed by circos plot are intalled on your PC,
#if this does not work, try 'chmod +x test.modules'
./test.modules
```
If all the packages are installed on your PC, it will look like this.
```
ok   Carp
ok   Config::General
ok   Cwd
ok   Data::Dumper
ok   Digest::MD5
ok   File::Basename
ok   File::Spec::Functions
ok   File::Temp
ok   FindBin
ok   GD
ok   GD::Image
ok   GD::Polyline
ok   Getopt::Long
ok   IO::File
ok   List::MoreUtils
ok   List::Util
ok   Math::Bezier
ok   Math::BigFloat
ok   Math::Round
ok   Math::VecStat
ok   Memoize
ok   Params::Validate
ok   Pod::Usage
ok   POSIX
ok   Readonly
ok   Regexp::Common
ok   Set::IntSpan
ok   Storable
ok   Sys::Hostname
ok   Text::Format
ok   Time::HiRes
```
Otherwise, it may look like this.
```
ok   Carp
fail Config::General is not usable (it or a sub-module is missing)
ok   Cwd
ok   Data::Dumper
ok   Digest::MD5
ok   File::Basename
ok   File::Spec::Functions
ok   File::Temp
ok   FindBin
fail GD is not usable (it or a sub-module is missing)
fail GD::Image is not usable (it or a sub-module is missing)
fail GD::Polyline is not usable (it or a sub-module is missing)
ok   Getopt::Long
ok   IO::File
fail List::MoreUtils is not usable (it or a sub-module is missing)
ok   List::Util
fail Math::Bezier is not usable (it or a sub-module is missing)
ok   Math::BigFloat
fail Math::Round is not usable (it or a sub-module is missing)
fail Math::VecStat is not usable (it or a sub-module is missing)
ok   Memoize
ok   Params::Validate
ok   Pod::Usage
ok   POSIX
fail Readonly is not usable (it or a sub-module is missing)
fail Regexp::Common is not usable (it or a sub-module is missing)
fail Set::IntSpan is not usable (it or a sub-module is missing)
ok   Storable
ok   Sys::Hostname
fail Text::Format is not usable (it or a sub-module is missing)
fail Time::HiRes is not usable (it or a sub-module is missing)
```

If this happens, try to install the missing packages by **cpan** (If you you don't have root previlege, please look at [here](#perlsetup) to set up your own perl environment).
```shell
#run in root if the packages will be installed in the system path
#like /usr/local/lib64/perl5
./test.modules |grep fail|cut -f2 -d" "|xargs -I {} cpan {}
``` 
<a name="step2"/>

## Step2: Install perl packages required by MitoSeek
In addition to the perl packages required by [circos](http://circos.ca/), there are several other packages needed to be installed on your PC.

```shell
#go the the folder where your MitoSeek is
#And test whether all the required modules have been installed.
./test.modules
```
The successful output would look like this
```shell
ok   Convert
ok   Cwd
ok   File::Basename
ok   File::Path
ok   FindBin
ok   GD
ok   GD::Graph::bars3d
ok   GD::Graph::boxplot
ok   GD::Graph::lines
ok   GD::Text::Wrap
ok   Getopt::Long
ok   List::Util
ok   Mitoanno
ok   Statistics::KernelEstimation
```
Otherwise, it may look like this
```shell
ok   Convert
ok   Cwd
ok   File::Basename
ok   File::Path
ok   FindBin
fail GD is not usable (it or a sub-module is missing)
fail GD::Graph::bars3d is not usable (it or a sub-module is missing)
fail GD::Graph::boxplot is not usable (it or a sub-module is missing)
fail GD::Graph::lines is not usable (it or a sub-module is missing)
fail GD::Text::Wrap is not usable (it or a sub-module is missing)
ok   Getopt::Long
ok   List::Util
ok   Mitoanno
fail Statistics::KernelEstimation is not usable (it or a sub-module is missing)
```
To install the missing packages is the same way as we did in [step1[)#step1) for those missing packages requried by [circos](http://circos.ca/).
```
#run in root if the packages will be installed in the system path
#like /usr/local/lib64/perl5
./test.modules |grep fail|cut -f2 -d" "|xargs -I {} cpan {}
``` 
<a name="step3"/>

## Step3: Build samtools
We include [samtools](http://samtools.sourceforge.net/) as part of MitoSeek, however, you need to build it before you use MitoSeek.

```
#go the samtools folder
cd Resources/samtools
tar jxvf samtools-0.1.18.tar.bz2
cd samtools-0.1.18
#build samtools and move the samtools program into its parental folder
make
mv samtools ../
#remove the samtools source code if you like
cd ../
rm -rf samtools-0.1.18
```

Or you can change the **mitoSeek.pl** script and tell the program where is the samtools on your PC. Here is the line you need to modify.
```perl
my $samtools = "$FindBin::Bin/Resources/samtools/samtools";  #Where is the samtools file
```


<a name="others"/>

Others
-----------------------------------

<a name="perlsetup"/>

## Configure your perl environment
Ususally you don't have root privilege and it will prevent you when you try to install perl packages by default. Here is a brief way to install perl packages without root previlege, details could be found at http://www.perl.com/pub/2002/04/10/mod_perl.html or  searching by google.

### Step 1, Setting your PERL5LIB environment variable
```
#If your own perl libary top folder is ~/perllib
#Then add the following lines in your ~/.bashrc
#because the subfolders could be different due to different versions of perl
#and linux, so to save time by adding subfolders in the PERL5LIB, we use 'for'
#to add all the subfolders. 
PERL5LIB=~/perllib
if [ -d ~/perllib ]; then
  for i in `find ~/perllib -type d`
    do
      PERL5LIB=$i:$PERL5LIB
    done
fi
export $PERL5LIB
```

### Step 2, Configure your cpan
```shell
#type cpan into cpan environment
cpan
```
Then
```shell
#In cpan environment
o conf makepl_arg PREFIX="~/perllib"
o conf commit

#exit the cpan environment
exit
```

<a name="mitogenome"/>

## Mitochondrial genome information (hg19/rCRS)
Mitochondrial information we used in **MitoSeek** includes
* [Mitochondrial genome reference](#mitoreference)
* [Mitochondrial genome annotation](#mitoanno)
* [Known pathogenic mutations](#pathogenic)

<a name="mitoreference"/>

### Mitochondrial genome reference
**result folder:** yourmitoseek/Resources/genome

**Files:**
   + hg19.fasta
   + rCRS.fasta

**Links**

```shell
#1.1) rCRS assembly (rCRS.fasta)
http://www.ncbi.nlm.nih.gov/nuccore/251831106

#1.2) hg19 assembly (hg19.fasta)
http://www.ncbi.nlm.nih.gov/nuccore/NC_001807.4?report=genbank
```

<a name="mitoanno"/>
### Mitochondrial genome annotation
Annotation is in **genbank** format, the class/package **Mitoanno** (Mitoanno.pm) will take the .gb file and annotate the given position on mitochondria.
   
**result folder:** yourmitoseek/Resources/

**Files:**
   + hg19_annotation.gb
   + rCRS_annotation.gb

**Links**

```shell
#1.1) rCRS assembly (download as genbank format)
http://www.ncbi.nlm.nih.gov/nuccore/251831106

#1.2) hg19 assembly (downlaod as genbank format)
http://www.ncbi.nlm.nih.gov/nuccore/NC_001807.4?report=genbank
```

<a name="pathogenic"/>

### Known pathogenic mutations
Known pathogenic mutations on mitochondria comes from [OMIM](http://omim.org/) and it will be used by class/package **Mitoanno** (Mitoanno.pm)

**File Name:** yourmitoseek/Resources/OMIMpathogenic.txt


<a name="exon"/>

## Whole genome exon coordinates
Whole genome exon coordinates are used to estimate the **relative mitochondrial copy numbers** when the input is exon/RNA-Seq sequencing. 

Here is the scripts we download and prepare the exon bed file

```shell
#Wkdir: yourmitseek/Resources/genome
#refGene annotation file (refGene.txt)
#description of refGene is here at 
#https://cgwb.nci.nih.gov/cgi-bin/hgTables?hgsid=79902&hgta_doSchemaDb=hg18&hgta_doSchemaTable=refGene
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
gunzip refGene.txt.gz
    
#get exon region on chromosome 1-22, and X,Y (with prefix chr or without)
perl -lane 'if($F[2]=~/chr\d+$|chrX|chrY/) {@start=split ",",$F[9]; @end=split ",",$F[10];  for($i=0;$i<$F[8];$i++) { print join "\t",$F[2],$start[$i],$end[$i],$F[12]."|".$F[1]}}' refGene.txt >refGeneExon_withChr.bed
   
perl -lane 'if($F[2]=~/chr\d+$|chrX|chrY/) {$F[2]=~s/chr//g;@start=split ",",$F[9]; @end=split ",",$F[10];  for($i=0;$i<$F[8];$i++) { print join "\t",$F[2],$start[$i],$end[$i],$F[12]."|".$F[1]}}' refGene.txt >refGeneExon_withoutChr.bed
   
#How to get the total bases from a bed file (getTotalBasesFromBed.R)
./getTotalBasesFromBed.R genome_withChr.bed
  
```






