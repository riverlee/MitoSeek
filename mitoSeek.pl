#!/usr/bin/perl
#############################################
#Author: Jiang Li
#email:  riverlee2008@gmail.com
#Creat Time: Tue 23 Oct 2012 01:37:54 PM CDT
#Vanderbilt Center for Quantitative Sciences
#############################################
use strict;
use warnings;
#========Begin loading necessary packages===========#
use GD;
use GD::Text::Wrap;
use GD::Graph::boxplot;
use Statistics::KernelEstimation;
use GD::Graph::lines;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use File::Path;
use Cwd 'abs_path';
use FindBin;
use List::Util qw(sum min max);
use Text::NSP::Measures::2D::Fisher::left;  # For fisher test
use Statistics::Multtest qw(BH);
use Math::SpecFun::Beta qw(beta);
use lib "$FindBin::RealBin";
## load our own packages
use Convert;
use Mitoanno;
use Circoswrap;
#========End loading necessary packages================#

#========Begin defining global variables===============#
my $starttime                     = time();                                          #Assign this value in the BEGIN section
my $commandline                   = "perl $0 ".join (" ",@ARGV);                     #Output you input command in the report to reproduce you result
my %mitogenome                    = (                                               #Will be used in the pileup to determine the reference allele
                                        "hg19" => $FindBin::Bin . "/Resources/genome/hg19.fasta",
                                        "rCRS" => $FindBin::Bin . "/Resources/genome/rCRS.fasta");
my %exonbed                       = (                                                #Human exon regions in bed format
                                        "withchr"=>$FindBin::Bin."/Resources/genome/refGeneExon_withChr.bed",
                                        "withoutchr"=>$FindBin::Bin."/Resources/genome/refGeneExon_withoutChr.bed");
my %genomebed                     = (                                                #Human chromosome regions in bed format
                                        "withchr"=>$FindBin::Bin."/Resources/genome/genome_withChr.bed",
                                        "withoutchr"=>$FindBin::Bin."/Resources/genome/genome_withoutChr.bed");
my $totalgenomebases              =3095677412;                                       #Total human genome size, used to calculate to average depth
my $totalexonbases                =70757781;                                         #Total human exon size, used to calculate the average depth
my %acceptedgenomelength          = ( "hg19" => 16571, "rCRS" => 16569 );            #The total mitochondrial genome size of hg19 and rCRS
my $samtools                      = "$FindBin::Bin/Resources/samtools/samtools";                                      #Where is the samtools file
my $mitomap                     = "$FindBin::Bin/mitomap.pl";                                      #Where is the samtools file

my $isbam                         = 1;                                               #Default is 1, will auto determined from the $inbam1 file
my $ischr                         = 0;                                               #Default is 0, means the chromosomes are named without prefix chr, will auto determined from the $inbam1's header
# sam/bam flag (will be used in the function _mito_qc_stat)
my $flag_paired                   = 0x0001;
my $flag_properly_paired          = 0x0002;
my $flag_read_unmapped            = 0x0004;
my $flag_next_read_unmapped       = 0x0008;
my $flag_reverse_strand           = 0x0010;
my $flag_next_read_reverse_strand = 0x0020;
my $flag_first_fragment           = 0x0040;
my $flag_second_fragment          = 0x0080;

my $debug = 0;


# Other variables will be used as output file name
my $folder                             = undef;                              #determine from the $inbam1 file, and will create this folder, all the report files will be put in this folder
my $mitobam1                           = "mito1.bam";                        #Reads aligned to mitochondria from inbam1, use this as output file name
my $mitobam2                           = "mito2.bam";                        #Reads aligned to mitochondria from inbam2m use this as output file name
my $reference                          = "mito.fasta";                       #copy $mitogenome{$inref} to $reference; will rename the chromosome name '>' when necessary
my $mitostart                          = 0;                                  #Mitochondrial region start for the analysis
my $mitoend                            = $acceptedgenomelength{'hg19'};      #Mitochondrial region end for the analysis

my $perbasequality_figure1             = "per_base_quality1.png";            #
my $perbasequality_table1              = "per_base_quality_table1.txt";     #
my $mappingquality_figure1             = "mapping_quality1.png";             #
my $mappingquality_table1              = "mapping_quality_table1.txt";       #
my $depthdistribution_figure1          = "depth_distribution1.png";          #
my $depthdistribution_table1           = "depth_distribution_table1.txt";    #
my $templatelengthdistribution_figure1 = "template_length_distribution1.png";#
my $templatelengthdistribution_table1  = "template_length_distribution_table1.txt";
my $perbasequality_figure2             = "per_base_quality2.png";             #
my $perbasequality_table2              = "per_base_quality_table2.txt";      #
my $mappingquality_figure2             = "mapping_quality2.png";              #
my $mappingquality_table2              = "mapping_quality_table2.txt";        #
my $depthdistribution_figure2          = "depth_distribution2.png";           #
my $depthdistribution_table2           = "depth_distribution_table2.txt";     #
my $templatelengthdistribution_figure2 = "template_length_distribution2.png";
my $templatelengthdistribution_table2  = "template_length_distribution_table2.txt";
my $percentofbasepairscovered_table1   = "percent_of_base_pairs_covered_table1.txt";
my $percentofbasepairscovered_table2   = "percent_of_base_pairs_covered_table2.txt";

my $mitopileup1                        = "mito1.pileup";
my $mitopileup2                        = "mito2.pileup";
my $mitobasecall1                      = "mito1_basecall.txt";
my $mitobasecall2                      = "mito2_basecall.txt";
my $mitoheteroplasmy1                  = "mito1_heteroplasmy.txt";
my $mitoheteroplasmy2                  = "mito2_heteroplasmy.txt";

my $mitostructure1                     = "mito1_structure_discordant_mates.txt";
my $mitostructure2                     = "mito2_structure_discordant_mates.txt";
my $mitostructuredeletion1             = "mito1_structure_large_deletion.sam";
my $mitostructuredeletion2             = "mito2_structure_large_deletion.sam";

my $mitocnv1                           = "mito1_cnv.txt";       #Result of cnv of mito1
my $mitocnv2                           = "mito2_cnv.txt";
my $mitodepth1                         = "mito1_depth.txt";     #Store depth of mito1
my $mitodepth2                         = "mito2_depth.txt";
my $sampledepthi                       = "sample_i_depth.txt";  #store depth of imput sample1 
my $sampledephtj                       = "sample_j_depth.txt";
my $mitosomatic                        = "mito_somatic_mutation.txt";
my $mitoreport                         = "mitoSeek.html";

#circos related
my $mitocircosheteroplasmyfigure1      = "mito1_heteroplasmy_circos.png";
my $mitocircosheteroplasmyconfig1      = "mito1_heteroplasmy_circos.conf";
my $mitoheteroplasmytextoutput1        = "mito1_heteroplasmy_circos.text.txt";
my $mitoheteroplasmyscatteroutput1     = "mito1_heteroplasmy_circos.scatter.txt";
my $mitocircosheteroplasmyfigure2      = "mito2_heteroplasmy_circos.png";
my $mitocircosheteroplasmyconfig2      = "mito2_heteroplasmy_circos.conf";
my $mitoheteroplasmytextoutput2        = "mito2_heteroplasmy_circos.text.txt";
my $mitoheteroplasmyscatteroutput2     = "mito2_heteroplasmy_circos.scatter.txt";
my $mitocircossomaticfigure            = "mito_somatic_mutation_circos.png";
my $mitocircossomaticconfig            = "mito_somatic_mutation_circos.config";
my $mitosomatictextoutput              = "mito_somatic_mutation_circos.text.txt";

my $mitooffset1                        = 33;                            #Not used currently
my $mitooffset2                        = 33;
# There three variables will be re-assigned during the main
my $mitobases                          = $acceptedgenomelength{'hg19'};  #Read from the regionbed file,default is the total bases of hg19
my $totalbases                         = $totalexonbases;                #This the genome length, the exon length is 70757781
my $totalbed                           = $exonbed{'withchr'};            #Default is exon or RNA-Seq sequences

#Our own defined class, will be used in the _determine_heteroplasmy
my $convert = Convert->new();                                           #Handle h19 position and rCRS position convertion
my $mitoanno = Mitoanno->new();                                         #Annotating a give position on mitochondria
my $circos   = Circoswrap->new();
#========End defining global variables===================#

#========Begin defining other variables==================#
my $inbam1            = undef;      #-i
my $inbam2            = undef;      #-j ,if this is provided, will conduct somatic mutation mining, and this will be assumed as normal sample
my $type              = 1;          #1=exome, 2=whole genome, 3= RNAseq, 4 = mitochondria only
#my $savebam           = 1;          #-b
#my $saveallelecount   = 1;          #-a
my $producecircosplot = 1;          #-ch, produce circos plot for heteroplasmic mutation
my $hp                = 5;          #-hp, heteroplasmy threshold using [int] percent alternatie allele observed, default=5;
my $ha                = 0;          #-ha, heteroplasmy threshold using [int] allele observed, default=0;
my $depth             = 50;         #The minimum recommended depth requirement for detecting heteroplasmy is 50. Lower depth will severely damage the confidence of heteroplasmy calling
my $isall             = 0;          #If - A is used, the total read count is the total allele count of all allele observed. Otherwise, the total read count is the sum of major and minor allele counts. Default = off
my $mmq               = 20;         #minimum map quality, default=20
my $mbq               = 20;         #minimum base quality, default=20
my $sb                = 10;         #remove all sites with strand bias score in the top [int] %, default=10;
my $cn                = 0;          #Estimate relative copy number of input bam(s),does not work with mitochondria targeted sequencing bam files
my $sp                = 5;          #somatic mutation detection threshold, [int]% of alternative allele observed in tumor, default=5;
my $sa                = 3;          #somatic mutation detection trheshold, int number of alternative allele observed in tumor.
my $cs                = 1;          #Produce circos plot input files and circos plot figure for somatic mutations
my $regionbed         = undef;      #A bed file that contains the regions mitoSeek will perform analysis on
my $inref             = 'hg19';     #The reference used in the input bam files
my $outref            = 'hg19';     #The output files used in the reference rCRS;
my $qc                = 1;          #Produce QC result
my $str               = 2;          #structure variants cutoff, this cutoff applied to >$str mates supporting this cross different chromosome mapping
my $strflagmentsize   = 500;        #structure variants cutoff for those abnormal large delete/insertion
my $A                 = 3.87;       #
my $B                 = 174.28;     #
my $advance          = 0;          #if set 1, need to remove those mito reads could be remapped to non-mitochondria human genome
my $bwaindex          = undef;
my $bwa               = "bwa";      #supporse you have bwa install in your $PATH
my $help              =0;
#========End defining other variables======================#

#========Begin assiging values to variables================#
unless (
    GetOptions(
        "i=s"   => \$inbam1,
        "j=s"   => \$inbam2,
        "t=i"   =>\$type,
#        "b!"    => \$savebam,
#        "a!"    => \$saveallelecount,
        "ch!"   => \$producecircosplot,
        "hp=i"  => \$hp,
        "ha=i"  => \$ha,
        "d=i"   => \$depth,
        "A!"    => \$isall,
        "mmq=i" => \$mmq,
        "mbq=i" => \$mbq,
        "sb=i"  => \$sb,
        "cn!"   => \$cn,
        "sp=i"  => \$sp,
        "sa=i"  => \$sa,
        "cs!"   => \$cs,
#        "L=s"   => \$regionbed,
        "r=s"   => \$inref,
        "R=s"   => \$outref,
        "QC!"   => \$qc,
        "t=i"   => \$type,
        "str=i" => \$str,
        "strf=i"=> \$strflagmentsize,
        "alpha=f"   => \$A,
        "beta=f"   => \$B,
        "bwa=s" =>\$bwa,
        "bwaindex=s"=>\$bwaindex,
        "samtools=s"=>\$samtools,
        "advance"=>\$advance,
        "h|help"=>\$help
    )
  )
{
    print $!, "\n";
    _usage(1);
}
#==========End assiging values to variables================#

#==========Begin Main Program==============================#
if ($help){
    _usage();
    exit(0);
}
_check();
_initialVal();
_print_analysis_steps();
_main();
_make_report();
#============End main part================================#

#============Begin END section============================#
END {
    if($?){
        print "Program exist with Error \n";
    }else{
        my $interval = time() - $starttime;
        #my $interval=192232;
        my $hours = int( $interval / 3600 );  # calculate precise, and then floor it
        my $minutes = int( ( $interval - $hours * 3600 ) / 60 );
        my $seconds = $interval % 60;
        my $time    = sprintf( "%dh:%02dm:%02ds", $hours, $minutes, $seconds );
        print "Total Running Time: $time\n";
    }
}
#===========End end section===============================#

#============Begin defining functions=====================#
# Checking input parameters
sub _check{
    #1) inbam1 necessary
    if ( !defined($inbam1) ) {
        _error("input bam file (-i) has not been defined yet\n");
        _usage(1);
    }
    else {
        if ( !-e $inbam1 ) {
            _error("input bam file (-i) '$inbam1' does not exists\n");
            _usage(1);
        }
        #Conver to abs_path
        $isbam = _is_file_binary($inbam1);
        unless($isbam){
              _error("input bam file (-i) '$inbam1' is not in bam format\n") ;
            _usage(1);
        }
         $inbam1 = abs_path($inbam1) or die $!;
    }

    #2) inbam2 checking
    if ( defined($inbam2) ) {
        if ( !-e $inbam2 ) {
            _error("input bam file2 (-j) '$inbam2' does not exists\n");
            _usage(1);
        }
         $isbam = _is_file_binary($inbam2);
        unless($isbam){
              _error("input bam file (-i) '$inbam1' is not in bam format\n");
            _usage(1);
        }
        $inbam2 = abs_path($inbam2) or die $!;
    }

    #3) region file checking
    if ( defined($regionbed) ) {
        if ( !-e $regionbed ) {
            _error("bed file (-L) does not exists\n");
            _usage(1);
        }
        $regionbed = abs_path($regionbed);
    }

    #4) inref and outref checking
    #print $inref,"\n";
    if ( $inref ne 'hg19' && $inref ne 'rCRS' ) {
        _error(
    "The reference used in the bam file (-r) should be either hg19 or rCRS, yours is '$inref'\n"
        );
        _usage(1);
    }
    if ( $outref ne 'hg19' && $outref ne 'rCRS' ) {
        _error(
    "The reference used in the output file (-R) should be either hg19 or rCRS, yours is '$outref'\n"
        );
        _usage(1);
    }

    if ($cn && $type==4){
        _error("Your input bam is mitochondria only (-t), could not conduct the relative cony number estimation (-cn)\n");
        _usage(1);
    }

    if($type !=1 && $type!=2 && $type !=3 && $type!=4){
        _error("Input bam file type (-t) should be 1,2,3 or 4\n");
        _usage(1);
    }
}

# Initial some variables after checking the parameters
sub _initialVal{
    ($folder, undef, undef) = fileparse( $inbam1, qr/\.[^.]*/ );
    
    if ( -d $folder ) {
        _warn("folder $folder already exists,will delete it first");
        rmtree($folder) or die $!;
    }
    mkdir $folder or die $!;
    chdir $folder;

    if ( !$regionbed ) {
        $regionbed = "m.bed";
    }
    
    $mitoanno->build($outref);
    $mitoend=$acceptedgenomelength{$inref};
}

#Print out steps will be analyzed.
sub _print_analysis_steps{
    print "=" x 50, "\n";
    print "Steps will be run:\n";
    my $index=1;
    
    if($inbam2){
        print "    ",$index,".1,Extracting reads in mitochondria from '$inbam1' (Output: $mitobam1)\n";
        print "    ",$index++,".2,Extracting reads in mitochondria from '$inbam2' (Output: $mitobam2)\n";
        
        print "    ",$index,".1,Checking Reads number in '$mitobam1'\n";
        print "    ",$index++,".2,Checking Reads number in '$mitobam2'\n";
        
        print "    ",$index++,",Moving mitochondrial genome '".$mitogenome{$inref}. "' into $folder (Output:$reference)\n";
        
        print "    ",$index,".1,Pileuping '$mitobam1' (Output:$mitopileup1)\n";
        print "    ",$index++,".2,Pileuping '$mitobam2' (Output:$mitopileup2)\n";
        
        print "    ",$index,".1,Parsing pileup file of '$mitopileup1' (Output: $mitobasecall1)\n";
        print "    ",$index++,".2,Parsing pileup file of '$mitopileup2' (Output: $mitobasecall2)\n";
        
        if($qc){
            print "    ",$index,".1,Getting quality metrics on '$mitobam1' \n";
            print "    ",$index++,".2,Getting quality metrics on '$mitobam2' \n";
        }
        
        print "    ",$index,".1,Detecting heteroplasmy from '$mitobam1' (Output: $mitoheteroplasmy1)\n";
        print "    ",$index++,".2,Detecting heteroplasmy from '$mitobam2' (Output: $mitoheteroplasmy2)\n";
        
        print "    ",$index,".1,Detecting structure variants from '$mitobam1' (Output: $mitostructure1 | $mitostructuredeletion1)\n";
        print "    ",$index++,".2,Detecting structure variants from '$mitobam2' (Output: $mitostructure2 | $mitostructuredeletion2)\n";
        
        print "    ",$index++,",Detecting somatic mutations (Output: $mitosomatic)\n";
        
         if($cn && $type !=4){
            print "    ",$index,".1,Estimating relative copy number of '$mitobam1' (Output: $mitocnv1)\n";
            print "    ",$index++,".2,Estimating relative copy number of '$mitobam2' (Output: $mitocnv2)\n";
        }
    }else{
        print "    ",$index++,",Extracting reads in mitochondria from '$inbam1' (Output: $mitobam1)\n";
        print "    ",$index++,",Checking Reads number in '$mitobam1'\n";
        print "    ",$index++,",Moving mitochondrial genome '".$mitogenome{$inref}. "' into $folder (Output:$reference)\n";
        print "    ",$index++,",Pileuping '$mitobam1' (Output:$mitopileup1)\n";
        print "    ",$index++,",Parsing pileup file of '$mitopileup1' (Output: $mitobasecall1)\n";
        if($qc){
            print "    ",$index++,",Getting quality metrics on '$mitobam1' \n";
        }
        print "    ",$index++,",Detecting heteroplasmy from '$mitobam1' (Output: $mitoheteroplasmy1)\n";
        print "    ",$index++,",Detecting structure variants from '$mitobam1' (Output: $mitostructure1 | $mitostructuredeletion1)\n";
        
        if($cn && $type!=4){
            print "    ",$index++,",Estimating relative copy number of '$mitobam1' (Output: $mitocnv1)\n";
        }
     }
    print "    ",$index++,",Generating report (Output: $mitoreport)\n";
    print "=" x 50, "\n";
}


sub _main{
    print "=" x 50, "\n";
    print "Start analyzing:\n";
    my $index=1;
    _get_mitochondrial_bed( $inbam1, $regionbed );
    
    $ischr=_determine_chr_from_bam($inbam1,$isbam);
    #Assign values to 
    if($type==1 || $type==3){
        $totalbases=$totalexonbases;
        if($ischr){
            $totalbed=$exonbed{'withchr'};
        }else{
            $totalbed=$exonbed{'withoutchr'};
        }
    }else{
        $totalbases=$totalgenomebases;
        if($ischr){
            $totalbed=$genomebed{'withchr'};
        }else{
            $totalbed=$genomebed{'withoutchr'};
        }
    }
    
    if($inbam2){
        _info($index.".1,Extracting reads in mitochondria from '$inbam1' (Output: $mitobam1)");
        
        if($advance){
            _get_mitochondrial_bam_advance($inbam1,$isbam,$regionbed,$mmq,$mitobam1);
        }else{
            _get_mitochondrial_bam( $inbam1, $isbam, $regionbed, $mmq, $mitobam1 );
        }   
        _info($index++.".2,Extracting reads in mitochondria from '$inbam2' (Output: $mitobam2)");
        if($advance){
            _get_mitochondrial_bam_advance( $inbam2, $isbam, $regionbed, $mmq, $mitobam2 );
        }else{
            _get_mitochondrial_bam( $inbam2, $isbam, $regionbed, $mmq, $mitobam2 );
        }
               
        _info($index.".1,Checking Reads number in '$mitobam1'",1);
        _print_read($mitobam1);
         _info($index++.".2,Checking Reads number in '$mitobam2'",1);
         _print_read($mitobam2);
        
        _info($index++.",Moving mitochondrial genome '".$mitogenome{$inref}. "' into $folder (Output:$reference)");
        _move_mitogenome( $mitogenome{$inref}, $regionbed, $reference );
        
        _info($index.".1,Pileuping '$mitobam1' (Output:$mitopileup1)");
        _pileup( $mitobam1, $isbam, $mbq, $regionbed, $reference, $mitopileup1 );
        _info($index++.".2,Pileuping '$mitobam2' (Output:$mitopileup2)");
        _pileup( $mitobam2, $isbam, $mbq, $regionbed, $reference, $mitopileup2 );       
        
        _info($index.".1,Parsing pileup file of '$mitopileup1' (Output: $mitobasecall1)");
        _parse_pileup( $mitopileup1, $mbq, $mitooffset1, $mitobasecall1 );
        _info($index++.".2,Parsing pileup file of '$mitopileup2' (Output: $mitobasecall2)");
        _parse_pileup( $mitopileup2, $mbq, $mitooffset2, $mitobasecall2 );
        
        if($qc){
            _info($index.".1,Getting quality metrics on '$mitobam1' ");
            $mitooffset1 = _mito_qc_stat(
            $mitobam1,                  $isbam,
            $mitopileup1,               $regionbed,
            $perbasequality_figure1,    $mappingquality_figure1,
            $depthdistribution_figure1, $templatelengthdistribution_figure1,
            $perbasequality_table1,     $mappingquality_table1,
            $depthdistribution_table1,  $templatelengthdistribution_table1,
            $percentofbasepairscovered_table1
            );
            _info($index++.".2,Getting quality metrics on '$mitobam2' ");
            $mitooffset2 = _mito_qc_stat(
            $mitobam2,                  $isbam,
            $mitopileup2,               $regionbed,
            $perbasequality_figure2,    $mappingquality_figure2,
            $depthdistribution_figure2, $templatelengthdistribution_figure2,
            $perbasequality_table2,     $mappingquality_table2,
            $depthdistribution_table2,  $templatelengthdistribution_table2,
            $percentofbasepairscovered_table2
            );
        }
        
        _info($index.".1,Detecting heteroplasmy from '$mitobam1' (Output: $mitoheteroplasmy1)");
        _determine_heteroplasmy( $mitobasecall1, $hp, $ha, $isall, $sb,$mitoheteroplasmy1 );
        _info($index++.".2,Detecting heteroplasmy from '$mitobam2' (Output: $mitoheteroplasmy2)");
        _determine_heteroplasmy( $mitobasecall2, $hp, $ha, $isall, $sb,$mitoheteroplasmy2 );
        
        #Plot heteroplasmy ciros plot
        if($producecircosplot){
            $circos->build($outref);
            $circos->circosoutput($mitocircosheteroplasmyfigure1);
            $circos->configoutput($mitocircosheteroplasmyconfig1);
            $circos->datafile($mitoheteroplasmy1);
            $circos->textoutput($mitoheteroplasmytextoutput1);
            $circos->scatteroutput($mitoheteroplasmyscatteroutput1);
            $circos->cwd(getcwd()."/circos");
            $circos->prepare("heteroplasmy");
            $circos->plot();
            
            #For mito2
            $circos->circosoutput($mitocircosheteroplasmyfigure2);
            $circos->configoutput($mitocircosheteroplasmyconfig2);
            $circos->datafile($mitoheteroplasmy2);
            $circos->textoutput($mitoheteroplasmytextoutput2);
            $circos->scatteroutput($mitoheteroplasmyscatteroutput2);
            $circos->prepare("heteroplasmy");
            $circos->plot();
            
        }
        
        _info($index.".1,Detecting structure variants from '$mitobam1' (Output: $mitostructure1 | $mitostructuredeletion1)");
        _structure_variants( $inbam1, $isbam, $mmq, $regionbed, $str,$strflagmentsize,$mitostructure1,$mitostructuredeletion1);
        _info($index++.".2,Detecting structure variants from '$mitobam2' (Output: $mitostructure2 | $mitostructuredeletion2)");
        _structure_variants( $inbam2, $isbam, $mmq, $regionbed, $str,$strflagmentsize,$mitostructure2,$mitostructuredeletion2);
        
        _info($index++.",Detecting somatic mutations (Output: $mitosomatic)");
         _determine_somatic( $mitobasecall1, $mitobasecall2, $sp, $sa, $isall,$mitosomatic );
         
         if($cs){
            $circos->build($outref);
            $circos->changeconfig("somatic");
            $circos->circosoutput($mitocircossomaticfigure);
            $circos->configoutput($mitocircossomaticconfig);
            $circos->datafile($mitosomatic);
            $circos->textoutput($mitosomatictextoutput);
            $circos->cwd(getcwd()."/circos");
            $circos->prepare("somatic");
            $circos->plot();
         }
         
          if($cn && $type!=4){
            _info($index++.",Estimating relative copy number of '$mitobam1' (Output: $mitocnv1)");
            _wrap_mito_cnv($mitobam1,$inbam1,$mitobases,$totalbases,$isbam,$mbq,$mmq,$totalbed,$mitodepth1,$sampledepthi,$mitocnv1);
            _info($index++.",Estimating relative copy number of '$mitobam2' (Output: $mitocnv2)");
            _wrap_mito_cnv($mitobam2,$inbam2,$mitobases,$totalbases,$isbam,$mbq,$mmq,$totalbed,$mitodepth2,$sampledephtj,$mitocnv2);
         }
        

    }else{
        _info($index++.",Extracting reads in mitochondria from '$inbam1' (Output: $mitobam1)");
        if($advance){
            _get_mitochondrial_bam_advance( $inbam1, $isbam, $regionbed, $mmq, $mitobam1 );
        }else{
            _get_mitochondrial_bam( $inbam1, $isbam, $regionbed, $mmq, $mitobam1 );
        }
        
        _info($index++.",Checking Reads number in '$mitobam1'",1);
        _print_read($mitobam1);
        
        _info($index++.",Moving mitochondrial genome '".$mitogenome{$inref}. "' into $folder (Output:$reference)");
        _move_mitogenome( $mitogenome{$inref}, $regionbed, $reference );
        
        _info($index++.",Pileuping '$mitobam1' (Output:$mitopileup1)");
        _pileup( $mitobam1, $isbam, $mbq, $regionbed, $reference, $mitopileup1 );
        
        _info($index++.",Parsing pileup file of '$mitopileup1' (Output: $mitobasecall1)");
        _parse_pileup( $mitopileup1, $mbq, $mitooffset1, $mitobasecall1 );
        if($qc){
            _info($index++.",Getting quality metrics on '$mitobam1' ");
            $mitooffset1 = _mito_qc_stat(
            $mitobam1,                  $isbam,
            $mitopileup1,               $regionbed,
            $perbasequality_figure1,    $mappingquality_figure1,
            $depthdistribution_figure1, $templatelengthdistribution_figure1,
            $perbasequality_table1,     $mappingquality_table1,
            $depthdistribution_table1,  $templatelengthdistribution_table1,
            $percentofbasepairscovered_table1
            );
        }
        _info($index++.",Detecting heteroplasmy from '$mitobam1' (Output: $mitoheteroplasmy1)");
        _determine_heteroplasmy( $mitobasecall1, $hp, $ha, $isall, $sb,$mitoheteroplasmy1 );
        if($producecircosplot){
            $circos->build($outref);
            $circos->circosoutput($mitocircosheteroplasmyfigure1);
            $circos->configoutput($mitocircosheteroplasmyconfig1);
            $circos->datafile($mitoheteroplasmy1);
            $circos->textoutput($mitoheteroplasmytextoutput1);
            $circos->scatteroutput($mitoheteroplasmyscatteroutput1);
            $circos->cwd(getcwd()."/circos");
            $circos->prepare("heteroplasmy");
            $circos->plot();
         }
        
        _info($index++.",Detecting structure variants from '$mitobam1' (Output: $mitostructure1 | $mitostructuredeletion1)");
        _structure_variants( $inbam1, $isbam, $mmq, $regionbed, $str,$strflagmentsize,$mitostructure1,$mitostructuredeletion1); 
    
         if($cn && $type!=4){
            _info($index++.",Estimating relative copy number of '$mitobam1' (Output: $mitocnv1)");
            
            _wrap_mito_cnv($mitobam1,$inbam1,$mitobases,$totalbases,$isbam,$mbq,$mmq,$totalbed,$mitodepth1,$sampledepthi,$mitocnv1);
         }
    }
    _info($index++.",Generating report (Output: $mitoreport)");
    print "=" x 50, "\n";
    
}


#To determine wether the chromosome is named with prefix chr or not
#Will use the top 100 of mapped reads to determine
#Return value 0 or 1; 0 means without prefix chr while 1 means with prefix chr
sub _determine_chr_from_bam{
    my ($inbam,$isbam) = @_;
    my $flag=0; #default is not with chr
    my $count=0;
    my $max_read=100;
    if ($isbam) {
        open( IN, "$samtools view $inbam |" )
          or die $!;
    }
    else {
        open( IN, "$samtools view -S $inbam |" )
          or die $!;
    }
    
    while(<IN>){
        my @a = split "\t";
        if($a[2]=~/chr/){
            $flag=1;
            last;
        }
        $count++;
        last if ($count>$max_read);
    }
    close IN;
    return $flag;
}

#Used in checking the number of reads in mitochondria
#print out the number or otherwise exist the program
sub _print_read{
    my($mitobam)=@_;
    my $num1 = _number_of_reads_in_bam($mitobam);
    if ( $num1 != 0 ) {
        print " n=", $num1, "\n";
    }
    else {
        print "\n";
        _error("No reads in $mitobam \n",1);
        exit(1);
    }
}

#Get total lines(reads) in a bam file
#Parameter
#$inbam
sub _number_of_reads_in_bam {
    my ($inbam) = @_;
    my $num = 0;
    $num = `$samtools view $inbam|wc -l`;
    $num =~ s/\r|\n//g;
    return $num;
}

#Accoring to the input mitochondrial bam file, detect structure variants from this bam file. 
#Parameters:
#$inmitobam: alignment result in mitochrondrial
#$isbam: whether the input is in bam, otherwise it's sam 
#$mmq: minimal mapping quality
#$regionbed: a bed file containing analyzed regions
#$str: structure cutoff, means at least 'str' spanning reads supporting this variants
#$strflagmentsize. large deletion cutoff
#$outfile: Store the structure variants in the output file (discordantly mapped mates)
#$outfile2: store large deletion variants in the output file
sub _structure_variants {
    my ( $inmitobam, $isbam, $mmq, $regionbed, $str,$strflagmentsize, $outfile,$outfile2 ) = @_;
    my %temphash;
    if ($isbam) {
        open( IN, "$samtools view -L $regionbed -q $mmq $inmitobam |" )
          or die $!;
    }
    else {
        open( IN, "$samtools view -S -L $regionbed -q $mmq $inmitobam |" )
          or die $!;
    }
    open( OUT, ">$outfile" ) or die $!;
    open(OUT2,">$outfile2") or die $!;
    print OUT join "\t",
      (
        "#MitoChr",       "MitoPos",
        "VarChr",         "VarPos",
        "SupportedReads", "MeanMappingQuality",
        "Mito.gene","Mito.genedetail\n"
      );
    while (<IN>) {
        my (
            $qname, $flag,  $rname, $pos, $mapq, $cigar,
            $rnext, $pnext, $tlen,  $seq, $qual, $others
        ) = split "\t";

        if ( $rnext ne "*" && $rnext ne "=" && $rnext ne $rname ) {
            #print OUT $_;
            my $k = join "\t", ( $rname, $pos, $rnext, $pnext );
            $temphash{$k}->{'count'}++;
            push @{ $temphash{$k}->{'mapping'} }, $mapq;
        }
        if($tlen=~/\d+/ && abs($tlen)>$strflagmentsize){ #output candidate large delete changes
            print OUT2 $_;
        }
    }
    close IN;
            
    foreach my $k ( sort keys %temphash ) {
        if ( $temphash{$k}->{count} >= $str ) {
            my($rnametmp,$loc,$rnexttmp,$pnexttmp)= split "\t",$k;
            my $convertloc=$loc;
            #my $convertref=$result{$loc}->{'ref'};
            if($inref ne $outref){ #need to convert genome location
                  if($inref eq 'hg19' && $outref eq 'rCRS'){
                        $convertloc=$convert->hg19TorCRS($loc);
                        #$convertref=$convert->rCRSref($convertloc);
                  }else{
                        $convertloc=$convert->rCRSTohg19($loc);
                        #$convertref=$convert->hg19ref($convertloc);
                  }
                  next if (!defined($convertloc));  #can't be mapped between different assembly
            }
            
            print OUT join "\t",
              (
                "MT",$convertloc,$rnexttmp,$pnexttmp,
                $temphash{$k}->{'count'},
                _mymean( $temphash{$k}->{'mapping'} )
              );
            my ($genetmp,$genedetailtmp,undef,undef) = $mitoanno->annotation($convertloc);
            print OUT "\t",$genetmp,"\t",$genedetailtmp;
            print OUT "\n";
        }
    }
    close OUT;
}

sub _mymean {
    my ($arrayref) = @_;
    if ( scalar( @{$arrayref} ) == 0 ) {
        return 0;
    }
    else {
        my $t = 0;
        foreach ( @{$arrayref} ) {
            $t += $_;
        }
        return _formatnumeric($t / scalar( @{$arrayref} ));
    }

}

#Use samtools flagstat to get the total reads, mapped reads, and other inofrmation
#Parameter: $inbam
=example output
46123 + 0 in total (QC-passed reads + QC-failed reads)
12181 + 0 duplicates
46110 + 0 mapped (99.97%:nan%)
46123 + 0 paired in sequencing
23191 + 0 read1
22932 + 0 read2
41608 + 0 properly paired (90.21%:nan%)
45508 + 0 with itself and mate mapped
602 + 0 singletons (1.31%:nan%)
2923 + 0 with mate mapped to a different chr
2923 + 0 with mate mapped to a different chr (mapQ>=5)
=cut
sub _samtools_flagstat{
    my ($inbam) = @_;
    my ($total,$mapped)=(0,0);
    open(IN,"$samtools flagstat $inbam |") or die $!;
    my @array=<IN>;close IN;
    @array=map {/(^\d+)/;$1} @array;
    $total=$array[0];
    $mapped=$array[2];
    return ($total,$mapped);
}

sub _samtools_flagstat_use_pipe{
    my($inbam,$isbam,$mmq)=@_;
    my ($total,$mapped)=(0,0);
    if($isbam){
        open(IN,"$samtools view -u -q $mmq $inbam|$samtools  flagstat - |") or die $!;
    }else{
        open(IN,"$samtools view -Su -q $mmq $inbam|$samtools  flagstat - |") or die $!;
    }
    my @array=<IN>;close IN;
    @array=map {/(^\d+)/;$1} @array;
    $total=$array[0];
    $mapped=$array[2];
    return ($total,$mapped);
}




#Call cnv and write out by wrap by read and by depth
sub _wrap_mito_cnv{
       my ($mitobam,$totalbam,$mitobases,$totalbases,$isbam,$mbq,$mmq,$bed,$mitodepthoutput,$totaldepthoutput,$mitocnv) = @_;
       my ($mitomapped,$totalmapped,$byreadcnv)=_mito_cnv_by_reads($mitobam,$totalbam,$isbam,$mmq);
       my ($mitoaveragedepth,$totalaveragedepth,$bydepthcnv)=_mito_cnv_by_depth($mitobam,$totalbam,$mitobases,$totalbases,$isbam,$mbq,$mmq,$bed,$mitodepthoutput,$totaldepthoutput);
       open(OUT,">$mitocnv") or die $!;
       print OUT "Method\tmito.reads/mito.ave.depth\ttotal.reads/total.ave.depth\tCNV\n";
       print OUT join "\t",("ByRead",$mitomapped,$totalmapped,$byreadcnv);
       print OUT "\n";
       print OUT join "\t",("ByDepth",$mitoaveragedepth,$totalaveragedepth,$bydepthcnv);
       print OUT "\n";
       close OUT;
}

#Relative Copy Number Estimate
#Parameter:
#$mitobam: alignment mapped to mitochondria
#$totalbam: alignment mapped to whole genome
#$isbam: bam or sam
#$mmq: mimimal mapping quality

sub _mito_cnv_by_reads{
    my ($mitobam,$totalbam,$isbam,$mmq) = @_;
    my ($mitototal,$mitomapped) = _samtools_flagstat_use_pipe($mitobam,$isbam,$mmq);
    my ($totaltotal,$totalmapped) = _samtools_flagstat_use_pipe($totalbam,$isbam,$mmq);
    #Remove the mapped reads in the mito
    $totalmapped=$totalmapped-$mitomapped;
    return ($mitomapped,$totalmapped,_divide($mitomapped,$totalmapped));
}


#For exom sequence bed is the region of exon while for whole genome, the bed is the region of all chromosome
sub _mito_cnv_by_depth{
    my ($mitobam,$totalbam,$mitobases,$totalbases,$isbam,$mbq,$mmq,$bed,$mitodepthoutput,$totaldepthoutput) = @_;
    my $com = "$samtools depth -q $mbq -Q mmq $mitobam > $mitodepthoutput";
    _run($com);
    $com = "$samtools depth -q $mbq -Q mmq -b $bed $totalbam >$totaldepthoutput";
    _run($com);
    
    #Read the 
    my $mitodepth = _read_and_sum_depth($mitodepthoutput);
    my $totaldepth = _read_and_sum_depth($totaldepthoutput);
    
    my $mitoaveragedepth = _formatnumeric( $mitodepth/$mitobases);
    my $totalaveragedepth = _formatnumeric($totaldepth/$totalbases);
    return ($mitoaveragedepth,$totalaveragedepth,_divide($mitoaveragedepth,$totalaveragedepth));
}

sub _read_and_sum_depth{
    my ($in) = @_;
    my $total=0.0;
    open(IN,$in) or die $!;
    while(<IN>){
        s/\r|\n//g;
        my (undef,undef,$a) = split /\t/;
        $total+=$a;
    }
    return $total;
}

sub _divide{
    my ($a,$b) = @_;
    if($b == 0){
        return 'NA';
    }else{
        return _formatnumeric($a/$b);
    }
}


#Move the mitochondrial genome into the result folder, change the header (>) when necessary.
sub _move_mitogenome {
    my ( $infasta, $regionbed, $outfasta ) = @_;
    my $chromosome = "";
    open( IN, $regionbed ) or die $!;
    my $l = <IN>;
    ( $chromosome, undef ) = split /\s+/, $l;
    close IN;
    open( IN,  $infasta )     or die $!;
    open( OUT, ">$outfasta" ) or die $!;
    <IN>;
    print OUT ">$chromosome\n";

    while (<IN>) {
        print OUT $_;
    }
    close IN;
    close OUT;
}


#Given a input file, determine whether it's a binary format or text format
#return 1 if binary; Default is 0 Means text
sub _is_file_binary {
    my ($infile) = @_;
    my $flag = 0;
    if ( -B $infile ) {
        #print "binary\n";
        $flag = 1;
    }

    #if(-T $infile){
    #    print "Text\n";
    #}
    return $flag;
}

#Print out usage,
#Parameters:
#$flag, if $flag is set and it will exit the program
sub _usage {
    my ($flag) = @_;
=h com
@----------------------------------------------------------@
|        mitoSeek     |     v1.01      |   02/15/2013      |
|----------------------------------------------------------|
|     (C) 2013 Jiang Li, GNU General Public License, v2    |
|----------------------------------------------------------|
|  For documentation, citation & bug-report instructions:  |
|         https://github.com/riverlee/MitoSeek             |
@----------------------------------------------------------@
=cut
    my $usage=<<USAGE;
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
-beta  [float]          Shape1 parameter of Beta prior distribution, default is 174.28, which is estimated from 600 BRCA samples
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
-bwa [bwa]              Tell where is the bwa program, default value is 'bwa' which is your \$PATH
-bwaindex [bwaindex]    Tell where is the bwa index of non-mitochondrial human genome, no default value
-advance                Will get mitochondrial reads in an advanced way, generally followed by 1) Initially extract mitochrodrial reads from 
                        a bam file, then 2) remove those could be remapped to non-mitochondrial human genome by bwa. Advanced extraction needs 
                        -bwaindex option. Default extraction without removing step.


\@------------------------------------------------------------@
|       Statistical framework for heteroplasmy detection     |
|------------------------------------------------------------|
|              Fisher test for heterplasmy                   |
|------------------------------------------------------------|
|                           major   minor                    |
|               observed    n11     n12 | n1p                |
|               expected    n21     n22 | n2p                |
|                          -----------------                 |
|                           np1     np2   npp                |
|   n21 = (n11+n12)*(1-hp/100) in which hp is defined by -hp |
|   n22 = (n11+n12)*hp/100  in which hp is defined by -hp    |
|------------------------------------------------------------|
|                Empirical Bayesian method                   |
|------------------------------------------------------------|
|                          _Inf                              |
|                         /                                  |
|           probability = | f(x)dx                           |
|                       _/p                                  |
|   probablity is the calculus of f(x) from p to Inf         |
|   x=hp/100 in which hp is defined by -hp                   |
|   f(x) = 1/beta(b+A,a+B)*x^(A+b-1)*(1-x)^(B+a-1)           |
|   in which a/b is the number of major/minor allele,        |
|   A and B are estimated from 600 BRCA samples.             |
|------------------------------------------------------------|
|    For documentation, citation & bug-report instructions:  |
|           https://github.com/riverlee/MitoSeek             |
\@------------------------------------------------------------@

USAGE

    print $usage;
    exit(1) if ($flag);
}

#Get a serials metrics to quality control alignment on mitochondrial(
#1)Insert size distribution if paired,
#2)depth distribution,
#3)base quality,
#4)mapping quality
#5)what is its pair which mapped to other chromosome (table))
#Parameters:
#$inmitobam  input bam or sam format file
#$isbam  0/1, 0 indicates it's in sam format while 1 indicates a bam format

sub _mito_qc_stat {
    my (
        $inmitobam,                $isbam,
        $pileupfile,               $regionbed,
        $perbasequality_figure,    $mappingquality_figure,
        $depthdistribution_figure, $templatelengthdistribution_figure,
        $perbasequality_table,     $mappingquality_table,
        $depthdistribution_table,  $templatelengthdistribution_table,
        $percentofbasepairscovered_table
    ) = @_;
    my $maxreads = 100000; #In case some alignment on mitochondrial is extremly large that will take too much memory and then crash ymy server.
    $maxreads = 10 if ($debug);

    if ($isbam) {
        open( IN, "$samtools view -L $regionbed $inmitobam |" ) or die $!;
    }
    else {
        open( IN, "$samtools view -S -L $regionbed $inmitobam |" ) or die $!;
    }

#Use the top 5 reads to determine whether it is pair-end reads and get the quality encoding
    my $isparied = 0;
    my $offset   = 33;
    my $encoding = "Sanger / Illumina 1.9";    #Default
    my $top      = 5;
    my $count    = 1;
    my @scores;    #store the qual ascii value to determine quality encoding
    while (<IN>) {
        s/\r|\n//g;
        $count++;
        my (
            $qname, $flag,  $rname, $pos, $mapq, $cigar,
            $rnext, $pnext, $tlen,  $seq, $qual, $others
        ) = split "\t";
        if ( $flag & $flag_paired ) {
            $isparied = 1;
        }
        foreach my $s ( split //, $qual ) {
            push @scores, ord($s);
        }
        last if ( $count > $top );
    }
    @scores = sort { $a <=> $b } @scores;
    my $lowestChar = $scores[0];
    if ( $lowestChar < 33 ) {
        #Generate a warnings
        _warn(
"No known encodings with chars < 33 (Yours was $lowestChar), However, we will use Sanger / Illumina 1.9 encoding instead"
        );
    }
    elsif ( $lowestChar < 59 ) {
        #it's Sanger/Illumina 1.9 encoding, don't need to change
    }
    elsif ( $lowestChar < 64 ) {
        $offset   = 59;
        $encoding = "Illumina <1.3";
    }
    elsif ( $lowestChar == 65 ) {
        $offset   = 64;
        $encoding = "Illumina 1.3";
    }
    elsif ( $lowestChar <= 126 ) {
        $offset   = 64;
        $encoding = "Illumina 1.5";
    }
    else {
        _warn(
"No known encodings with chars > 126 (Yours was $lowestChar), However, we will Sanger / Illumina 1.9 encoding instead"
        );
    }

    #move to the start of file and begin to get the stat
    seek( IN, 0, 0 );
    $count = 1;

    my @tlen
      ; #if paired, the Insert size between two pairs is the tlen-read1_len-read2_len, if single, the tlen is equal to 0
    my @perbasequality;
    my @persequencequality;
    my @mapquality;

    while (<IN>) {
        s/\r|\n//g;
        $count++;
        my (
            $qname, $flag,  $rname, $pos, $mapq, $cigar,
            $rnext, $pnext, $tlen,  $seq, $qual, $others
        ) = split "\t";
        last if ( $count > $maxreads );

        #Store template length
        if ($isparied) {
            #Only store leftmost fragment's
            if ( $tlen =~ /^\d/ && $tlen != 0 ) {
                push @tlen, $tlen;
            }
        }

        #Store mapping quality
        push @mapquality, $mapq if ( $mapq =~ /^\d/ );

        #Store perbasequality and persequencequality
        my @qual = split //, $qual;
        my $totalscore = 0;
        for ( my $i = 0 ; $i < @qual ; $i++ ) {
            my $s = ord( $qual[$i] ) - $offset;
            $totalscore += $s;
            push @{ $perbasequality[1]->[$i] }, $s;
        }
        if (@qual) {
            push @persequencequality, $totalscore / scalar(@qual);
        }
    }
    @{ $perbasequality[0] } = 1 .. scalar( @{ $perbasequality[1] } );

    #Plot the figures;
    #1) Template length density distribution
    _density(
        \@tlen, 800, 600, 'Length in bp', 'Density',
        'Template length distribution',
        $templatelengthdistribution_figure,
        0, undef, 30,$templatelengthdistribution_table
    ) if ($isparied);

    #2) persequencequality density distribution
    _density(
        \@persequencequality, 800, 600, 'Quality score',
        'Density', "Per sequence quality socre distribution ($encoding)",
        $mappingquality_figure, 0, undef, 30,$mappingquality_table
    );

    #3) perbasequality boxplot
    _boxplot(
        \@perbasequality, 800, 600, 'Bases',
        'Quality score',
        "Per base quality socre distribution ($encoding)",
        $perbasequality_figure, undef, undef, 30,$perbasequality_table
    );

    #Plot depth distribution
    my %position;
    open( IN, $regionbed ) or die $!;
    while (<IN>) {
        s/\r|\n//g;
        my ( $chr, $start, $end ) = split /\s+/;
        for ( my $i = $start + 1 ; $i <= $end ; $i++ ) {
            $position{$i} = 0;
            last if ( $inref eq 'hg19' && $i > 16571 );
            last if ( $inref eq 'rCRS' && $i > 16569 );
        }
    }
    close IN;

    #print scalar(keys %position),"\n";
    my @depth;
    open( IN, $pileupfile ) or die $!;
    while (<IN>) {
        my ( $chr, $pos, $ref, $depth, @others ) = split "\t";
        push @depth, $depth;
        $position{$pos} = $depth;
    }

    # print scalar(@depth),"\n";
    for ( my $i = 0 ;
        $i < ( scalar( keys %position ) - scalar(@depth) ) ; $i++ )
    {
        push @depth, 0;
        print $i, "\n";
    }
    _density( \@depth, 800, 600, 'Depth', 'Density', 'Depth distribution',
        $depthdistribution_figure, 0, undef, 30 ,$depthdistribution_table);

    open(PER,">$percentofbasepairscovered_table") or die $!;
    print PER "#Covered bases\tTotal bases\tCoverage\n";
    my $covered = grep {$_>0} @depth;
    my $total = scalar(@depth);
    print PER join "\t",($covered,$total,_formatnumeric(_divide($covered,$total)));
    print PER "\n";
    close PER;
    return ($offset);
}

#Function to plot boxplot
#Parameters:
sub _boxplot {
    my (
        $data,  $width,  $height, $xlabel, $ylabel,
        $title, $output, $xmin,   $xmax,   $nbin, $outputtable
    ) = @_;
    my $graph = new GD::Graph::boxplot( $width, $height );
    $graph->set(
        x_label           => $xlabel,
        y_label           => $ylabel,
        title             => $title,
        x_labels_vertical => 1,
        transparent       => 0,
        x_tick_number     => $nbin,
        x_number_format   => \&_number_format,
        x_label_position  => 1 / 2,
        upper_percent     => 75,
        lower_percent     => 25,
        do_stats          => 1,
        box_spacing       => 2,
        y_min_value       => 0,
        fov_const         => 0,
        t_margin          => 10,
        b_margin          => 10,
        l_margin          => 10,
        r_margin          => 20,
        do_stats              => 0,
        box_fill            =>1,
        dclrs             => ['lgreen','lred'],  #will be the color in the box of 25%-75%
        fgclr             => 'dblue' ,        #color for the outlier, etc
    ) or warn $graph->error;
    $graph->set_title_font(gdGiantFont);
    $graph->set_x_label_font(gdMediumBoldFont); 
    $graph->set_y_label_font(gdMediumBoldFont);
    #$graph->set_x_axis_fobbnnnnt(gdMediumBoldFont);
    $graph->set_y_axis_font(gdMediumBoldFont);
    $graph->set_values_font(gdMediumBoldFont);
    
    #Customed calculate the stat for boxplot
    my $newdata=[$data->[0],];
    for(my $i=0;$i<@{$data->[1]};$i++){
          my $stat = Statistics::Descriptive::Full->new();
          my $value=$data->[1]->[$i];
          my $j;	# declaration required for comparison below
          for($j=0; defined $value->[$j]; $j++){
            $stat->add_data($value->[$j]);
          }
         my $upper = $stat->percentile(75);
         my $lower = $stat->percentile(25 );
         my $meanv = $stat->mean();
         my $medianv = $stat->median();
         my $highest = $stat->max();
         my $lowest  = $stat->min();
         $newdata->[1]->[$i]=[$meanv,$lowest,$lower,$medianv,$upper,$highest];
    }
    
    # Output data for the plot
    open(TMP,">$outputtable") or die $!;
    print TMP "#Data for the boxplot\n";
    print TMP join "\t",("Stat",@{$newdata->[0]});
    print TMP "\n";
    print TMP join "\t",("Mean",map {$_->[0]} @{$newdata->[1]});
    print TMP "\n";
    print TMP join "\t",("Lowest",map {$_->[1]} @{$newdata->[1]});
    print TMP "\n";
    print TMP join "\t",("Per25",map {$_->[2]} @{$newdata->[1]});
    print TMP "\n";
    print TMP join "\t",("Median",map {$_->[3]} @{$newdata->[1]});
    print TMP "\n";
    print TMP join "\t",("Per75",map {$_->[4]} @{$newdata->[1]});
    print TMP "\n";
    print TMP join "\t",("Highest",map {$_->[5]} @{$newdata->[1]});
    print TMP "\n";
    close TMP;
    
    my $gd = $graph->plot($newdata) or die $graph->error;
    open( IMG, ">$output" ) or die $!;
    binmode IMG;
    print IMG $gd->png;
    close IMG;
}

#Functions to estimate the density using Statistics::KernelEstimation
#Parameters:
#$data --- array reference [1,24,5,7]
#$width ---width for the figure
#$height--- height for the figure
#$xlabel----x-axis label
#$ylabel----y-axis label
#$title-----title for the figure
#$output ---output figure name
#$xmin -----minimun value allowed at x-axis
#$xmax -----maximun value allowed at x-axis
#$bin  -----How many bins at x-axis,Default is 30
sub _density {
    my (
        $data,  $width,  $height, $xlabel, $ylabel,
        $title, $output, $xmin,   $xmax,   $nbin,$outputtable
    ) = @_;
    my $s = Statistics::KernelEstimation->new();
    foreach my $x ( @{$data} ) {
        $s->add_data($x);
    }
    my $w = $s->default_bandwidth();
    my ( $min, $max ) = $s->extended_range();

    #my ($min,$max) = (min(@{$data}), max(@{$data}));
    if ( defined($xmin) && $min < $xmin ) {
        $min = $xmin;
    }
    if ( defined($xmax) && $max > $xmax ) {
        $max = $xmax;
    }
    $nbin = 30 if ( !defined($nbin) || $nbin <= 0 );

    my @plot_data;
    

#If nbins is very large, this may cause a lots of bins if the range of $data is large, leading to ugly figure
    for ( my $x = $min ; $x <= $max ; $x += ( $max - $min ) / $nbin ) {

        # print $x, "\t", $s->pdf( $x, $w ), "\t", $s->cdf( $x, $w ), "\n";
        push @{ $plot_data[0] }, $x;
        push @{ $plot_data[1] }, $s->pdf( $x, $w );
    }
  
    #print $max,"\n";
    my $graph = new GD::Graph::lines( $width, $height );
    $graph->set(
        x_label           => $xlabel,
        y_label           => $ylabel,
        title             => $title,
        x_labels_vertical => 1,
        transparent       => 0,
        x_tick_number     => $nbin,
        x_min_value       => $min,
        x_max_value       => int($max) + 1,
        x_number_format   => \&_number_format,
        x_label_position  => 1 / 2,
        legend_placement  => 'RL',
        t_margin          => 10,
        b_margin          => 10,
        l_margin          => 10,
        r_margin          => 20,
        line_width        => 3
    ) or warn $graph->error;
    $graph->set_title_font(gdGiantFont);
    $graph->set_x_label_font(gdMediumBoldFont); 
    $graph->set_y_label_font(gdMediumBoldFont);
    $graph->set_x_axis_font(gdMediumBoldFont);
    $graph->set_y_axis_font(gdMediumBoldFont);
    $graph->set_values_font(gdMediumBoldFont);
   
    
    my $gd = $graph->plot( \@plot_data ) or die $graph->error;

    #Calculate the mean, median, sd from the $data
    my $stat = Statistics::Descriptive::Full->new();
    foreach my $tmp ( @{$data} ) {
        $stat->add_data($tmp);
    }
    my $mean   = sprintf( "%.3f", $stat->mean() );
    my $median = sprintf( "%.3f", $stat->median() );
    my $var    = sprintf( "%.3f", $stat->variance() );
    my $v25    = sprintf( "%.3f", $stat->percentile(25) );
    my $v75    = sprintf( "%.3f", $stat->percentile(75) );
    my $minv   = sprintf( "%.3f", min( @{$data} ) );
    my $maxv   = sprintf( "%.3f", max( @{$data} ) );
    my $text   = <<EOSTR;
Min. = $minv
25%   = $v25
Mean = $mean
Median = $median
75% = $v75
Max. = $maxv
Var = $var
EOSTR

    # Print out the data for the plot
    open(TMP,">$outputtable") or die $!;
    print TMP "#Summary of the data\n";
    print TMP join "\t",("#Min.","25%","Mean","Median","75%","Max.","Var\n");
    print TMP join "\t",("#$minv",$v25,$mean,$median,$v75,$maxv,$var);
    print TMP "\n";
    print TMP "#Data for density plot\n";
    print TMP "#X\tY(density)\n";
    for (my $i=0;$i<@{$plot_data[0]};$i++){
        print TMP join "\t",($plot_data[0]->[$i],$plot_data[1]->[$i]);
        print TMP "\n";
    }
    close TMP;
    
    my $gdat = GD::Text::Wrap->new(
        $gd,
        text       => $text,
        line_space => 4,
        align      => 'left',
        width      => 200,
        preserve_nl=>1
    );
    $gdat->set_font(gdMediumBoldFont);

    #$gdat->draw(10,140);
    $gdat->draw( $width / 2, $height / 5 );
    open( IMG, ">$output" ) or die $!;
    binmode IMG;
    print IMG $gd->png;
    close IMG;
}

#Function used to format value into integer (invoked by $graph->set)
#Parameters:
#$value-- input number
sub _number_format {
    my $value = shift;
    my $ret;

    if ( $value >= 0 ) {
        $ret = sprintf( "\%d", $value );
    }
    else {
        $ret = sprintf( "-\%d", abs($value) );
    }
    return $ret;
}

#Functions to plot histgram using GD::Graph::histogram; (Will cause problem when the ranges of $data is large, not too much space for the plot)
#Parameters:
#

sub _histogram {
    my ( $data, $height, $width, $xlabel, $ylabel, $title, $output ) = @_;
    my $graph = new GD::Graph::histogram( $width, $height );
    $graph->set(
        x_label           => $xlabel,
        y_label           => $ylabel,
        title             => $title,
        x_labels_vertical => 1,
        bar_spacing       => 0,
        shadow_depth      => 1,
        shadowclr         => 'dred',
        transparent       => 0,
    ) or warn $graph->error;
    my $gd = $graph->plot($data) or die $graph->error;
    open( IMG, ">$output" ) or die $!;
    binmode IMG;
    print IMG $gd->png;
    close IMG;
}

#Use the top 5 reads in the alignment result to determine whether its paired
sub _determine_paired_from_bam {
}

#Use samtools view -H read the header, and prepare the mithochondrial bed file.
sub _get_mitochondrial_bed {
    my ( $in, $outbed ) = @_;
    open( IN, "samtools view -H $in|" ) or die $!;
    my $flag = 0;
    my $m    = "";
    my $len  = "";
    while (<IN>) {
        if (/SN:M|SN:chrM/i) {
            if (/SN:(.*?)\s+LN:(\d+)/) {
                $m    = $1;
                $len  = $2;
                $flag = $1;
                last;
            }
        }
    }
    unless ($flag) {
        _error("No mithochondrial chromosome detected from the header");
        exit(1);
    }
   
    if ( $len == 16569 && $inref eq "hg19" ) {
        _warn(
"It seems like the mitochondrial reference used in the bam is rCRS, not hg19"
        );
    }

    if ( $len == 16571 && $inref eq "rCRS" ) {
        _warn(
"It seems like the mitochondrial reference used in the bam is hg19, not rCRS"
        );
    }
    
    if($len != 16569 && $len !=16571){
       _error("Mithochondrial genome length ($len) is invalid");
       exit(1);
    }
    open( OUT, ">$outbed" );
    print OUT join "\t", ( $m, 0, $len );
    close OUT;
}

#         word2   ~word2
#word1    n11      n12 | n1p
#~word1   n21      n22 | n2p
#         --------------
#         np1      np2   npp
#a number of major allele
#b number of minor allele
#p percentage of determining a heteroplasmy
sub fisher_left{
    my ($a,$b,$p) = @_;
    if(@_ !=3){
        print  "usage fisher_left major minor percetage\n";
        exit(1);
    }
    my $npp = $a+$b+$a+$b;
    my $n1p = $a+$b;
    my $np1 = $a+($a+$b)*(1-$p);
    my $n11 = $a;
    my $left_value = calculateStatistic( n11=>$n11,
         n1p=>$n1p,
         np1=>$np1,
         npp=>$npp);
    
   if(getErrorCode()){
        $left_value=1;
   }
   $left_value=0 if($left_value<0);
   $left_value=1 if($left_value>1);
   return($left_value);
}

#$major number of major allele
#$minor number of minor allele
#$x     default 0.01, equal to hp/100
#$len   to separate 0 to x into len
#$A     A=3.87
#$B     B=174.28
sub empirical_bayesian{
    my ($major,$minor,$x,$len,$A,$B) = @_;
    my $width=$x/($len-1);
    my @values;
    for my $i (1..$len){
        my $tmp = ($i-1)*$width;
        my $v = (1/beta($minor+$A,$major+$B)) * ($tmp**($A+$minor-1)) * ((1-$tmp)**($B+$major-1)) ;
        push @values,$v*$width;
    }
    #sum @values
    my $sum = 0;
    foreach (@values){
        $sum+=$_;
    }
    my $r=(1-$sum);
    $r=0 if ($r<0);
    $r=1 if ($r>1);
    return $r;
}

#convert a probablity into phred score
sub phred{
    my ($p) = @_;
    if($p<=0){
        return 255;
    }else{
        return -10*log10($p);
    }
}

sub log10{
    my $v=shift;
    return log($v)/log(10);
}

#Determine heteroplasmy mutations from a basecall format file with given parameters
#Parameters:
#$inbase  input file of allele count parsed from pileup file
#$hp      (-hp) Heteroplasmy threshold using percent alternative allele observed
#$ha      (-ha) Heteroplasmy threshold using allele observed
#$isall   (-A)  0/1, 1 denotes that the total read count is the total allele count of all allele observed,
#         while 0 indicates the total read count is the sum of major and minor allele counts, default=0
#$sb      (-sb) Remove all sites with strand bias score in the top %, default = 10
sub _determine_heteroplasmy {
    my ( $inbase, $hp, $ha, $isall, $sb, $outheteroplasmy ) = @_;
    open( IN, $inbase ) or die $!;
    <IN>;    #skip the header
             #my @result;
    my %result;
    my %sb;
    my %rawpvalue;
    my %empirical_probality;
    while (<IN>) {
        s/\r|\n//g;
        my (
            $chr,       $loc,       $ref,       $forward_A,
            $forward_T, $forward_C, $forward_G, $reverse_A,
            $reverse_T, $reverse_C, $reverse_G
        ) = split "\t";

        # Get the major allele and minor allele.
        my %atcg;
        $atcg{A} = $forward_A + $reverse_A;
        $atcg{T} = $forward_T + $reverse_T;
        $atcg{C} = $forward_C + $reverse_C;
        $atcg{G} = $forward_G + $reverse_G;
        my @atcg         = sort { $atcg{$b} <=> $atcg{$a} } keys %atcg; #large to small
        my $major_allele = $atcg[0];
        my $minor_allele = $atcg[1];
        my $totaldepth= $atcg{A} + $atcg{C} + $atcg{T} + $atcg{G} ;
        
        next if ($totaldepth <= $depth);  #only have sufficiant depth will conduct heteroplasmy detection

        # Added on 2013-02-12, fisher test (left)
        my $heteroplasmy_fisher_pvalue = fisher_left($atcg{$major_allele},$atcg{$minor_allele},$hp/100);
        $rawpvalue{$loc}=$heteroplasmy_fisher_pvalue;
        #Calculate heteroplasmy ratio
        my $heteroplasmy = 0;
        if ($isall) {
            $heteroplasmy =
              _formatnumeric($atcg{$minor_allele} /$totaldepth);
        }
        else {
            $heteroplasmy =
            _formatnumeric($atcg{$minor_allele} /
              ( $atcg{$major_allele} + $atcg{$minor_allele} ));
        }

        #Calculate 95% confidence interval
        my ( $lower, $upper ) = ( 0, 0 );
        if ( $heteroplasmy > 0 ) {
            my $n = 0;
            if ($isall) {
                $n = $totaldepth;
            }
            else {
                $n = $atcg{$major_allele} + $atcg{$minor_allele};
            }
            ( $lower, $upper ) = _ci( $heteroplasmy, $n );
        }

        #Stat values could be
        #0 (not a heteroplasmy mutation)
        #1 (show heteroplasmy, but not pass the cutoff) (hp,ha,depth)
        #2 (heteroplasmy mutation pass cutoff and not have strong strand bias)
        #3 (heteroplasmy mutation pass cutoff but show strong strand bias)
        my $stat = 0;
        if ( $heteroplasmy > 0 ) {

#determine 1 or 2 at this stage (for value 3, need to first calculate strand bias across all the site, and then modify those with stat=2, only if when $sb is not equal to 0)
            #if ( $heteroplasmy > $hp / 100 && $atcg{$minor_allele} > $ha ) {
            if ( $heteroplasmy > $hp / 100 && $atcg{$minor_allele} > $ha ) {  #use $depth which is not passed by function
                $stat = 2;
            }
            else {
                $stat = 1;
            }
        }

#If provided strand bias cutoff, need to re-loop the @result, change state=2 to 3 when necessary.
#        push @result,
#          [
#            $chr,          $loc,          $ref,
#            $forward_A,    $forward_T,    $forward_C,
#            $forward_G,    $reverse_A,    $reverse_T,
#            $reverse_C,    $reverse_G,    $stat,
#            $heteroplasmy, $lower,        $upper,
#            $major_allele, $minor_allele, $atcg{$major_allele},
#            $atcg{$minor_allele}
#          ];
        $result{$loc}->{'chr'}                = $chr;
        $result{$loc}->{'ref'}                = $ref;
        $result{$loc}->{'A'}                  = $forward_A;
        $result{$loc}->{'T'}                  = $forward_T;
        $result{$loc}->{'C'}                  = $forward_C;
        $result{$loc}->{'G'}                  = $forward_G;
        $result{$loc}->{'a'}                  = $reverse_A;
        $result{$loc}->{'t'}                  = $reverse_T;
        $result{$loc}->{'c'}                  = $reverse_C;
        $result{$loc}->{'g'}                  = $reverse_G;
        $result{$loc}->{'stat'}               = $stat;
        $result{$loc}->{'heteroplasmy'}       = $heteroplasmy;
        $result{$loc}->{'lower'}              = $lower;
        $result{$loc}->{'upper'}              = $upper;
        $result{$loc}->{'major_allele'}       = $major_allele;
        $result{$loc}->{'minor_allele'}       = $minor_allele;
        $result{$loc}->{'major_allele_count'} = $atcg{$major_allele};
        $result{$loc}->{'minor_allele_count'} = $atcg{$minor_allele};

        #calculate strand bias if $sb>0
        #if ( $sb > 0 ) {
            my (
                $forward_primary_allele_count,
                $forward_non_primary_allele_count,
                $reverse_primary_allele_count,
                $reverse_non_primary_allele_count
            ) = ( 0, 0, 0, 0 );
            my %tmp;
            $tmp{'A'}->{'forward'} = $forward_A;
            $tmp{'A'}->{'reverse'} = $reverse_A;
            $tmp{'T'}->{'forward'} = $forward_T;
            $tmp{'T'}->{'reverse'} = $reverse_T;
            $tmp{'C'}->{'forward'} = $forward_C;
            $tmp{'C'}->{'reverse'} = $reverse_C;
            $tmp{'G'}->{'forward'} = $forward_G;
            $tmp{'G'}->{'reverse'} = $reverse_G;

#Now  no matter isall=0 or 1, strand bias is calculated based on major allele and minor allele
            (
                $forward_primary_allele_count,
                $forward_non_primary_allele_count,
                $reverse_primary_allele_count,
                $reverse_non_primary_allele_count
              )
              = (
                $tmp{$major_allele}->{'forward'},
                $tmp{$minor_allele}->{'forward'},
                $tmp{$major_allele}->{'reverse'},
                $tmp{$minor_allele}->{'reverse'}
              );

            $sb{$loc} = _sb(
                $forward_primary_allele_count,
                $forward_non_primary_allele_count,
                $reverse_primary_allele_count,
                $reverse_non_primary_allele_count
            );
            
            $result{$loc}->{'sb'} = $sb{$loc};
        #} #sb calculate
    }

    #Re-loop the @result and write out
    #    if ( $sb > 0 ) {
    #        my @sbvalue = sort { $b <=> $a } values %sb;
    #        my $index = int( scalar(@sbvalue) * $sb / 100 ) - 1;
    #        $index = 0         if ( $index < 0 );
    #        $index = $#sbvalue if ( $index > $#sbvalue );
    #        my $cutoff = $sbvalue[$index];
    #        foreach my $arrayref (@result) {
    #            if (   $sb{ $arrayref->[1] } >= $cutoff
    #                && $arrayref->[11] == 2 )
    #            {
    #                $arrayref->[11] = 3;
    #            }
    #        }
    #    }
    if ( $sb > 0 ) { 
        my @sbvalue = sort { $a <=> $b } values %sb;   #small to large
        my $index = int( scalar(@sbvalue) * $sb / 100 ) - 1;
        $index = 0         if ( $index < 0 );
        $index = $#sbvalue if ( $index > $#sbvalue );
        my $cutoff = $sbvalue[$index];

        foreach my $loc ( keys %result ) {
            if ( $sb{$loc} >= $cutoff && $result{$loc}->{'stat'} == 2 ) {
                $result{$loc}->{'stat'} = 3;
            }
        }
    }

    #Added on 2013-02-12, multiple test
    my $adjustpref=BH(\%rawpvalue);

    #write out
    open( OUT, ">$outheteroplasmy" ) or die $!;
    print OUT join "\t",
      (
        "#chr",               "pos",
        "ref",                "forward_A",
        "forward_T",          "forward_C",
        "forward_G",          "reverse_A",
        "reverse_T",          "reverse_C",
        "reverse_G",          "heteroplasmy",
        "95\%ci_lower",       "95\%ci_upper",
        "major_allele",       "minor_allele",
        "major_allele_count", "minor_allele_count",
        "gene","genedetail","exonic_function","aminochange",
        "strand_bias","pathogenic_variants","diseases","links",
        "fisher.pvalue","fisher.adjust.pvalue","fisher.phred.score",
        "empirical.probability","empirical.phred.score"
      );
    
    print OUT "\n";

    foreach my $loc ( sort { $a <=> $b } keys %result ) {
        if ( $result{$loc}->{'stat'} == 2 ) {
            my $convertloc=$loc;
            my $convertref=$result{$loc}->{'ref'};
            if($inref ne $outref){ #need to convert genome location
                  if($inref eq 'hg19' && $outref eq 'rCRS'){
                        $convertloc=$convert->hg19TorCRS($loc);
                        $convertref=$convert->rCRSref($convertloc);
                  }else{
                        $convertloc=$convert->rCRSTohg19($loc);
                        $convertref=$convert->hg19ref($convertloc);
                  }
                  next if (!defined($convertloc));  #can't be mapped between different assembly
            }
            print OUT join "\t",
              (
                #$result{$loc}->{'chr'},
                "MT",
                $convertloc,
                #$result{$loc}->{'ref'},
                $convertref,
                $result{$loc}->{'A'},
                $result{$loc}->{'T'},
                $result{$loc}->{'C'},
                $result{$loc}->{'G'},
                $result{$loc}->{'a'},
                $result{$loc}->{'t'},
                $result{$loc}->{'c'},
                $result{$loc}->{'g'},
                $result{$loc}->{'heteroplasmy'},
                $result{$loc}->{'lower'},
                $result{$loc}->{'upper'},
                $result{$loc}->{'major_allele'},
                $result{$loc}->{'minor_allele'},
                $result{$loc}->{'major_allele_count'},
                $result{$loc}->{'minor_allele_count'}
              );
              print OUT "\t";
              
              my $alt_allele=$result{$loc}->{'major_allele'};
              $alt_allele=$result{$loc}->{'minor_allele'} if ($alt_allele eq $convertref); 
              
              #need to determine whether rCRS  or hg19
              print OUT join "\t",($mitoanno->annotation($convertloc,$convertref,$alt_allele));
            #if ( $sb > 0 ) {
                print OUT "\t", $result{$loc}->{'sb'};
            #}
             print OUT "\t";
             print OUT join "\t",($mitoanno->pathogenic($convertloc));

             # Added on 2013-02-12
             print OUT "\t";
             my $empirical=empirical_bayesian($result{$loc}->{'major_allele_count'},$result{$loc}->{'minor_allele_count'},$hp/100,1000,$A,$B);
             print OUT join "\t",($rawpvalue{$loc},$adjustpref->{$loc},phred($rawpvalue{$loc}),$empirical,phred(1-$empirical));
            
            print OUT "\n";
        }
    }
    return %result;
}

#Calculate 95% confidence interval
#Parameters:
#$p  rate range from 0-1
#$n  total number of allele
sub _ci {
    my ( $p, $n ) = @_;
    my $tmp   = 1.96 * sqrt( ( $p * ( 1 - $p ) ) / $n );
    my $lower = 0;
    my $upper = 1;
    $lower = $p - $tmp if ( ( $p - $tmp ) > 0 );
    $upper = $p + $tmp if ( ( $p + $tmp ) < 1 );
    #Format the numeric
    $lower=_formatnumeric($lower);
    $upper=_formatnumeric($upper);
    return ( $lower, $upper );
}

sub _formatnumeric{
    my ($i) = @_;
    my $r=sprintf("%.3f",$i);
    if($r eq '0.000' && $i!=0){#use scientific notation
        $r=sprintf("%.3e",$i);
    }
    return $r;
}

#Give a bam/sam file($inbam) and a mithochondrial bed file($mbed),
#fetch the alignment reads in mithochondrial genome
#and output in bam format
#Parameters:
#$inbam  input bam or sam format file
#$isbam  0/1, 0 indicates it's in sam format while 1 indicates a bam format
#$mbed   bed file of mithochondrial genome
#$mmq    minimum map quality
sub _get_mitochondrial_bam {
    my ( $inbam, $isbam, $mbed, $mmq, $outbam ) = @_;
    my $comm = "";
    if ($isbam) {
        $comm = "$samtools view -b -q $mmq -L $mbed $inbam > $outbam";
    }
    else {
        $comm = "$samtools view -bS -q $mmq -L $mbed $inbam > $outbam";
    }
    _run($comm);
}

#Give a bam/sam file($inbam) and a mithochondrial bed file($mbed),
#fetch the alignment reads in mithochondrial genome
#and output in bam format
#Parameters:
#$inbam  input bam or sam format file
#$isbam  0/1, 0 indicates it's in sam format while 1 indicates a bam format
#$mbed   bed file of mithochondrial genome
#$mmq    minimum map quality
sub _get_mitochondrial_bam_advance {
    my ( $inbam, $isbam, $mbed, $mmq, $outbam ) = @_;
    my $comm = "perl $mitomap -r $bwaindex -i $inbam -b $mbed -mmq $mmq -o $outbam -samtools $samtools -bwa $bwa";
    _run($comm);
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
sub _determine_somatic {
    my ( $tumorbase, $normalbase, $sp, $sa, $isall, $somaticoutput ) = @_;
    my %tumor  = _read_basecall($tumorbase);
    my %normal = _read_basecall($normalbase);

    my %result;
    open( OUT, ">$somaticoutput" ) or die $!;
    print OUT join "\t",
      (
        "#chr",       "pos",            "ref", "tumorGenotype",
        "tumorDepth", "normalGenotype", "normalDepth",
        "Mito.gene","Mito.genedetail\n"
      );

#Get the overlapped locations and then loop them to determine whether it is a somatic mutation
    my @commloc =
      grep { exists( $tumor{$_} ) }
      sort { $a <=> $b } keys %normal;    #numeric sort
    foreach my $loc (@commloc) {

        #Only if normal is not heteroplasmy
        my $normalGeno = "";
        my $normalDP   = "";
        my $tumorGeno  = "";
        my $tumorDP    = "";
        my $issomatic  = 0;

        if ( $normal{$loc}->{'minor_allele_count'} == 0 ) { #though the minor allele in normal is 0, but the major allele in normal and turmor may different
            if ( $tumor{$loc}->{'minor_allele_count'} == 0 ) {
#If at this site, both normal and tumor are not heteroplasmy, then determine whether their
#major_allele is the same
                if ( $normal{$loc}->{'major_allele'} ne
                    $tumor{$loc}->{'major_allele'} )
                {
                    #a somatic mutation
                    $issomatic = 1;
                }
            }else {
                if ( $normal{$loc}->{'major_allele'} ne
                    $tumor{$loc}->{'major_allele'} )
                {
#No matter the minor_allele_count meet the cutoff of not, it is already a somatic mutation
                    $issomatic = 1;
                }
                else {

       #Only the minor_allele_count meet the cutoff and it is a somatic mutation
                    my $ratio = 0;
                    if ($isall) {
                        $ratio =
                          $tumor{$loc}->{'minor_allele_count'} /
                          ( sum @{ $tumor{$loc}->{'alleles_count'} } );
                    }
                    else {
                        $ratio =
                          $tumor{$loc}->{'minor_allele_count'} /
                          ( $tumor{$loc}->{'major_allele_count'} +
                              $tumor{$loc}->{'minor_allele_count'} );
                    }

                    if (   $ratio > $sp / 100
                        && $tumor{$loc}->{'minor_allele_count'} > $sa )
                    {
                        $issomatic = 1;
                    }
                }
            }
        }

        #Assign this loc to %result
        if ($issomatic) {
            $normalGeno = $normal{$loc}->{'major_allele'};
            $normalDP   = $normal{$loc}->{'major_allele_count'};

            my @tmp = grep { $_ != 0 } @{ $tumor{$loc}->{'alleles_count'} };
            $tumorDP   = join "|", @tmp;
            $tumorGeno = join "|", @{ $tumor{$loc}->{'alleles'} }[ 0 .. $#tmp ];
            $result{$loc}->{'ref'}        = $normal{$loc}->{'reference_allele'};
            $result{$loc}->{'normalGeno'} = $normalGeno;
            $result{$loc}->{'normalDP'}   = $normalDP;
            $result{$loc}->{'tumorGeno'}  = $tumorGeno;
            $result{$loc}->{'tumorDP'}    = $tumorDP;
            
            my $convertloc=$loc;
            my $convertref=$normal{$loc}->{'reference_allele'};
            #my $convertref=$result{$loc}->{'ref'};
            if($inref ne $outref){ #need to convert genome location
                  if($inref eq 'hg19' && $outref eq 'rCRS'){
                        $convertloc=$convert->hg19TorCRS($loc);
                        $convertref=$convert->rCRSref($convertloc);
                  }else{
                        $convertloc=$convert->rCRSTohg19($loc);
                        $convertref=$convert->hg19ref($convertloc);
                  }
                  next if (!defined($convertloc));  #can't be mapped between different assembly
            }
                        
            print OUT join "\t",
              (
                "MT", $convertloc, $convertref,
                $tumorGeno, $tumorDP, $normalGeno, $normalDP
              );
              my ($genetmp,$genedetailtmp,undef,undef) = $mitoanno->annotation($convertloc);
            print OUT "\t",$genetmp,"\t",$genedetailtmp;
            print OUT "\n";
        }
    }
    close OUT;
    return %result;
}

#Given a basecall format file, parse it and return a hash table
#Parameters:
#$in  input basecall file
sub _read_basecall {
    my ($in) = @_;
    my %return;
    open( IN, $in ) or die $!;
    <IN>;
    while (<IN>) {
        s/\r|\n//g;
        my (
            $chr,       $loc,       $ref,       $forward_A,
            $forward_T, $forward_C, $forward_G, $reverse_A,
            $reverse_T, $reverse_C, $reverse_G
        ) = split "\t";

        # Get the major allele and minor allele.
        my %atcg;
        $atcg{A} = $forward_A + $reverse_A;
        $atcg{T} = $forward_T + $reverse_T;
        $atcg{C} = $forward_C + $reverse_C;
        $atcg{G} = $forward_G + $reverse_G;
        my @atcg         = sort { $atcg{$b} <=> $atcg{$a} } keys %atcg;
        my $major_allele = $atcg[0];
        my $minor_allele = $atcg[1];

        $return{$loc}->{'major_allele'}           = $atcg[0];
        $return{$loc}->{'major_allele_count'}     = $atcg{ $atcg[0] };
        $return{$loc}->{'minor_allele'}           = $atcg[1];
        $return{$loc}->{'minor_allele_count'}     = $atcg{ $atcg[1] };
        $return{$loc}->{'reference_allele'}       = $ref;
        $return{$loc}->{'reference_allele_count'} = $atcg{$ref};
        $return{$loc}->{'alleles'}                = [@atcg];
        $return{$loc}->{'alleles_count'} =
          [ @atcg{@atcg} ];    #Another ways to fetch value from a hash
        $return{$loc}->{'basecallline'} = [
            $chr,       $loc,       $ref,       $forward_A,
            $forward_T, $forward_C, $forward_G, $reverse_A,
            $reverse_T, $reverse_C, $reverse_G
        ];
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
sub _pileup {
    my ( $inbam, $isbam, $mbq, $regionbed, $refseq, $outpileup ) = @_;
    my $comm = "";
    if ($isbam) {
        $comm =
"$samtools mpileup -l $regionbed -Q $mbq -f $refseq $inbam > $outpileup 2>/dev/null";
    }
    else {
        $comm =
"$samtools view -Su $inbam | $samtools mpileup -l $regionbed -Q $mbq -f $refseq - > $outpileup 2>/dev/null";
    }
    _run($comm);
}

#Given a pileup file, parse the pileup file into basecall format and write it out
#Parameters:
#$inpileup  one-sample pileup format file
#$mbq       minimum base quality
#$offset    quality encoding offset, Default is 33
#$outbase   output allele count file parsed from pileup file

sub _parse_pileup {
    my ( $inpileup, $mbq, $offset, $outbase ) = @_;
    open( IN,  $inpileup )   or die $!;
    open( OUT, ">$outbase" ) or die $!;

    #Uppeer case is from forward strand, while lower case is from reverse strand
    print OUT "Chr\t" . "Loc\t" . "Ref\t" . "A\t" . "T\t" . "C\t" . "G\t"
      . "a\t" . "t\t" . "c\t" . "g\n";
    while (<IN>) {
        s/\r|\n//g;
        my ( $chr, $loc, $ref, $dp, $bases, $bq ) = split /\s+/;
        $ref = uc($ref);

#do some modificaton on $bases to remove additional characters
#1,remove the ^. pattern (marks the start of a read segment and the ASCII of the character following `^' minus 33 gives the mapping quality)
        $bases =~ s/\^.//g;

        #2,remove the $ pattern (marks the end of a read segment)
        $bases =~ s/\$//g;

#3,remove -[0-9]+[ACGTNacgtn]+ pattern (denotes a deletion of one or more bases)
        my %hash = ();
        while ( $bases =~ /-(\d+)/g ) {
            $hash{$1} = 1;
        }
        foreach my $k ( keys %hash ) {
            $bases =~ s/-$k[ACGTNacgtn]{$k}//g;
        }

#4,remove +[0-9]+[ACGTNacgtn]+ pattern (denotes a insertion of one or more bases)
        %hash = ();
        while ( $bases =~ /\+(\d+)/g ) {
            $hash{$1} = 1;
        }
        foreach my $k ( keys %hash ) {
            $bases =~ s/\+$k[ACGTNacgtn]{$k}//g;
        }

#Now @base and @bq have the same length (Note that: the < or > in the $bases denote a gap)
        my @base = split( //, $bases );
        my @bq   = split( //, $bq );
        my $A    = 0;
        my $T    = 0;
        my $C    = 0;
        my $G    = 0;
        my $a    = 0;
        my $t    = 0;
        my $c    = 0;
        my $g    = 0;

        #start the loop
        for ( my $i = 0 ; $i < @base ; $i++ ) {
            my $ch    = $base[$i];
            my $score = ord( $bq[$i] ) - $offset;  #Need to be more robust later
            if ( $score >= $mbq ) {
                if ( $ch eq "A" ) {
                    $A++;
                }
                elsif ( $ch eq "T" ) {
                    $T++;
                }
                elsif ( $ch eq "C" ) {
                    $C++;
                }
                elsif ( $ch eq "G" ) {
                    $G++;
                }
                elsif ( $ch eq "a" ) {
                    $a++;
                }
                elsif ( $ch eq "t" ) {
                    $t++;
                }
                elsif ( $ch eq "c" ) {
                    $c++;
                }
                elsif ( $ch eq "g" ) {
                    $g++;
                }
                elsif ( $ch eq "." ) {
                    if ( $ref eq "A" ) {
                        $A++;
                    }
                    elsif ( $ref eq "T" ) {
                        $T++;
                    }
                    elsif ( $ref eq "C" ) {
                        $C++;
                    }
                    elsif ( $ref eq "G" ) {
                        $G++;
                    }
                }
                elsif ( $ch eq "," ) {
                    if ( $ref eq "A" ) {
                        $a++;
                    }
                    elsif ( $ref eq "T" ) {
                        $t++;
                    }
                    elsif ( $ref eq "C" ) {
                        $c++;
                    }
                    elsif ( $ref eq "G" ) {
                        $g++;
                    }
                }
            }    #end the condition  $score>=$mbq
        }    #end the loop

        if ( $A + $T + $C + $G + $a + $t + $c + $g > 0 ) {
            print OUT "$chr\t$loc" . "\t" 
              . $ref . "\t" 
              . $A . "\t" 
              . $T . "\t"
              . $C . "\t"
              . $G . "\t"
              . $a . "\t"
              . $t . "\t"
              . $c . "\t"
              . $g . "\n";
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
sub _sb {
    my ( $a, $b, $c, $d ) = @_;
    if ( $a + $b == 0 || $c + $d == 0 || $b + $d == 0 ) {
        return 0;
    }
    else {
        return
          _formatnumeric(abs( $b / ( $a + $b ) - $d / ( $c + $d ) ) /
          ( ( $b + $d ) / ( $a + $b + $c + $d ) ));
    }
}

#Check whether samtools exist in ymy $PATH, If not, the program will exit
sub _check_samtools {
    my $r = `which samtools`;
    if ($r) {
        chomp($r);
        _info("Checking: samtools in $r");
    }
    else {
        _info("Error: samtools not installed yet");
        exit(1);
    }
}

#Customed print with time included, if providing $flag then there will be no "\n"
sub _info {
    my ( $s, $flag ) = @_;
    print "[", scalar(localtime), "] $s";
    print "\n" unless ($flag);
}

sub _error {
    my ( $s, $flag ) = @_;
    print STDERR"[", scalar(localtime), "] [ERROR] $s";
    print STDERR "\n" unless ($flag);
}

sub _warn {
    my ( $s, $flag ) = @_;
    print STDERR"[", scalar(localtime), "] [WARN] $s";
    print STDERR "\n" unless ($flag);
}

sub _run {
    my ($comm) = @_;
    if($comm !~/>/){  #suppress the output
        $comm.=">/dev/null";
    }
    system($comm) == 0 or die $!;
}




##Generate the report
sub _make_report {
    my $createdtime=localtime;
    open( OUT, ">$mitoreport" ) or die $!;
    print OUT <<HTML;
<html>
<head><title>MitoSeek Report</title>
<style type="text/css">
   \@media screen {
  div.summary {
    width: 18em;
    position:fixed;
    top: 3em;
    margin:1em 0 0 1em;
  }
  
  div.main {
    display:block;
    position:absolute;
    overflow:auto;
    height:auto;
    width:auto;
    top:4.5em;
    bottom:2.3em;
    left:18em;
    right:0;
    border-left: 1px solid #CCC;
    padding:0 0 0 1em;
    background-color: white;
    z-index:1;
  }
  
  div.header {
    background-color: #EEE;
    border:0;
    margin:0;
    padding: 0.5em;
    font-size: 200%;
    font-weight: bold;
    position:fixed;
    width:100%;
    top:0;
    left:0;
    z-index:2;
  }

  div.footer {
    background-color: #EEE;
    border:0;
    margin:0;
    padding:0.5em;
    height: 1.3em;
    overflow:hidden;
    font-size: 100%;
    font-weight: bold;
    position:fixed;
    bottom:0;
    width:100%;
    z-index:2;
  }
  
  img.indented {
    margin-left: 3em;
  }
 }
 
 \@media print {
    img {
        max-width:100% !important;
        page-break-inside: avoid;
    }
    h2, h3 {
        page-break-after: avoid;
    }
    div.header {
      background-color: #FFF;
    }
    
 }
 
 body {    
  font-family: sans-serif;   
  color: #000;   
  background-color: #FFF;
  border: 0;
  margin: 0;
  padding: 0;
  }
  
  div.header {
  border:0;
  margin:0;
  padding: 0.5em;
  font-size: 200%;
  font-weight: bold;
  width:100%;
  }    
  
  #header_title {
  display:inline-block;
  float:left;
  clear:left;
  }
  #header_filename {
  display:inline-block;
  float:right;
  clear:right;
  font-size: 12px;
  margin-right:2em;
  text-align: right;
  margin-top:3em;
  }

  div.header h3 {
  font-size: 50%;
  margin-bottom: 0;
  }
  
  div.summary ul {
  padding-left:0;
  list-style-type:none;
  font-weight:bold;
  }
  
  div.summary ul li img {
  margin-bottom:-0.5em;
  margin-top:0.5em;
  }
      
  div.main {
  background-color: white;
  }
      
  div.module {
  padding-bottom:1.5em;
  padding-top:1.5em;
  }
      
  div.footer {
  background-color: #EEE;
  border:0;
  margin:0;
  padding: 0.5em;
  font-size: 100%;
  font-weight: bold;
  width:100%;
  }


  a {
  color: #000080;
  }

  a:hover {
  color: #800000;
  }
      
  h2 {
  color: #800000;
  padding-bottom: 0;
  margin-bottom: 0;
  clear:left;
  }

 
  img {
  padding-top: 0;
  margin-top: 0;
  border-top: 0;
  }

  
  p {
  padding-top: 0;
  margin-top: 0;
  }
  
  #image-container { 
min-width:1600px; /* width of 2 images (1600px) plus images' padding and border (0px) plus images' margin (0px) */ 
height:600px; 
padding:0px 0px; 
margin:0; 
border:0px solid #ffffff; 
background-color:#ffffff; 
list-style-type:none; 
} 
  
#image-container li { 
width:50%; 
float:left; 
text-align: center;
margin-bottom:10px
} 

#image-container img { 
display:block; 
width:800px; 
height:600px; 
padding:0px; 
border:1px solid #d3d2d2; 
margin:auto; 
background-color:#fff; 
} 

#image-container1 { 
min-width:800px; /* width of 5 images (625px) plus images' padding and border (60px) plus images' margin (50px) */ 
height:600px; 
padding:0px 0px; 
margin:0; 
border:0px solid #ffffff; 
background-color:#ffffff; 
list-style-type:none; 
} 
  
#image-container1 li { 
width:100%; 
float:left; 
text-align: center;
margin-bottom:10px
} 

#image-container1 img { 
display:block; 
width:800px; 
height:600px; 
padding:5px; 
border:1px solid #d3d2d2; 
#margin:auto; 
background-color:#fff; 
}


#image-container2 { /*for cirsocs*/ 
min-width:1600px; /* width of 5 images (625px) plus images' padding and border (60px) plus images' margin (50px) */ 
height:800px; 
padding:0px 0px; 
margin:0; 
border:0px solid #ffffff; 
background-color:#ffffff; 
list-style-type:none; 
} 
  
#image-container2 li { 
width:50%; 
float:left; 
text-align: center;
margin-bottom:10px
} 

#image-container2 img { 
display:block; 
width:800px; 
height:800px; 
padding:5px; 
border:1px solid #d3d2d2; 
#margin:auto; 
background-color:#fff; 
}

#image-container3 { /*for cirsocs*/ 
min-width:1000px; /* width of 5 images (625px) plus images' padding and border (60px) plus images' margin (50px) */ 
height:1000px; 
padding:0px 0px; 
margin:0; 
border:0px solid #ffffff; 
background-color:#ffffff; 
list-style-type:none; 
} 
  
#image-container3 li { 
width:100%; 
float:left; 
text-align: center;
margin-bottom:10px
} 

#image-container3 img { 
display:block; 
width:1000px; 
height:1000px; 
padding:5px; 
border:1px solid #d3d2d2; 
#margin:auto; 
background-color:#fff; 
}



table {
    margin-left: 3em;
    text-align: center;
    border-width: 1px;
    border-spacing: 2px;
    border-style: outset;
    border-color: gray;
    border-collapse: collapse;
    background-color: rgb(209, 255, 255);
}
table th {
    border-width: 2px;
    padding: 2px;
    border-style: solid;
    border-color: gray;
    background-color: rgb(209, 255, 255);
    text-align: center;
    padding: 0.4em;
}
table td {
    border-width: 2px;
    padding: 2px;
    border-style: solid;
    border-color: gray;
    background-color: white;
    font-family: monospace; 
  text-align: center;
  color: #000;
  padding: 0.4em;
}
pre {
    padding: 1em;
    border: 1px dashed #2f6fab;
    color: black;
    background-color: #f9f9f9;
    line-height: 1.1em;
}
</style>
</head>
<body>
<div class="header">
    <div id="header_title">MitoSeek Report</div>
    <div id="header_filename">Created Time: $createdtime</div>
</div>

<div class="summary">
    <h2>Table of Contents</h2>
        <ul>
            <li><a href="#M001">Command</a></li>
            <li> <a href="#M0">QC</a></li>
                <ol>
                    <li> <a href="#M000">Percent of base pairs covered</a></li>
                    <li> <a href="#M00">Per base quality</a></li>
                    <li> <a href="#M01">Mapping quality</a></li>
                    <li> <a href="#M02">Depth distribution</a></li>
                    <li> <a href="#M03">Template length distribution</a></li>
                </ol>
            <li> <a href="#M1">Heteroplasmy</a></li>
            <li> <a href="#M4">Structural Changes</a></li>
            <li> <a href="#M2">Somatic Mutation</a></li>
            <li> <a href="#M3">Relative Copy Number Estimation</a></li>
            <li> <a href="#M5">File List</a></li>
        </ul>
</div>
<div class="main">
    <div class="module"><h2 id="M001">Command</h2>
        <p>Your command for generating this report, keep it in order to reproduce the result.</p>
        <pre>$commandline</pre>
    </div>
    <div class="module"><h2 id="M0">QC</h2>
        <p>
        Quality control (QC) is applied to the mapped reads on mitochondria. If -L [bed] is provided, the QC report following is constrained to your provided region. The QC report contains the following 4 parts.
        </p>
        <h3 id="M000">Percent of base pairs covered</h3>
HTML
        
       if($inbam2){
           print OUT "<h3>Tumor</h3>",  _generate_html_table($percentofbasepairscovered_table1);
           print OUT "<h3>Normal</h3>", _generate_html_table($percentofbasepairscovered_table2);
        }else{
            print OUT _generate_html_table($percentofbasepairscovered_table1);
        }
       print OUT <<HTML;

        <h3 id="M00">Per base quality</h3>
        <p>
        The y-axis on the graph shows the quality scores and the x-axis on the graph shows the positions in the fastq file. For each position a BoxWhisker type plot is drawn. The elements of the plot are as follows:

       <ul>
        <li>The central blue line is the median value</li>
        <li>The green box represents the inter-quartile range (25-75%)</li>
        <li>The upper and lower whiskers represent the 10% and 90% points</li>
        <li>The blue '+' mark represents the mean quality</li>
        <li>The blue '*' mark outside the upper and lower whiskers represents the outlier points</li>
        </ul>

        </p>
HTML

if($inbam2){
    print OUT <<HTML;
    <ul id='image-container'>
        <li><img src='$perbasequality_figure1' alt='Per base quality figure'>Tumor</li>
        <li><img src='$perbasequality_figure2' alt='Per base quality figure'>Normal</li>
    </ul>
HTML
}else{
    print OUT <<HTML;
    <ul id='image-container1'>
        <li><img src='$perbasequality_figure1' alt='Per base quality figure'></li>
    </ul>
HTML
}

print OUT <<HTML;
        <h3 id="M01">Mapping quality</h3>
        <p>
        The y-axis on the graph shows the density and the x-axis on the graph shows the mapping quality scores. The mapping quality score is stored in the 5th column of a SAM/BAM file. It equals to <b>-10log10Pr{mapping position is wrong}</b>, rounded to the nearest integer. A summary of the mapping score is also included in the middle of the plot.
        </p>
HTML
if($inbam2){
    print OUT <<HTML;
    <ul id='image-container'>
        <li><img src='$mappingquality_figure1' alt='Mapping quality figure'>Tumor</li>
        <li><img src='$mappingquality_figure2' alt='Mapping quality figure'>Normal</li>
    </ul>
HTML
}else{
    print OUT <<HTML;
    <ul id='image-container1'>
        <li><img src='$mappingquality_figure1' alt='Mapping quality figure'></li>
    </ul>
HTML
}

print OUT <<HTML;
        <h3 id="M02">Depth distribution</h3>
        <p>
        The y-axis on the graph shows the density and the x-axis on the graph shows the depth. A summary of the depth is also included in the middle of the plot. 
        </p>
HTML
if($inbam2){
    print OUT <<HTML;
    <ul id='image-container'>
        <li><img src='$depthdistribution_figure1' alt='Depth distribution figure'>Tumor</li>
        <li><img src='$depthdistribution_figure2' alt='Depth distribution figure'>Normal</li>
    </ul>
HTML
}else{
    print OUT <<HTML;
    <ul id='image-container1'>
        <li><img src='$depthdistribution_figure1' alt='Depth distribution figure'></li>
    </ul>
HTML
}
print OUT <<HTML;
        <h3 id="M03">Template length distribution</h3>
        <p>
        The y-axis on the graph shows the density and the x-axis on the graph shows the template length. The template length is stored in the 9th column of a SAM/BAM file. It is set as 0 for single-segment template or when the information is unavailable.
        A summary of the template depth is also included in the middle of the plot. 
        </p>
HTML
if(-e $templatelengthdistribution_figure1){
    if($inbam2){
        print OUT <<HTML;
        <ul id='image-container'>
            <li><img src='$templatelengthdistribution_figure1' alt='Template length distribution figure'>Tumor</li>
            <li><img src='$templatelengthdistribution_figure2' alt='Template length distribution figure'>Normal</li>
        </ul>
HTML
}else{
        print OUT <<HTML;
        <ul id='image-container1'>
            <li><img src='$templatelengthdistribution_figure1' alt='Template length distribution figure'>Tumor</li>
        </ul>
HTML
    }
}else{
    print OUT "<b>Not Avaiable</b><br/>";
}
print OUT <<HTML;
    </div>
    <div class="module"><h2 id="M1">Heteroplasmy</h2>
        <p>Heteroplasmy detection threshold is defined on two scales: read count (-ha) and read percentage (-hp). Read count denotes the number of reads we must observe to support heteroplasmy while read percentage denotes the percentage of reads we must observe to support heteroplasmy. Both scales can be used together or individually. The minimum recommended depth requirement (-d) for detecting heteroplasmy is 50. Lower depth will severely damage the confidence of heteroplasmy calling.</p>
HTML

    if ($inbam2) {
        if($producecircosplot){
            print OUT <<HTML;
        <ul id='image-container2'>
            <li><img src='circos/$mitocircosheteroplasmyfigure1' alt='Circos plot of heteroplasmy on tumor'>Tumor</li>
            <li><img src='circos/$mitocircosheteroplasmyfigure2' alt='Circos plot of heteroplasmy on normal'>Normal</li>
        </ul>
HTML
        }
        print OUT "<h3>Tumor</h3>",  _generate_html_table($mitoheteroplasmy1);
        print OUT "<h3>Normal</h3>", _generate_html_table($mitoheteroplasmy2);
    }
    else {
    if($producecircosplot){
            print OUT <<HTML;
        <ul id='image-container3'>
            <li><img src='circos/$mitocircosheteroplasmyfigure1' alt='Circos plot of heteroplasmy on tumor'></li>
        </ul>
HTML
        }
        print OUT _generate_html_table($mitoheteroplasmy1);
    }


    print OUT <<HTML;
    </div>
    <div class="module"><h2 id="M4">Structural Changes</h2>
    <p>
    MitoSeek reports mitochondria structural changes when pair-end sequencing data is given as input. During alignment, a portion of the read-pairs will be discordantly mapped, meaning one read of the pair is aligned to mitochondria and the mate pair is aligned elsewhere. Such reads are like to be the results of alignment errors from homologous regions between mitochondria genome and other genomes. However, they could also indicate mitochondria integration into other genomes which has been reported to be possible by multiple studies. Only those structural changes with >=n (-str) spanning reads support will be outputed.
    </p>
HTML
    if ($inbam2) {
        print OUT "<h3>Tumor</h3>",  _generate_html_table($mitostructure1);
        print OUT "<h3>Normal</h3>", _generate_html_table($mitostructure2);
    }
    else {
        print OUT  _generate_html_table($mitostructure1);
    }

print OUT <<HTML;
    </div>
    <div class="module"><h2 id="M2">Somatic Mutation</h2>
    <p>
    MitoSeek takes the input bam provided by <b>-i</b> as tumor while the other input bam by <b>-j</b> as its control normal. We propose to compare the empirical allele counts between tumor and normal control directly instead of using a genotype caller. MitoSeek can extract empirical allele count for every mitochondria position then compare the allele counts between tumor and normal to determine somatic mutation status. Two parameters (-sp and -sa) control the somatic mutation detection.
    </p>
HTML
    if ($inbam2) {
        if($cs){
            print OUT <<HTML;
        <ul id='image-container3'>
            <li><img src='circos/$mitocircossomaticfigure' alt='Circos plot of somatic mutations'></li>
        </ul>
HTML
        }
        print OUT _generate_html_table($mitosomatic);
    }
    else {
        print OUT "NA\n";
    }

    print OUT <<HTML;
    </div>
    <div class="module"><h2 id="M3">Relative Copy Number Estimation</h2>
    <p>
    Two methods are implemented in MitoSeek to estimate the relative copy number of mithochondria, namely '<b>byRead</b>' and '<b>byDepth</b>'.<br/>
    <h4>byRead</h4>
    <b>CN=Rm/Rt</b>,  where <b>Rm</b> is the reads aligned to mitochondria and passed quality filter and <b>Rt</b> is the total reads passed quality filter. 
    <h4>byDepth</h4>
    <b>CN=Dm/Dt</b>, where <b>Dm</b> is the average depth of mitochondria, and <b>Dt</b> is the average of exome or whole genome, depending on your input data.
    </p>
HTML

    if ($inbam2) {
        if($cn && $type !=4){
            print OUT "<h3>Tumor</h3>",  _generate_html_table($mitocnv1);
            print OUT "<h3>Normal</h3>", _generate_html_table($mitocnv2);
        }
    }
    else {
        if($cn && $type !=4){
            print OUT  _generate_html_table($mitocnv1);
        }else{
            print OUT "NA\n";
        }
    }


    print OUT <<HTML;
    </div>
    <div class="module"><h2 id="M5">File list</h2>
    <p>File list description</p>
HTML
    
    print OUT _print_file_list()," </div>\n";

    print OUT <<HTML;
    </div>
    </div>
<div class="footer"  style="background-color: #EEE;font-weight:" align='center'>Produced by <a href="https://github.com/riverlee/MitoSeek">MitoSeek</a> (Jiang Li)</div>
</body>
</html>
HTML

    close OUT;
}


#Print out file list into the report
sub _print_file_list(){
    my $html="<table border='1'>\n<tr>".
             "<th>Category</th>".
             "<th>File Name</th>".
             "<th>Description</th></tr>\n";
    if($inbam2){
        #QC related table
        if($qc){
            $html.="<tr>".
                    "<td rowspan='18'>QC</td>\n".
                    "<td>"._html_link($percentofbasepairscovered_table1)."</td>\n".
                    "<td>Percent of base pairs covered (Tumor)</td>\n".
               "</tr>\n".
               "<tr>".
                    "<td>"._html_link($percentofbasepairscovered_table2)."</td>\n".
                    "<td>Percent of base pairs covered (Normal)</td>\n".
               "</tr>\n".
               "<tr>".
                    "<td>"._html_link($perbasequality_figure1)."</td>\n".
                    "<td> Per base quality plot (Tumor)</td>\n".
               "</tr>\n".
               "<tr>".
                    "<td>"._html_link($perbasequality_table1)."</td>\n".
                    "<td> Data for per base quality plot (Tumor)</td>\n".
                "</tr>".
                "<tr>".
                    "<td>"._html_link($perbasequality_figure2)."</td>\n".
                    "<td> Per base quality plot (Normal)</td>\n".
               "</tr>\n".
               "<tr>".
                    "<td>"._html_link($perbasequality_table2)."</td>\n".
                    "<td> Data for per base quality plot (Normal)</td>\n".
                "</tr>".
                
                "<tr>".
                    "<td>"._html_link($mappingquality_figure1)."</td>\n".
                    "<td> Mapping quality plot (Tumor)</td>\n".
                "</tr>".
                "<tr>".
                    "<td>"._html_link($mappingquality_table1)."</td>\n".
                    "<td> Data for mapping quality plot (Tumor)</td>\n".
                "</tr>".
                 "<tr>".
                    "<td>"._html_link($mappingquality_figure2)."</td>\n".
                    "<td> Mapping quality plot (Normal)</td>\n".
                "</tr>".
                "<tr>".
                    "<td>"._html_link($mappingquality_table2)."</td>\n".
                    "<td> Data for mapping quality plot (Normal)</td>\n".
                "</tr>".
                
                "<tr>".
                    "<td>"._html_link($depthdistribution_figure1)."</td>\n".
                    "<td> Depth distribution plot (Tumor)</td>\n".
                "</tr>".
                "<tr>".
                    "<td>"._html_link($depthdistribution_table1)."</td>\n".
                    "<td> Data for depth distribution plot (Tumor)</td>\n".
                "</tr>".
                "<tr>".
                    "<td>"._html_link($depthdistribution_figure2)."</td>\n".
                    "<td> Depth distribution plot (Normal)</td>\n".
                "</tr>".
                "<tr>".
                    "<td>"._html_link($depthdistribution_table2)."</td>\n".
                    "<td> Data for depth distribution plot (Normal)</td>\n".
                "</tr>".
                
                "<tr>".
                    "<td>"._html_link($templatelengthdistribution_figure1)."</td>\n".
                    "<td> Template length distribution plot (Tumor)</td>\n".
                "</tr>".
                "<tr>".
                    "<td>"._html_link($templatelengthdistribution_table1)."</td>\n".
                    "<td> Data for template length distribution plot (Tumor)</td>\n".
                "</tr>".
                "<tr>".
                    "<td>"._html_link($templatelengthdistribution_figure2)."</td>\n".
                    "<td> Template length distribution plot (Normal)</td>\n".
                "</tr>".
                "<tr>".
                    "<td>"._html_link($templatelengthdistribution_table2)."</td>\n".
                    "<td> Data for template length distribution plot (Normal)</td>\n".
                "</tr>";
          }
         #heteroplasmy
        if($producecircosplot){
             $html.="<tr>".
                    "<td rowspan='10'>Heteroplasmy</td>\n".
                    "<td>"._html_link($mitoheteroplasmy1)."</td>".
                    "<td> Heteroplasmy result (Tumor) </td>".
               "</tr>".
               "<tr>".
                    "<td>"._html_link($mitoheteroplasmy2)."</td>".
                    "<td> Heteroplasmy result (Normal) </td>".
               "</tr>".
              
               "<tr>".
                    "<td>"._html_link("circos/".$mitocircosheteroplasmyfigure1)."</td>".
                    "<td> Circos plot of heteroplasmy result (Tumor)</td>".
               "</tr>".
               "<tr>".
                    "<td>"._html_link("circos/".$mitocircosheteroplasmyfigure2)."</td>".
                    "<td> Circos plot of heteroplasmy result (Normal)</td>".
               "</tr>".
               
               "<tr>".
                    "<td>"._html_link("circos/".$mitocircosheteroplasmyconfig1)."</td>".
                    "<td> Configure file for circos plot of heteroplasmy result (Tumor)</td>".
               "</tr>".
               "<tr>".
                    "<td>"._html_link("circos/".$mitocircosheteroplasmyconfig2)."</td>".
                    "<td> Configure file for circos plot of heteroplasmy result (Normal)</td>".
               "</tr>".
               
               "<tr>".
                    "<td>"._html_link("circos/".$mitoheteroplasmytextoutput1)."</td>".
                    "<td> Data file (text labels) for circos plot of heteroplasmy result (Tumor)</td>".
               "</tr>".
               "<tr>".
                    "<td>"._html_link("circos/".$mitoheteroplasmytextoutput2)."</td>".
                    "<td> Data file (text labels) for circos plot of heteroplasmy result (Normal)</td>".
               "</tr>".
               
               "<tr>".
                    "<td>"._html_link("circos/".$mitoheteroplasmyscatteroutput1)."</td>".
                    "<td> Data file (scatter plots) for circos plot of heteroplasmy result (Tumor)</td>".
               "</tr>".
                "<tr>".
                    "<td>"._html_link("circos/".$mitoheteroplasmyscatteroutput2)."</td>".
                    "<td> Data file (scatter plots) for circos plot of heteroplasmy result (Normal)</td>".
               "</tr>";
        }else{
            $html.="<tr>".
                    "<td rowspan='2'>Heteroplasmy</td>\n".
                    "<td>"._html_link($mitoheteroplasmy1)."</td>".
                    "<td> Heteroplasmy result (Tumor)</td>".
               "</tr>".
               "<tr>".
                    "<td>"._html_link($mitoheteroplasmy2)."</td>".
                    "<td> Heteroplasmy result (Normal)</td>".
               "</tr>";
          }
          
        #Structural variants
        #Structural variants
         $html.="<tr>".
                    "<td rowspan='4'>Structural Variants</td>\n".
                    "<td>"._html_link($mitostructure1)."</td>".
                    "<td> Structural changes of those discordantly mapped mate reads (Tumor)</td>".
               "</tr>".
               "<tr>".
                    "<td>"._html_link($mitostructure2)."</td>".
                    "<td> Structural changes of those discordantly mapped mate reads (Normal)</td>".
               "</tr>".
               
               "<tr>".
                    "<td>"._html_link($mitostructuredeletion1)."</td>".
                    "<td> Structural changes of those reads mapped with candidate large deletions (Tumor)</td>".
               "</tr>".
                "<tr>".
                    "<td>"._html_link($mitostructuredeletion2)."</td>".
                    "<td> Structural changes of those reads mapped with candidate large deletions (Normal)</td>".
               "</tr>";
        #somatic mutation
        if($cs){
             $html.="<tr>".
                        "<td rowspan='4'>Somatic Mutation</td>".
                        "<td>"._html_link($mitosomatic)."</td>".
                        "<td>Somatic mutation result </td>".
                   "</tr>".
                   "<tr>".
                        "<td>"._html_link("circos/".$mitocircossomaticfigure)."</td>".
                        "<td>Circos plot of somatic mutation result </td>".
                   "</tr>".
                    "<tr>".
                        "<td>"._html_link("circos/".$mitocircossomaticconfig)."</td>".
                        "<td>Configure file for circos plot of somatic mutation result </td>".
                   "</tr>".
                    "<tr>".
                        "<td>"._html_link("circos/".$mitosomatictextoutput)."</td>".
                        "<td>Data file (text label) for circos plot of somatic mutation result </td>".
                   "</tr>";
        }else{
            $html.="<tr>".
                        "<td>Somatic Mutation</td>".
                        "<td>"._html_link($mitosomatic)."</td>".
                        "<td>Somatic mutation result </td>".
                   "</tr>";
        }
        
        #Relative Copy number
        if($cn && $type !=4){
            $html.="<tr>".
                        "<td rowspan='6'> CNV</td>".
                        "<td>"._html_link($mitocnv1)."</td>".
                        "<td> Relative copy number estimation result (Tumor)</td>".
                    "</td>".
                    "<tr>".
                        "<td>"._html_link($mitocnv2)."</td>".
                        "<td> Relative copy number estimation result (Normal)</td>".
                    "</td>".
                    
                    "<tr>".
                        "<td>"._html_link($mitodepth1)."</td>".
                        "<td> Result of depth on mithochondrial genome for CNV estimation (Tumor)</td>".
                    "</tr>".
                    "<tr>".
                        "<td>"._html_link($mitodepth2)."</td>".
                        "<td> Result of depth on mithochondrial genome for CNV estimation (Normal)</td>".
                    "</tr>".
                    
                    "<tr>".
                        "<td>"._html_link($sampledepthi)."</td>".
                        "<td> Result of depth on whole genome for CNV estimation (-i) (Tumor)</td>".
                    "</tr>".
                    "<tr>".
                        "<td>"._html_link($sampledephtj)."</td>".
                        "<td> Result of depth on whole genome for CNV estimation (-j) (Normal)</td>".
                    "</tr>";
        }
        
           #Other
        $html.="<tr>".
                 "<td rowspan='8'> Others </td>".
                 "<td>"._html_link($regionbed)."</td>".
                 "<td> User provided region (bed format) or mitoSeek detected mitochondrial genome region</td>".
               "</tr>".
               "<tr>".
                 "<td>"._html_link($reference)."</td>".
                 "<td> Mitochondrial geneome reference file in fasta format </td>".
               "</tr>".
               
               "<tr>".
                 "<td>"._html_link($mitobam1)."</td>".
                 "<td> Mapped reads in mitochondrial genome (Tumor)</td>".
               "</tr>".
                "<tr>".
                 "<td>"._html_link($mitobam2)."</td>".
                 "<td> Mapped reads in mitochondrial genome (Normal) </td>".
               "</tr>".
               
               "<tr>".
                 "<td>"._html_link($mitopileup1)."</td>".
                 "<td> Pileup file on mithochondrial genome (Tumor)</td>".
               "</tr>".
               "<tr>".
                 "<td>"._html_link($mitopileup2)."</td>".
                 "<td> Pileup file on mithochondrial genome (Normal)</td>".
               "</tr>".
               
               "<tr>".
                 "<td>"._html_link($mitobasecall2)."</td>".
                 "<td> Parsed result of pileup file on mithochondrial genome (Tumor)</td>".
               "</tr>".
               "<tr>".
                 "<td>"._html_link($mitobasecall2)."</td>".
                 "<td> Parsed result of pileup file on mithochondrial genome (Normal)</td>".
               "</tr>";
        
        
    }else{
        #QC related table
        if($qc){
            $html.="<tr>".
                    "<td rowspan='9'>QC</td>\n".
                    "<td>"._html_link($percentofbasepairscovered_table1)."</td>\n".
                    "<td>Percent of base pairs covered</td>\n".
                "</tr>\n".
                "<tr>".
                    "<td>"._html_link($perbasequality_figure1)."</td>\n".
                    "<td> Per base quality plot</td>\n".
               "</tr>\n".
               "<tr>".
                    "<td>"._html_link($perbasequality_table1)."</td>\n".
                    "<td> Data for per base quality plot</td>\n".
                "</tr>".
                "<tr>".
                    "<td>"._html_link($mappingquality_figure1)."</td>\n".
                    "<td> Mapping quality plot</td>\n".
                "</tr>".
                "<tr>".
                    "<td>"._html_link($mappingquality_table1)."</td>\n".
                    "<td> Data for mapping quality plot</td>\n".
                "</tr>".
                "<tr>".
                    "<td>"._html_link($depthdistribution_figure1)."</td>\n".
                    "<td> Depth distribution plot</td>\n".
                "</tr>".
                "<tr>".
                    "<td>"._html_link($depthdistribution_table1)."</td>\n".
                    "<td> Data for depth distribution plot</td>\n".
                "</tr>".
                "<tr>".
                    "<td>"._html_link($templatelengthdistribution_figure1)."</td>\n".
                    "<td> Template length distribution plot</td>\n".
                "</tr>".
                "<tr>".
                    "<td>"._html_link($templatelengthdistribution_table1)."</td>\n".
                    "<td> Data for template length distribution plot</td>\n".
                "</tr>";
          }
        #heteroplasmy
        if($producecircosplot){
            $html.="<tr>".
                    "<td rowspan='5'>Heteroplasmy</td>\n".
                    "<td>"._html_link($mitoheteroplasmy1)."</td>".
                    "<td> Heteroplasmy result </td>".
               "</tr>".
               "<tr>".
                    "<td>"._html_link("circos/".$mitocircosheteroplasmyfigure1)."</td>".
                    "<td> Circos plot of heteroplasmy result </td>".
               "</tr>".
               "<tr>".
                    "<td>"._html_link("circos/".$mitocircosheteroplasmyconfig1)."</td>".
                    "<td> Configure file for circos plot of heteroplasmy result </td>".
               "</tr>".
               "<tr>".
                    "<td>"._html_link("circos/".$mitoheteroplasmytextoutput1)."</td>".
                    "<td> Data file (text labels) for circos plot of heteroplasmy result </td>".
               "</tr>".
               "<tr>".
                    "<td>"._html_link("circos/".$mitoheteroplasmyscatteroutput1)."</td>".
                    "<td> Data file (scatter plots) for circos plot of heteroplasmy result </td>".
               "</tr>";
               
        }else{
            $html.="<tr>".
                    "<td>Heteroplasmy</td>\n".
                    "<td>"._html_link($mitoheteroplasmy1)."</td>".
                    "<td> Heteroplasmy result </td>".
               "</tr>";
        }
        
        #Structural variants
         $html.="<tr>".
                    "<td rowspan='2'>Structural Variants</td>\n".
                    "<td>"._html_link($mitostructure1)."</td>".
                    "<td> Structural changes of those discordantly mapped mate reads </td>".
               "</tr>".
               "<tr>".
                    "<td>"._html_link($mitostructuredeletion1)."</td>".
                    "<td> Structural changes of those reads mapped with candidate large deletions</td>".
               "</tr>";
               
        
        #somatic mutation (No somatic mutation here)
        
        #Relative Copy number
        if($cn && $type !=4){
            $html.="<tr>".
                        "<td rowspan='3'> CNV</td>".
                        "<td>"._html_link($mitocnv1)."</td>".
                        "<td> Relative copy number estimation result</td>".
                    "</td>".
                    "<tr>".
                        "<td>"._html_link($mitodepth1)."</td>".
                        "<td> Result of depth on mithochondrial genome for CNV estimation</td>".
                    "</tr>".
                    "<tr>".
                        "<td>"._html_link($sampledepthi)."</td>".
                        "<td> Result of depth on whole genome for CNV estimation</td>".
                    "</tr>";
        }
        
        #Other
        $html.="<tr>".
                 "<td rowspan='5'> Others </td>".
                 "<td>"._html_link($regionbed)."</td>".
                 "<td> User provided region (bed format) or mitoSeek detected mitochondrial genome region</td>".
               "</tr>".
               "<tr>".
                 "<td>"._html_link($reference)."</td>".
                 "<td> Mitochondrial geneome reference file in fasta format </td>".
               "</tr>".
               "<tr>".
                 "<td>"._html_link($mitobam1)."</td>".
                 "<td> Mapped reads in mitochondrial genome </td>".
               "</tr>".
               "<tr>".
                 "<td>"._html_link($mitopileup1)."</td>".
                 "<td> Pileup file on mithochondrial genome </td>".
               "</tr>".
               "<tr>".
                 "<td>"._html_link($mitobasecall1)."</td>".
                 "<td> Parsed result of pileup file on mithochondrial genome</td>".
               "</tr>";
    }
             
    $html.="</table>\n";
    return $html;
}

sub _html_link{
    my ($in) = @_;
    return "<a href='$in'>$in</a>";
}




#Parse an input file into html table
#
sub _generate_html_table {
    my ($infile, $seperator ) = @_;
    $seperator = "\t" unless ( defined($seperator) );
    my $html ="<div align='right'><a href='$infile'>Download</a>&nbsp;&nbsp;<br/></div>\n";
    
    $html.= "<table border='1'>\n<tr>";
    open( IN, $infile ) or die $!;
    my $line = <IN>;
    $line =~ s/\r|\n//g;
    my @array = split /$seperator/, $line;
    foreach (@array) {
        $html .= "<th>" . $_ . "</th>";
    }
    $html .= "</tr>\n";

    while (<IN>) {
       # s/\r|\n//g;
        @array = split /$seperator/;
        $html .= "<tr>";
        foreach (@array) {
            $html .= "<td>" . $_ . "</td>";
        }
        $html .= "</tr>\n";
    }
    $html .= "</table>";
    return $html;
}


