#!/usr/bin/perl
#############################################
#Author: Jiang Li
#email: riverlee2008@gmail.com
#Creat Time: Tue 23 Oct 2012 01:37:54 PM CDT
#Vanderbilt Center for Quantitative Sciences
#############################################
use strict;
use warnings;

#========Begin loading necessary packages===========
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

#=======End loading necessary packages===============

#========Begin defining global variables==============
our $starttime;    #Assign this value in the BEGIN
our %mitogenome = (
	"hg19" => $FindBin::Bin . "/Resources/genome/hg19.fasta",
	"rCRS" => $FindBin::Bin . "/Resources/genome/rCRS.fasta"
);                 #Will be used in the pileup
our %acceptedgenomelength = ( "hg19" => 16571, "rCRS" => 16569 )
  ;                #The total mitochondrial genome length of hg19 and rCRS

our $samtools = "samtools";

#sam/bam flag (will be used in the function _mito_qc_stat)
our $flag_paired                   = 0x0001;
our $flag_properly_paired          = 0x0002;
our $flag_read_unmapped            = 0x0004;
our $flag_next_read_unmapped       = 0x0008;
our $flag_reverse_strand           = 0x0010;
our $flag_next_read_reverse_strand = 0x0020;
our $flag_first_fragment           = 0x0040;
our $flag_second_fragment          = 0x080;

my $debug = 1;

#print $ENV{PATH};
#========End defining global variables=======================

#========Begin defining other variables=============
#Step 1), Assign values to variables
#Variables Assigned by input
my $inbam1 = undef;    #-i
my $inbam2 =
  undef;    #-j ,if this is provided, will conduct somatic mutation mining
my $type              = 1;
my $savebam           = 1;  #-b
my $saveallelecount   = 1;  #-a
my $producecircusplot = 0;  #-ch, produce circus plot for heteroplasmic mutation
my $hp                = 5
  ; #-hp, heteroplasmy threshold using [int] percent alternatie allele observed, default=5;
my $ha = 0; #-ha, heteroplasmy threshold using [int] allele observed, default=0;
my $isall = 0
  ; #If - A is used, the total read count is the total allele count of all allele observed. Otherwise, the total read count is the sum of major and minor allele counts. Default = off
my $mmq = 20;    #minimum map quality, default=20
my $mbq = 20;    #minimum base quality, default=20
my $sb =
  10;   #remove all sites with strand bias score in the top [int] %, default=10;
my $cn = 0
  ; #Estimate relative copy number of input bam(s),does not work with mitochondria targeted sequencing bam files
my $sp = 5
  ; #somatic mutation detection threshold, [int]% of alternative allele observed in tumor, default=5;
my $sa = 3
  ; #somatic mutation detection trheshold, int number of alternative allele observed in tumor.
my $cs = 0
  ; #Produce circus plot input files and circus plot figure for somatic mutations
my $regionbed =
  undef; #A bed file that contains the regions mitoSeek will perform analysis on

#my $inref='hg19';         #The reference used in the input bam files
my $inref  = 'rCRS';
my $outref = 'hg19';    #The output files used in the reference rCRS;
my $qc     = 1;         #Produce QC result
my $str    = 2;         #structure variants cutoff

my $isbam = 1;    #Default is 1, will auto determined from the $inbam1 file

#========End this part=============================

#========Begin assiging values to variables==============
unless (
	GetOptions(
		"i=s"   => \$inbam1,
		"j=s"   => \$inbam2,
		"b!"    => \$savebam,
		"a!"    => \$saveallelecount,
		"ch!"   => \$producecircusplot,
		"hp=i"  => \$hp,
		"ha=i"  => \$ha,
		"A!"    => \$isall,
		"mmq=i" => \$mmq,
		"mbq=i" => \$mbq,
		"sb=i"  => \$sb,
		"cn!"   => \$cn,
		"sp=i"  => \$sp,
		"sa=i"  => \$sa,
		"cs!"   => \$cs,
		"L=s"   => \$regionbed,
		"r=s"   => \$inref,
		"R=s"   => \$outref,
		"QC!"   => \$qc,
		"t=i"   => \$type
	)
  )
{
	print $!, "\n";
	_usage(1);
}

#==========End this part===============================

#=====================Begin Check==================
#1) inbam1 necessary
if ( !defined($inbam1) ) {
	_error("input bam file (-i) has not been defined yet\n");
	_usage();
	exit(1);
}
else {
	if ( !-e $inbam1 ) {
		_error("input bam file (-i) '$inbam1' does not exists\n");
		_usage();
		exit(1);
	}

	#Conver to abs_path
	$inbam1 = abs_path($inbam1) or die $!;
	$isbam = _is_file_binary($inbam1);
}

#2) inbam2 checking
if ( defined($inbam2) ) {
	if ( !-e $inbam2 ) {
		_error("input bam file2 (-j) '$inbam2' does not exists\n");
		_usage();
		exit(1);
	}
	$inbam2 = abs_path($inbam2) or die $!;
}

#3) region file checking
if ( defined($regionbed) ) {
	if ( !-e $regionbed ) {
		_error("bed file (-L) does not exists\n");
		_usage();
		exit(1);
	}
	$regionbed = abs_path($regionbed);
}

#4) inref and outref checking
#print $inref,"\n";
if ( $inref ne 'hg19' && $inref ne 'rCRS' ) {
	_error(
"The reference used in the bam file (-r) should be either hg19 or rCRS, yours is '$inref'\n"
	);
	_usage();
	exit(1);
}
if ( $outref ne 'hg19' && $outref ne 'rCRS' ) {
	_error(
"The reference used in the output file (-R) should be either hg19 or rCRS, yours is '$outref'\n"
	);
	_usage();
	exit(1);
}

#=======================End Check========================

#===================Begin variables will be used==================
my $folder = undef
  ; #determine from the inbam1 file, and will creat this folder, all the report files will be put in this folder
my $mitobam1  = "mito1.bam";    #Reads aligned to mitochondria from inbam1
my $mitobam2  = "mito2.bam";    #Reads aligned to mitochondria from inbam2
my $reference = "mito.fasta"
  ; #$mitogenome{$inref};            #reference sequence in fasta format, used to pileup file

my $mitostart = 1;       #Mitochondrial region start for the analysis
my $mitoend   = undef;

our $perbasequality_figure1             = "per_base_quality1.png";             #
our $perbasequality_table1              = "perl_base_quality_table1.txt";      #
our $mappingquality_figure1             = "mapping_quality1.png";              #
our $mappingquality_table1              = "mapping_quality_table1.txt";        #
our $depthdistribution_figure1          = "depth_distribution1.png";           #
our $depthdistribution_table1           = "depth_distribution_table1.png";     #
our $templatelengthdistribution_figure1 = "template_length_distribution1.png";
our $templatelengthdistribution_table1 =
  "template_length_distribution_table1.txt";
our $perbasequality_figure2             = "per_base_quality2.png";             #
our $perbasequality_table2              = "perl_base_quality_table2.txt";      #
our $mappingquality_figure2             = "mapping_quality2.png";              #
our $mappingquality_table2              = "mapping_quality_table2.txt";        #
our $depthdistribution_figure2          = "depth_distribution2.png";           #
our $depthdistribution_table2           = "depth_distribution_table2.png";     #
our $templatelengthdistribution_figure2 = "template_length_distribution2.png";
our $templatelengthdistribution_table2 =
  "template_length_distribution_table2.txt";

our $mitopileup1   = "mito1.pileup";
our $mitopileup2   = "mito2.pileup";
our $mitobasecall1 = "mito1_basecall.txt";
our $mitobasecall2 = "mito2_basecall.txt";

our $mitoheteroplasmy1 = "mito1_heteroplasmy.txt";
our $mitoheteroplasmy2 = "mito2_heteroplasmy.txt";

our $mitosomatic    = "mito_somatic_mutation.txt";
our $mitostructure1 = "mito1_structure.txt";
our $mitostructure2 = "mito2_structure.txt";

our $mitoreport = "mitoSeek.html";

our $mitooffset1 = 33;
our $mitooffset2 = 33;

( $folder, undef, undef ) = fileparse( $inbam1, qr/\.[^.]*/ );

#print $reference,"\n";
#=================End this part================

#===========Begin Main Part==========================
if ($folder) {
	if ( -d $folder ) {
		_warn("folder $folder already exists,will delete it first");
		rmtree($folder) or die $!;
	}
	mkdir $folder or die $!;
	chdir $folder;
}

my $step = 1;

#No mether its mitochondrial sequences or not, we will to extract reads mapping to mithochondrial first
if ( !$regionbed ) {
	$regionbed = "m.bed";
	_info( $step++
		  . "), Generating mithochondrial bed file from bam file's header (Output: $regionbed)"
	);
	_get_mitochondrial_bed( $inbam1, $regionbed );
}

if ($inbam2) {
	_info( $step
		  . ".1), Extracting reads in mitochondria from $inbam1(tumor) (Output:$mitobam1)"
	);
	_get_mitochondrial_bam( $inbam1, $isbam, $regionbed, $mmq, $mitobam1 );
	_info( $step++
		  . ".2), Extracting reads in mitochondria from $inbam2 (normal) (Output:$mitobam2)"
	);
	_get_mitochondrial_bam( $inbam2, $isbam, $regionbed, $mmq, $mitobam2 );
}
else {
	_info( $step++
		  . "), Extracting reads in mitochondria from $inbam1 (Output: $mitobam1) "
	);
	_get_mitochondrial_bam( $inbam1, $isbam, $regionbed, $mmq, $mitobam1 );
}

# If not reads in the bam, exits
if ($inbam2) {
	_info( $step . ".1), Checking Reads number in $mitobam1(tumor)...", 1 );
	my $num1 = _number_of_reads_in_bam($mitobam1);
	if ( $num1 != 0 ) {
		print $num1, "\n";
	}
	else {
		print "\n";
		_error("No reads in $mitobam1 (tumor)\n");
		exit(1);
	}
	_info( $step++ . ".2), Checking Reads number in $mitobam2(normal)...", 1 );
	my $num2 = _number_of_reads_in_bam($mitobam2);
	if ( $num2 != 0 ) {
		print $num2, "\n";
	}
	else {
		print "\n";
		_error("No reads in $mitobam2 (normal)\n");
		exit(1);
	}
}
else {
	_info( $step++ . "), Checking Reads number in $mitobam1...", 1 );
	my $num1 = _number_of_reads_in_bam($mitobam1);
	if ( $num1 != 0 ) {
		print $num1, "\n";
	}
	else {
		print "\n";
		_error("No reads in $mitobam1 \n");
		exit(1);
	}
}

# Prepare mitochondrial reference, change header in fasta when necessary
_info(  $step++
	  . "), Moving mitochondrial genome '"
	  . $mitogenome{$inref}
	  . "' into $folder (Output:$reference)" );
_move_mitogenome( $mitogenome{$inref}, $regionbed, $reference );

# Pileup on bam files
if ($inbam2) {
	_info( $step . ".1), Pileup on $mitobam1 (tumor)(Output:$mitopileup1)" );
	_pileup( $mitobam1, $isbam, $mbq, $regionbed, $reference, $mitopileup1 );
	_info( $step++ . ".2), Pileup on $mitobam2 (normal)(Output:$mitopileup2)" );
	_pileup( $mitobam2, $isbam, $mbq, $regionbed, $reference, $mitopileup2 );
}
else {
	_info( $step++ . "), Pileup on $mitobam1 (Output:$mitopileup1)" );
	_pileup( $mitobam1, $isbam, $mbq, $regionbed, $reference, $mitopileup1 );
}

# Running QC on mitobam
if ($qc) {
	if ($inbam2) {
		_info( $step . ".1), Getting quality metrics on $mitobam1 (tumor) " );
		$mitooffset1 = _mito_qc_stat(
			$mitobam1,                  $isbam,
			$mitopileup1,               $regionbed,
			$perbasequality_figure1,    $mappingquality_figure1,
			$depthdistribution_figure1, $templatelengthdistribution_figure1
		);
		_info(
			$step++ . ".2), Getting quality metrics on $mitobam2 (normal) " );
		$mitooffset2 = _mito_qc_stat(
			$mitobam2,                  $isbam,
			$mitopileup2,               $regionbed,
			$perbasequality_figure2,    $mappingquality_figure2,
			$depthdistribution_figure2, $templatelengthdistribution_figure2
		);

	}
	else {
		_info( $step++ . "), Getting quality metrics on $mitobam1  " );
		$mitooffset1 = _mito_qc_stat(
			$mitobam1,                  $isbam,
			$mitopileup1,               $regionbed,
			$perbasequality_figure1,    $mappingquality_figure1,
			$depthdistribution_figure1, $templatelengthdistribution_figure1
		);
	}
}

# Parse pileup file
if ($inbam2) {
	_info( $step
		  . ".1), Parsing pileup file of $mitopileup1 (tumor) (Output: $mitobasecall1)"
	);
	_parse_pileup( $mitopileup1, $mbq, $mitooffset1, $mitobasecall1 );
	_info( $step++
		  . ".2), Parsing pileup file of $mitopileup2 (normal) (Output: $mitobasecall2)"
	);
	_parse_pileup( $mitopileup2, $mbq, $mitooffset2, $mitobasecall2 );
}
else {
	_info( $step++
		  . "), Parsing pileup file of $mitopileup1 (Output: $mitobasecall1)" );
	_parse_pileup( $mitopileup1, $mbq, $mitooffset1, $mitobasecall1 );
}

# heteroplasmy detection
if ($inbam2) {
	_info( $step
		  . ".1), Heteroplasmy detection from $inbam1 (tumor) (Output: $mitoheteroplasmy1)"
	);
	_determine_heteroplasmy( $mitobasecall1, $hp, $ha, $isall, $sb,
		$mitoheteroplasmy1 );
	_info( $step++
		  . ".2), Heteroplasmy detection from $inbam2 (normal) (Output: $mitoheteroplasmy2)"
	);
	_determine_heteroplasmy( $mitobasecall2, $hp, $ha, $isall, $sb,
		$mitoheteroplasmy2 );
}
else {
	_info( $step++
		  . "), Heteroplasmy detection from $inbam1 (Output: $mitoheteroplasmy1)"
	);
	my %hash = _determine_heteroplasmy( $mitobasecall1, $hp, $ha, $isall, $sb,
		$mitoheteroplasmy1 );
}

# somatic mutation
if ($inbam2) {
	_info( $step++ . ") Somatic mutation detection (Output: $mitosomatic)" );
	_determine_somatic( $mitobasecall1, $mitobasecall2, $sp, $sa, $isall,
		$mitosomatic );
}

# Structure variation
if ($inbam2) {
	_info( $step
		  . ".1) Structure variation detection from $mitobam1 (Output: $mitostructure1)"
	);
	_structure_variants( $inbam1, $isbam, $mmq, $regionbed, $str,
		$mitostructure1 );
	_info( $step++
		  . ".2) Structure variation detection from $mitobam2 (Output: $mitostructure2)"
	);
	_structure_variants( $inbam2, $isbam, $mmq, $regionbed, $str,
		$mitostructure2 );
}
else {
	_info( $step++
		  . ") Structure variation detection from $mitobam1 (Output: $mitostructure1)"
	);
	_structure_variants( $inbam1, $isbam, $mmq, $regionbed, $str,
		$mitostructure1 );
}

_info( $step++ . ") Generating report (Output: $mitoreport)" );
_make_report();

#============End main part=============================

#===========Start BEGIN and END section===================
BEGIN {
	$starttime = time();
}

END {
	my $interval = time() - $starttime;

	#my $interval=192232;
	my $hours = int( $interval / 3600 );  # calculate precise, and then floor it
	my $minutes = int( ( $interval - $hours * 3600 ) / 60 );
	my $seconds = $interval % 60;
	my $time    = sprintf( "%dh:%02dm:%02ds", $hours, $minutes, $seconds );
	print "Total Running Time: $time\n";
}

#===========End begin and end section======================

#============Begin defining functions=====================

sub _structure_variants {
	my ( $inmitobam, $isbam, $mmq, $regionbed, $str, $outfile ) = @_;
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
	print OUT join "\t",
	  (
		"#MitoChr",       "MitoPos",
		"VarChr",         "VarPos",
		"SupportedReads", "MeanMappingQuality\n"
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
	}
	close IN;
	foreach my $k ( sort keys %temphash ) {
		if ( $temphash{$k}->{count} >= $str ) {
			print OUT join "\t",
			  (
				$k,
				$temphash{$k}->{'count'},
				_mymean( $temphash{$k}->{'mapping'} )
			  );
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
		return $t / scalar( @{$arrayref} );
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

#Get bam file number
#Parameter
#$inbam
sub _number_of_reads_in_bam {
	my ($inbam) = @_;
	my $num = 0;
	$num = `$samtools view $inbam|wc -l`;
	$num =~ s/\r|\n//g;
	return $num;
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
	#	print "Text\n";
	#}
	return $flag;
}

#Print out usage,
#Parameters:
#$flag, if $flag is set and it will exit the program
sub _usage {
	my ($flag) = @_;
	my $usage = <<USAGE;
Usage: perl mitoSeek.pl -i inbam.pl
-i [bam]                Input bam file
-j [bam]                Input bam file2, if this file is provided, it will conduct somatic mutation mining, and it will be 
                        taken as normal tissue.        
-t [input type]         Type of the bam files, the possible choices are 1=exome, 2=whole genome, 3= RNAseq, 4 = mitochondria only
-nob                    Save mitochondria only bam files, default =on, if this is provided,it will turn off
-noa                    Produce allele count file for input bam(s), default =on, if this is provided,it will turn off
-ch                     Produce circus plot input files and circus plot figure for heteroplasmic mutation, default = off
-hp [int]               Heteroplasmy threshold using [int] percent alternative allele observed, default = 5
-ha [int]               Heteroplasmy threshold using [int] allele observed, default = 0
-A                      If - A is used, the total read count is the total allele count of all allele observed. 
                        Otherwise, the total read count is the sum of major and minor allele counts. Default = off
-mmq [int]              Minimum map quality, default =20
-mbq [int]              Minimum base quality, default =20
-sb [int]               Remove all sites with strand bias score in the top [int] %, default = 10 
-cn                     Estimate relative copy number of input bam(s), does not work with mitochondria targeted sequencing bam files.
-s [bam1] [bam2]        Compute somatic mutation between bam1 and bam2
-sp [int1][int2]        Somatic mutation detection threshold, int1 = percent of alternative allele observed in normal, 
                         int2 = percent of alternative allele observed in tumor, default int1=0, int2=5
-sa [int1][int2]        Somatic mutation detection threshold, int1 = number of alternative allele observed in normal,
                         int2 = number of alternative allele observed in tumor, default int1=0, int2=3
-cs                     Produce circus plot input files and circus plot figure for somatic mutation, default = off
-L [bed]                A bed file that contains the regions MitoSeek will perform analysis on
-r [ref]                The reference used in the bam file, the possible choices are HG19 and rCRS, default=HG19
-R [ref]                The reference used in the output files, the possible choices are HG19 and rCRS, default=HG19
-QC                     Produce QC result,default=on
   
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
		$depthdistribution_table,  $templatelengthdistribution_table
	) = @_;
	my $maxreads = 1000000
	  ; #In case some alignment on mitochondrial is extremly large that will take too much memory and then crash your server.
	$maxreads = 1000 if ($debug);

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
		0, undef, 30
	) if ($isparied);

	#2) persequencequality density distribution
	_density(
		\@persequencequality, 800, 600, 'Quality score',
		'Density', "Per sequence quality socre distribution ($encoding)",
		$mappingquality_figure, 0, undef, 30
	);

	#3) perbasequality boxplot
	_boxplot(
		\@perbasequality, 800, 600, 'Bases',
		'Quality score',
		"Per base quality socre distribution ($encoding)",
		$perbasequality_figure, undef, undef, 30
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
		$depthdistribution_figure, 0, undef, 30 );

	return ($offset);
}

#Function to plot boxplot
#Parameters:
sub _boxplot {
	my (
		$data,  $width,  $height, $xlabel, $ylabel,
		$title, $output, $xmin,   $xmax,   $nbin
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
		r_margin          => 20
	) or warn $graph->error;
	my $gd = $graph->plot($data) or die $graph->error;
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
#$bin  -----How many bins at x-axis,Default is 20
sub _density {
	my (
		$data,  $width,  $height, $xlabel, $ylabel,
		$title, $output, $xmin,   $xmax,   $nbin
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
		r_margin          => 20
	) or warn $graph->error;
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
Min.=$minv
25%=$v25
Mean=$mean
Median=$median
75%=$v75
Max.=$maxv
Var=$var
EOSTR
	my $gdat = GD::Text::Wrap->new(
		$gd,
		text       => $text,
		line_space => 4,
		align      => 'left',
		width      => 5
	);

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
	open( OUT, ">$outbed" );
	print OUT join "\t", ( $m, 1, $len );
	close OUT;
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

		#Calculate heteroplasmy ratio
		my $heteroplasmy = 0;
		if ($isall) {
			$heteroplasmy =
			  $atcg{$minor_allele} /
			  ( $atcg{A} + $atcg{C} + $atcg{T} + $atcg{G} );
		}
		else {
			$heteroplasmy =
			  $atcg{$minor_allele} /
			  ( $atcg{$major_allele} + $atcg{$minor_allele} );
		}

		#Calculate 95% confidence interval
		my ( $lower, $upper ) = ( 0, 0 );
		if ( $heteroplasmy > 0 ) {
			my $n = 0;
			if ($isall) {
				$n = $atcg{A} + $atcg{C} + $atcg{T} + $atcg{G};
			}
			else {
				$n = $atcg{$major_allele} + $atcg{$minor_allele};
			}
			( $lower, $upper ) = _ci( $heteroplasmy, $n );
		}

		#Stat values could be
		#0 (not a heteroplasmy mutation)
		#1 (show heteroplasmy, but not pass the cutoff)
		#2 (heteroplasmy mutation pass cutoff and not have strong strand bias)
		#3 (heteroplasmy mutation pass cutoff but show strong strand bias)
		my $stat = 0;
		if ( $heteroplasmy > 0 ) {

#determine 1 or 2 at this stage (for value 3, need to first calculate strand bias across all the site, and then modify those with stat=2, only if when $sb is not equal to 0)
			if ( $heteroplasmy > $hp / 100 && $atcg{$minor_allele} > $ha ) {
				$stat = 2;
			}
			else {
				$stat = 1;
			}
		}

#If provided strand bias cutoff, need to re-loop the @result, change state=2 to 3 when necessary.
#		push @result,
#		  [
#			$chr,          $loc,          $ref,
#			$forward_A,    $forward_T,    $forward_C,
#			$forward_G,    $reverse_A,    $reverse_T,
#			$reverse_C,    $reverse_G,    $stat,
#			$heteroplasmy, $lower,        $upper,
#			$major_allele, $minor_allele, $atcg{$major_allele},
#			$atcg{$minor_allele}
#		  ];
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
		if ( $sb > 0 ) {
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

#Now whether isall=0 or 1, strand bias is calculated based on major allele and minor allele
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
		}
	}

	#Re-loop the @result and write out
	#	if ( $sb > 0 ) {
	#		my @sbvalue = sort { $b <=> $a } values %sb;
	#		my $index = int( scalar(@sbvalue) * $sb / 100 ) - 1;
	#		$index = 0         if ( $index < 0 );
	#		$index = $#sbvalue if ( $index > $#sbvalue );
	#		my $cutoff = $sbvalue[$index];
	#		foreach my $arrayref (@result) {
	#			if (   $sb{ $arrayref->[1] } >= $cutoff
	#				&& $arrayref->[11] == 2 )
	#			{
	#				$arrayref->[11] = 3;
	#			}
	#		}
	#	}
	if ( $sb > 0 ) {
		my @sbvalue = sort { $b <=> $a } values %sb;
		my $index = int( scalar(@sbvalue) * $sb / 100 ) - 1;
		$index = 0         if ( $index < 0 );
		$index = $#sbvalue if ( $index > $#sbvalue );
		my $cutoff = $sbvalue[$index];

		foreach my $loc ( keys %result ) {
			$result{$loc}->{'sb'} = $sb{$loc};
			if ( $sb{$loc} >= $cutoff && $result{$loc}->{'stat'} == 2 ) {
				$result{$loc}->{'stat'} = 3;
			}
		}
	}

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
		"major_allele_count", "minor_allele_count"
	  );
	if ($sb) {
		print OUT "\tstrand_bias";
	}
	print OUT "\n";

	foreach my $loc ( sort { $a <=> $b } keys %result ) {
		if ( $result{$loc}->{'stat'} == 2 ) {
			print OUT join "\t",
			  (
				$result{$loc}->{'chr'},
				$loc,
				$result{$loc}->{'ref'},
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
			if ( $sb > 0 ) {
				print OUT "\t", $result{$loc}->{'sb'};
			}
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
	return ( $lower, $upper );
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
		"tumorDepth", "normalGenotype", "normalDepth\n"
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

		if ( $normal{$loc}->{'minor_allele_count'} == 0 ) {
			if ( $tumor{$loc}->{'minor_allele_count'} == 0 ) {

#If at this site, both normal and tumor are not heteroplasmy, then determine whether their
#major_allele is the same
				if ( $normal{$loc}->{'major_allele'} ne
					$tumor{$loc}->{'major_allele'} )
				{

					#a somatic mutation
					$issomatic = 1;
				}
			}
			else {
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
			print OUT join "\t",
			  (
				"MT", $loc, $normal{$loc}->{'reference_allele'},
				$tumorGeno, $tumorDP, $normalGeno, $normalDP
			  );
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
"$samtools mpileup -l $regionbed -Q $mbq -f $refseq $inbam > $outpileup";
	}
	else {
		$comm =
"$samtools view -Su $inbam | $samtools mpileup -l $regionbed -Q $mbq -f $refseq - > $outpileup";
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
		  abs( $b / ( $a + $b ) - $d / ( $c + $d ) ) /
		  ( ( $b + $d ) / ( $a + $b + $c + $d ) );
	}
}

#Check whether samtools exist in your $PATH, If not, the program will exit
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
	system($comm) == 0 or die $!;
}




##Generate the report
sub _make_report {
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
  font-size: 25%;
  margin-right:2em;
  text-align: right;
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
min-width:800px; /* width of 5 images (625px) plus images' padding and border (60px) plus images' margin (50px) */ 
height:47px; 
padding:10px 5px; 
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
width:400px; 
height:300px; 
padding:5px; 
border:1px solid #d3d2d2; 
margin:auto; 
background-color:#fff; 
} 

#image-container1 { 
min-width:600px; /* width of 5 images (625px) plus images' padding and border (60px) plus images' margin (50px) */ 
height:47px; 
padding:10px 5px; 
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
width:600px; 
height:450px; 
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
</style>
</head>
<body>
<div class="header">
	<div id="header_title">MitoSeek Report</div>
	<!--div id="header_filename"></div-->
</div>

<div class="summary">
	<h2>Table of Contents</h2>
		<ul>
			<li> <a href="#M0">QC</a></li>
				<ol>
					<li> <a href="#M00">Per base quality</a></li>
					<li> <a href="#M01">Mapping quality</a></li>
					<li> <a href="#M02">Depth distribution</a></li>
					<li> <a href="#M03">Template length distribution</a></li>
				</ol>
			<li> <a href="#M1">Heteroplasmy</a></li>
			<li> <a href="#M2">Somatic Mutation</a></li>
			<li> <a href="#M3">Relative Copy Number Estimation</a></li>
			<li> <a href="#M4">Structural Changes</a></li>
			<li> <a href="#M5">File List</a></li>
		</ul>
</div>
<div class="main">
	<div class="module"><h2 id="M0">QC</h2>
		<p>QC description here</p>
		<h3 id="M00">Per base quality</h3>
		<p>Per base quality description here</p>
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
		<p>Mapping quality here</p>
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
		<p>Depth distribution here</p>
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
		<p>Template length distribution here</p>
HTML
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

print OUT <<HTML;
	</div>
	<div class="module"><h2 id="M1">Heteroplasmy</h2>
		<p>Heteroplasmy description</p>
HTML

	if ($inbam2) {
		print OUT "<h3>Tumor</h3>",  _generate_html_table($mitoheteroplasmy1);
		print OUT "<h3>Normal</h3>", _generate_html_table($mitoheteroplasmy2);
	}
	else {
		print OUT "<h3>AAA</h3>", _generate_html_table($mitoheteroplasmy1);
	}

print OUT <<HTML;
	</div>
	<div class="module"><h2 id="M2">Somatic Mutation</h2>
	<p>Somatic mutation description</p>
HTML
	if ($inbam2) {
		print OUT _generate_html_table($mitosomatic);
	}
	else {
		print OUT "NA\n";
	}

	print OUT <<HTML;
	</div>
	<div class="module"><h2 id="M3">Relative Copy Number Estimation</h2>
	<p>Relative copy number estimation description</p>
	<b>Not implemented yet</b>
	</div>
HTML

	print OUT <<HTML;
	<div class="module"><h2 id="M4">Structural Changes</h2>
	<p>Structural changes description</p>
HTML
	if ($inbam2) {
		print OUT "<h3>Tumor</h3>",  _generate_html_table($mitostructure1);
		print OUT "<h3>Normal</h3>", _generate_html_table($mitostructure2);
	}
	else {
		print OUT "<h2>AAA</h2>", _generate_html_table($mitostructure1);
	}


	print OUT <<HTML;
	</div>
	<div class="module"><h2 id="M5">File list</h2>
	<p>File list description</p>
	<b>Not implemented yet</b>
	</div>
HTML

	print OUT <<HTML;
	</div>
	</div>
<div class="footer"  style="background-color: #EEE;font-weight:" align='center'>Produced by <a href="https://github.com/riverlee/MitoSeek">MitoSeek</a> (Jiang Li)</div>
</body>
</html>
HTML


	close OUT;
}


##Generate the report
sub _make_report_old {
	open( OUT, ">$mitoreport" ) or die $!;
	print OUT <<HTML;
<html>
<head><title>MitoSeek Report</title>
<style type="text/css">
table {
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
}
table td {
	border-width: 2px;
	padding: 2px;
	border-style: solid;
	border-color: gray;
	background-color: white;
}
</style>
</head>
<body>
<div  style="background-color:#EEE;font-weight: bold;  bold;font-size:32px" align='left'>MitoSeek Report</div>
<h1>QC</h1>
<h2>Per base quality<h2>
<img src="$perbasequality_figure1" />
<h2>Mapping quality</h2>
<img src="$mappingquality_figure1" />
<h2>Depth distribution</h2>
<img src="$depthdistribution_figure1" />

<h1>Heteroplasmy</h1>
<p>description later</p>
HTML
	if ($inbam2) {
		print OUT "<h2>Tumor</h2>",  _generate_html_table($mitoheteroplasmy1);
		print OUT "<h2>Normal</h2>", _generate_html_table($mitoheteroplasmy2);
	}
	else {
		print OUT "<h2>AAA</h2>", _generate_html_table($mitoheteroplasmy1);
	}

	print OUT <<HTML;
<h1>Somatic Mutation</h1>
HTML

	if ($inbam2) {
		print OUT _generate_html_table($mitosomatic);
	}
	else {
		print OUT "NA\n";
	}

	print OUT <<HTML;
<h1>Relative Copy Number Estimation</h1>
HTML

	print OUT <<HTML;
<h1>Structural Change</h1>
HTML
	if ($inbam2) {
		print OUT "<h2>Tumor</h2>",  _generate_html_table($mitostructure1);
		print OUT "<h2>Normal</h2>", _generate_html_table($mitostructure2);
	}
	else {
		print OUT "<h2>AAA</h2>", _generate_html_table($mitostructure1);
	}

	print OUT <<HTML;
</div>
<div style="background-color: #EEE;font-weight:" align='center'>Produced by <a href="https://github.com/riverlee/MitoSeek">MitoSeek</a> (Jiang Li)</div>
</body>
</html>

HTML

	close OUT;
}

#Parse an input file into html table
#
sub _generate_html_table {
	my ($infile, $seperator ) = @_;
	$seperator = "\t" unless ( defined($seperator) );
	my $html = "<table border='1'>\n<tr>";
	open( IN, $infile ) or die $!;
	my $line = <IN>;
	$line =~ s/\r|\n//g;
	my @array = split /$seperator/, $line;
	foreach (@array) {
		$html .= "<th>" . $_ . "</th>";
	}
	$html .= "</tr>\n";

	while (<IN>) {
		s/\r|\n//g;
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

