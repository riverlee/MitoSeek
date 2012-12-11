#!/usr/bin/perl
#############################################
#Author: Jiang Li
#email: riverlee2008@gmail.com
#Creat Time: 
#Vanderbilt Center for Quantitative Sciences
#############################################
use strict;
use warnings;
use FindBin;
use Cwd;
package Circoswrap;

#Circos plot for heterplaysmy (somatic mutations etc, not implemented yet).

#Constructor
sub new{
    my($class) = @_;
    my $currentdir=Cwd::getcwd;
    my $self={
        _build=>"hg19",
        _circosbin=>"$FindBin::Bin/Resources/circos-0.56/bin/circos",     #where is the circos program
        _karyotype=>"$FindBin::Bin/Resources/hg19_karyotype_MT.txt",   #Tell where is the karyotype file of mitochidrial genome
        _genehighlight=>"$FindBin::Bin/Resources/hg19_genes_MT.highlights.txt",  #genes in mitochondrial, by highlight
        _genetext=>"$FindBin::Bin/Resources/hg19_genes_MT.text.txt",             #genes plot in text
        _configtemplate=>"$FindBin::Bin/Resources/circos-heteroplasmy.template.conf",         #default is the template for heteroplasmy plot
        _circosoutput=>"circos-mitoseek.png",           #circos output picture file name
        _configoutput=>"circos-mitoseek.conf",                       #circos conf file output
        _textoutput=>"circos-mitoseek.text.txt",   #circos text data
        _scatteroutput=>"circos-mitoseek.scatter.txt", #circos scatter data
        _datafile=>undef,                               #The file store the heterplaysmy result, will parse it into circos data format
        _cwd=>$currentdir                           #The output dir
    };
    bless $self,$class;
    return $self;
}


sub changeconfig{
    my ($self,$var) = @_;
    if(defined($var)){
        if($var eq "somatic"){
            $self->{_configtemplate} = "$FindBin::Bin/Resources/circos-somatic.template.conf";
        }elsif($var eq "heteroplasmy"){
            $self->{_configtemplate} = "$FindBin::Bin/Resources/circos-heteroplasmy.template.conf";
        }else{
            print ERROR "Undefined template, values should be either 'somatic' or 'heteroplasmy'\n";
            exit(1);
        }
    }
    return $self->{_configtemplate};
}

sub plot{
    my ($self)=@_;
    my $command="cd ".$self->{_cwd}.";".
                $self->{_circosbin}." -conf ".$self->{_configoutput}. ">/dev/null 2>&1";
   # print $command;
   system($command);
}

#Copy the templateconf and change the  into .text and .scatter
#
sub prepare{
    my ($self,$var) = @_;
    
    if(!defined($var)){
        print ERROR "Undefined value, values should be either 'somatic' or 'heteroplasmy'\n";
        exit(1);
    }elsif ( $var ne "somatic" && $var ne "heteroplasmy" ){
        print ERROR "'$var' is not acceptable, values should be either 'somatic' or 'heteroplasmy'\n";
        exit(1);
    }
    
    #copy confi file
    open(IN,$self->{_configtemplate}) or die $!;
    open(OUT,">".$self->{_cwd}."/".$self->{_configoutput}) or die $!;
    my $str="";
    while(<IN>){
        $str.=$_;
    }
    close IN;
    #pattern replace
    $str=~s/replacecircosoutput/$self->{_circosoutput}/g;
    $str=~s/replacegenehighlight/$self->{_genehighlight}/g;
    $str=~s/replacegenetext/$self->{_genetext}/g;
    $str=~s/replacescatter/$self->{_scatteroutput}/g;
    $str=~s/replacetext/$self->{_textoutput}/g;
    print OUT $str;
    close OUT;
    
    #Convert  result into circos result
    if($var eq "heteroplasmy"){
        open(IN,$self->{_datafile}) or die $!; <IN>; #skip header
        open(TEXT,">".$self->{_cwd}."/".$self->{_textoutput}) or die $!;
        open(SCATTER,">".$self->{_cwd}."/".$self->{_scatteroutput}) or die $!;
        while(<IN>){
            my($chr,$pos,$ref,$forward_A,$forward_T,$forward_C,$forward_G,$reverse_A,$reverse_T,$reverse_C,$reverse_G,
               $hetero,$lower,$upper,$major_allele,$minor_allele,$major_count,$minor_count,$gene,$genedetail,$exonic,$aminochange,undef)=split "\t";
            
            my $text="$ref:$major_allele->$minor_allele";
            if(defined($exonic) && ($exonic eq 'non-synonymous' || $exonic eq 'stopgain' || $exonic eq 'stoploss')){
                my (undef,undef,$tmp) = split ":",$aminochange;
                $text.="|$tmp";
            }
            print SCATTER join "\t",("MT",$pos,$pos,$hetero);
            print SCATTER "\n";
            
            print TEXT join "\t",("MT",$pos,$pos,$text);
            print TEXT "\n";
        }
        close IN;
        close TEXT;
        close SCATTER;
    }else{
        #Convert  result into circos result
        open(IN,$self->{_datafile}) or die $!; <IN>; #skip header
        open(TEXT,">".$self->{_cwd}."/".$self->{_textoutput}) or die $!;
        while(<IN>){
            my($chr,$pos,$ref,$tumorGeno,$tumorDP,$normalGeno,$normalDP,$mitoGene,$mitoGeneD)=split "\t";
            
            my $text="$ref:$tumorGeno->$normalGeno";
           
            print TEXT join "\t",("MT",$pos,$pos,$text);
            print TEXT "\n";
        }
        close IN;
        close TEXT;
    }
}

#get or change circos program
sub circosbin{
   my ($self,$var) = @_;
   $self->{_circosbin} = $var if defined($var);
   return $self->{_circosbin};
}

sub cwd{
   my ($self,$var) = @_;
   if(defined($var)){
        mkdir $var or die $! if (! -d  $var);
        $self->{_cwd} = $var;
   }
   return $self->{_cwd};
}

sub circosoutput{
   my ($self,$var) = @_;
   $self->{_circosoutput} = $var if defined($var);
   return $self->{_circosoutput};
}

sub configoutput{
   my ($self,$var) = @_;
   $self->{_configoutput} = $var if defined($var);
   return $self->{_configoutput};
}
sub textoutput{
   my ($self,$var) = @_;
   $self->{_textoutput} = $var if defined($var);
   return $self->{_textoutput};
}
sub scatteroutput{
   my ($self,$var) = @_;
   $self->{_scatteroutput} = $var if defined($var);
   return $self->{_scatteroutput};
}


#set the build or get the build info
#Will auto change the karyotype and gene files
sub build{
    my ($self,$var) = @_;
    if(defined($var)){
        if($var ne "hg19" && $var ne "rCRS"){
            warn "Only 'hg19' and 'rCRS' are accepted for mitochondria genome, will use 'hg19' as default\n";
        }else{
            $self->{_build}=$var;
            $self->{_karyotype}="$FindBin::Bin/Resources/".$var."_karyotype_MT.txt";
            $self->{_genehighlight}="$FindBin::Bin/Resources/".$var."_genes_MT.highlights.txt";
            $self->{_genetext}="$FindBin::Bin/Resources/".$var."_genes_MT.text.txt";
        }
    }
    return $self->{_build};
}

sub karyotype{
    my ($self,$var) = @_;
    $self->{_karyotype}=$var if (defined($var));
    return $self->{_karyotype};
}

sub genehighlight{
    my ($self,$var) = @_;
    $self->{_genehighlight}=$var if (defined($var));
    return $self->{_genehighlight};
}

sub genetext{
     my ($self,$var) = @_;
    $self->{_genetext}=$var if (defined($var));
    return $self->{_genetext};
}

sub configtemplate{
    my ($self,$var) = @_;
    $self->{_configtemplate}=$var if (defined($var));
    return $self->{_configtemplate}; 
}


sub datafile {
    my ($self,$var) = @_;
    $self->{_datafile}=$var if (defined($var));
    return $self->{_datafile}; 
}






 
















1;
