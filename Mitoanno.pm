#!/usr/bin/perl
#############################################
##Author: Jiang Li
##email: riverlee2008@gmail.com
##Creat Time: Mon Dec  3 10:51:37 CST 2012
##Vanderbilt Center for Quantitative Sciences
##############################################
package Mitoanno;
use strict;
use warnings;
use Convert;
use FindBin;
#This class handles the annotation for a given position in mitochondria genome
#The mitochnodrial annotation comes from NCBI in genbank format, stored at ./Resoures
#For hg19, http://www.ncbi.nlm.nih.gov/nuccore/NC_001807.4?report=genbank
#for rCRS, http://www.ncbi.nlm.nih.gov/nuccore/251831106

#Constructor
sub new{
    my ($class,$var) = @_;
    my $self={
        _build=>"hg19",  #only hg19 and rCRS are accpeted
        _pos2gene=>{},   #position refer to its located gene name
        _gene2pos=>{},   #gene information, includeing start, end,strand, gene name
        _hg19TorCRS=>undef #an Convert object, which is used to access reference allele for a given position
    };
    
    if(defined($var) ){
        if($var eq 'rCRS' || $var eq 'hg19'){
            $self->{_build}=$var;
        }else{
            warn "Only 'hg19' and 'rCRS' are accepted for mitochondria genome, will use 'hg19' as default\n";
        }
    }
    $self->{_gbfile}=$FindBin::RealBin."/Resources/".$self->{_build}."_annotation.gb";

    bless $self,$class;
    return $self;
}

sub genes{
    my ($self) = @_;
    $self->_initial unless (keys %{$self->{_gene2pos}});
    return $self->{_gene2pos};
}

sub build{
    my ($self,$var) = @_;
    if(defined($var)){ 
        if($var eq 'rCRS' || $var eq 'hg19'){
            $self->{_build}=$var;
            $self->{_gbfile}=$FindBin::RealBin."/Resources/".$self->{_build}."_annotation.gb";
            $self->_initial;
        }else{
            warn "Only 'hg19' and 'rCRS' are accepted for mitochondria genome, will use 'hg19' as default\n";
        }
    }
    return $self->{_build};
}

sub gbfile{
    my ($self,$var) = @_;
    if (defined($var)){
        $self->{_gbfile}=$var;
        $self->_initial;
    }
    return $self->{_gbfile};
}


#Given a position, annotation this position
sub annotation{
    my($self,$var,$ref,$alt) = @_;
    
    $self->_initial unless (keys %{$self->{_gene2pos}});

    if(!defined($var) || $var!~/^\d+$/ || $var<1 || $var>16571){
        warn "Please provide the right genome position of mitochondria\n";
    }else{
        #return value
        my ($gene,$genedetail,$aminochange,$aminochangetype);
        if(exists($self->{_pos2gene}->{$var})){
            my @g;
            my @gd;
            my @aminochange;
            my @aminochangetype;
            foreach my $g (@{$self->{_pos2gene}->{$var}}){
                push @g,$g;
                my $tmp=$self->{_gene2pos}->{$g}->{_gene_synonym}.":".
                        $self->{_gene2pos}->{$g}->{_type}.":".
                        $self->{_gene2pos}->{$g}->{_start}."-".
                        $self->{_gene2pos}->{$g}->{_end}.":(".
                        $self->{_gene2pos}->{$g}->{_strand}.")";
                push @gd,$tmp;
                if(defined($ref) && defined($alt) && $self->{_gene2pos}->{$g}->{_type} eq "mRNA"){
                    #Get the codon
                    my $refcodon="";
                    my $altcodon="";
                    my $pos1;
                    my $pos2;
                    if($self->{_gene2pos}->{$g}->{_strand} eq "+"){
                        my $mode = ($var-$self->{_gene2pos}->{$g}->{_start})%3;
                        if($mode==0){
                            $pos1=$var+1;
                            $pos2=$var+2;
                            if($self->{_build} eq "hg19"){
                                $refcodon=$ref.$self->{_hg19TorCRS}->hg19ref($pos1).$self->{_hg19TorCRS}->hg19ref($pos2);
                                $altcodon=$alt.$self->{_hg19TorCRS}->hg19ref($pos1).$self->{_hg19TorCRS}->hg19ref($pos2);
                            }else{
                                $refcodon=$ref.$self->{_hg19TorCRS}->rCRSref($pos1).$self->{_hg19TorCRS}->rCRSref($pos2);
                                $altcodon=$alt.$self->{_hg19TorCRS}->rCRSref($pos1).$self->{_hg19TorCRS}->rCRSref($pos2);
                            }
                        }elsif($mode ==1){
                            $pos1=$var-1;
                            $pos2=$var+1;
                            if($self->{_build} eq "hg19"){
                                $refcodon=$self->{_hg19TorCRS}->hg19ref($pos1).$ref.$self->{_hg19TorCRS}->hg19ref($pos2);
                                $altcodon=$self->{_hg19TorCRS}->hg19ref($pos1).$alt.$self->{_hg19TorCRS}->hg19ref($pos2);
                            }else{
                                $refcodon=$self->{_hg19TorCRS}->rCRSref($pos1).$ref.$self->{_hg19TorCRS}->rCRSref($pos2);
                                $altcodon=$self->{_hg19TorCRS}->rCRSref($pos1).$alt.$self->{_hg19TorCRS}->rCRSref($pos2);
                            }
                        }else{
                            $pos1=$var-2;
                            $pos2=$var-1;
                            if($self->{_build} eq "hg19"){
                                $refcodon=$self->{_hg19TorCRS}->hg19ref($pos1).$self->{_hg19TorCRS}->hg19ref($pos2).$ref;
                                $altcodon=$self->{_hg19TorCRS}->hg19ref($pos1).$self->{_hg19TorCRS}->hg19ref($pos2).$alt;
                            }else{
                                $refcodon=$self->{_hg19TorCRS}->rCRSref($pos1).$self->{_hg19TorCRS}->rCRSref($pos2).$ref;
                                $altcodon=$self->{_hg19TorCRS}->rCRSref($pos1).$self->{_hg19TorCRS}->rCRSref($pos2).$alt;
                            }
                        }
                    }else{
                        my $mode = ($self->{_gene2pos}->{$g}->{_end}-$var)%3;
                       
                        if($mode==0){
                            $pos1=$var-1;
                            $pos2=$var-2;
                            if($self->{_build} eq "hg19"){
                                $refcodon=$ref.$self->{_hg19TorCRS}->hg19ref($pos1).$self->{_hg19TorCRS}->hg19ref($pos2);
                                $altcodon=$alt.$self->{_hg19TorCRS}->hg19ref($pos1).$self->{_hg19TorCRS}->hg19ref($pos2);
                            }else{
                                $refcodon=$ref.$self->{_hg19TorCRS}->rCRSref($pos1).$self->{_hg19TorCRS}->rCRSref($pos2);
                                $altcodon=$alt.$self->{_hg19TorCRS}->rCRSref($pos1).$self->{_hg19TorCRS}->rCRSref($pos2);
                            }
                        }elsif($mode ==1){
                            $pos1=$var+1;
                            $pos2=$var-1;
                            if($self->{_build} eq "hg19"){
                                $refcodon=$self->{_hg19TorCRS}->hg19ref($pos1).$ref.$self->{_hg19TorCRS}->hg19ref($pos2);
                                $altcodon=$self->{_hg19TorCRS}->hg19ref($pos1).$alt.$self->{_hg19TorCRS}->hg19ref($pos2);
                            }else{
                                $refcodon=$self->{_hg19TorCRS}->rCRSref($pos1).$ref.$self->{_hg19TorCRS}->rCRSref($pos2);
                                $altcodon=$self->{_hg19TorCRS}->rCRSref($pos1).$alt.$self->{_hg19TorCRS}->rCRSref($pos2);
                            }
                        }else{
                            $pos1=$var+2;
                            $pos2=$var+1;
                            if($self->{_build} eq "hg19"){
                                $refcodon=$self->{_hg19TorCRS}->hg19ref($pos1).$self->{_hg19TorCRS}->hg19ref($pos2).$ref;
                                $altcodon=$self->{_hg19TorCRS}->hg19ref($pos1).$self->{_hg19TorCRS}->hg19ref($pos2).$alt;
                            }else{
                                $refcodon=$self->{_hg19TorCRS}->rCRSref($pos1).$self->{_hg19TorCRS}->rCRSref($pos2).$ref;
                                $altcodon=$self->{_hg19TorCRS}->rCRSref($pos1).$self->{_hg19TorCRS}->rCRSref($pos2).$alt;
                            }
                        }
                        $refcodon=~tr/ATCG/TAGC/;
                        $altcodon=~tr/ATCG/TAGC/;
                    }
                    
                    my $refamino=codon2aa($refcodon);
                    my $altamino=codon2aa($altcodon);
                    
                    my ($min,undef,$max) = sort {$a<=>$b} ($var,$pos1,$pos2);
                    
                    push @aminochange,$min."-".$max.":".$refcodon."->".$altcodon.":".$refamino."->".$altamino;
                    my $changetype="unknown";
                    if($refamino ne "unknown" && $altamino ne "unknown"){
                        if($refamino eq $altamino){
                            $changetype = "synonymous";
                        }elsif($refamino eq "STOP"){
                            $changetype = "stoploss";
                        }elsif($altamino eq "STOP"){
                            $changetype ="stopgain";
                        }else{
                            $changetype="non-synonymous";
                        }
                    }
                    push @aminochangetype,$changetype;
                }else{
                    push @aminochange,"";
                    push @aminochangetype,"";
                }
            }
            return (join("|",@g),join("|",@gd),join("|",@aminochangetype),join("|",@aminochange));
        }else{
            return("","","","");
        }
    }
}

#if _pos2gene has not been inititalized, load the gbfile 
#For genbank file, we only fetch features of tRNA, rRNA and CDS(mRNA)
sub _initial{
    my ($self) = @_;
    
    #remove the previous annotation
    $self->{_gene2pos}={};
    $self->{_pos2gene}={};
    open(IN,$self->{_gbfile}) or die "The genbank annotation is not exists \n";
    my $feature="";
    my $flag=0;
    while(<IN>){
        if(/^FEATURES/){
            $flag=1;
            next;
        }
        
        last  if(/^ORIGIN/); #The following lines display the sequences
        $feature.=$_ if ($flag);
    }
    close IN;
    
    #Begin parse the $feature;
    # print $feature;
    my @features;
    while($feature=~/^ {5}\S.*\n(^ {21}\S.*\n)*/gm){
        push @features,$&;
    }
    foreach my $f (@features){
        if($f=~/^ {5}(\S+).*\n(^ {21}\S.*\n)*/){
            if($1 eq "tRNA" || $1 eq "rRNA" || $1 eq "CDS"){
                # print $f,"=========\n";
                my ($start,$end,$strand,$type,$geneid,$gene,$gene_synonym)=("","","+",$1,"","","");
                $type="mRNA" if ($type eq "CDS");
                $strand="-" if($f=~/complement\(/);
                if($f=~/^ {5}\S.*?(\d+)\.\.(\d+)/){
                    ($start,$end)=($1,$2);
                }
                if($f=~/\/db_xref="GeneID:(\d+)"\n/){
                    $geneid=$1;
                }
                if($f=~/\/gene="(.*?)"\n/){
                    $gene=$1;
                }
                if($f=~/\/gene_synonym="(.*?)"\n/){
                    $gene_synonym=$1;
                }
                #print join "\t", ($start,$end,$strand,$type,$geneid,$gene,$gene_synonym,"\n");
                $self->{_gene2pos}->{$gene}={_start=>$start,
                                             _end=>$end,
                                             _strand=>$strand,
                                             _type=>$type,
                                             _geneid=>$geneid,
                                             _gene_synonym=>$gene_synonym
                                         };
                #Some of the gene has overlap,
                for(my $i=$start;$i<=$end;$i++){
                    push @{$self->{_pos2gene}->{$i}},$gene;
                }
            }
        }
    }

    $self->{_hg19TorCRS}=Convert->new();
}

sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => 'STOP',    # Stop
    'TAG' => 'STOP',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => 'STOP',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }else{
       print STDERR "Bad codon \"$codon\"!!\n";
            #exit;
      return "unknown";
    }
}
1;

=test
#Test Code
package test;
use Mitoanno;

my $m=Mitoanno->new();
my $ref="T";
my $alt="C";
foreach my $pos (14670..14700){
    print join "\t",($pos,$m->annotation($pos,$ref,$alt));
    print "\n";
}
1;

=cut
