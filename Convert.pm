#!/usr/bin/perl
#############################################
##Author: Jiang Li
##email: riverlee2008@gmail.com
##Creat Time: Mon Dec  3 10:51:37 CST 2012
##Vanderbilt Center for Quantitative Sciences
##############################################
package Convert;
use strict;
use warnings;
use FindBin;

#This class is used to mitochondrial genome position covertion between hg19 and rCRS
#Information comes from the table ./Resources/hg19_to_rCRS_table.txt
#See the script: hg19_rCRS_convert.pl in the ./Resources see how this table is generated

#Constructor
sub new{
    my ($class) = @_;
    my $self={
        _convertfile=>$FindBin::Bin."/Resources/hg19_to_rCRS_table.txt",
        _hg19TorCRS=>{},
        _rCRSTohg19=>{},
        _hg19ref=>{},
        _rCRSref=>{}
    };
    bless $self,$class;
    return $self;
}

#accessor method for Convert convertfile
sub convertfile{
    my ($self,$var) = @_;
    $self->{_convertfile} = $var if defined($var);
    return $self->{_convertfile};
}

#Initital _hg19TorCRS and _rCRSTohg19 
sub _initial{
    my ($self) = @_;
    open(IN,$self->{_convertfile}) or die "The convert file is not exists\n";
    <IN>;
    while(<IN>){
        s/\r|\n//g;
        my($hg19pos,$hg19ref,$rCRSpos,$rCRSref,$isEqual) = split "\t";
       
        $self->{_hg19TorCRS}->{$hg19pos}=$rCRSpos;
        $self->{_rCRSref}->{$rCRSpos}=$rCRSref;
        $self->{_rCRSTohg19}->{$rCRSpos}=$hg19pos;
        $self->{_hg19ref}->{$hg19pos}=$hg19ref;

    }
    close IN;
}

#Provide a hg19 position, convert it into rCRS position
sub hg19TorCRS{
    my ($self,$var) = @_;
    $self->_initial unless(keys %{$self->{_hg19TorCRS}});
    my $r=undef;
    return $r if (!defined($var) || $var!~/^\d+$/);
    if(exists($self->{_hg19TorCRS}->{$var})){
        $r=$self->{_hg19TorCRS}->{$var};
        $r=undef if ($r eq 'N');
    }
    return $r;
}


sub rCRSTohg19{
    my ($self,$var) = @_;
    $self->_initial unless(keys %{$self->{_hg19TorCRS}});
    my $r=undef;
    return $r if (!defined($var) || $var!~/^\d+$/);
    if(exists($self->{_rCRSTohg19}->{$var})){
        $r=$self->{_rCRSTohg19}->{$var};
        $r=undef if ($r eq 'N');
    }
    return $r;
}

#get the reference allele at a give hg19 position
sub hg19ref{
    my ($self,$var) = @_;
    $self->_initial unless(keys %{$self->{_hg19TorCRS}});
    my $r=undef;
    return $r if (!defined($var) || $var!~/^\d+$/);
    if(exists($self->{_hg19ref}->{$var})){
        $r=$self->{_hg19ref}->{$var};
    }
    return $r;
}

#get the reference allele at a give rCRS9 position
sub rCRSref{
    my ($self,$var) = @_;
    $self->_initial unless(keys %{$self->{_hg19TorCRS}});
    my $r=undef;
    return $r if (!defined($var) || $var!~/^\d+$/);
    if(exists($self->{_rCRSref}->{$var})){
        $r=$self->{_rCRSref}->{$var};
    }
    return $r;
}

1;
=test
#test
package main;
use Convert;
my $c=Convert->new();

my $pos=752;
my $dpos=$c->hg19TorCRS($pos);
print "hg19:pos=$pos\t","rCRS:pos=",$dpos,"\n";
print "hg19:ref=",$c->hg19ref($pos),"\trCRS:ref=",$c->rCRSref($dpos),"\n";

#$pos=752;
#$dpos=$c->rCRSTohg19($pos);
#print "hg19:pos=$dpos\t","rCRS:pos=",$pos,"\n";
#print "hg19:ref=",$c->hg19ref($dpos),"\trCRS:ref=",$c->rCRSref($pos),"\n";
=cut
