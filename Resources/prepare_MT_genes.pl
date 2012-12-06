#!/usr/bin/perl
use strict;
use warnings;
use lib "../";
use Mitoanno;

my $m=Mitoanno->new();
print $m->gbfile("hg19_annotation.gb"),"\n";
#Output hg19 gene annotation, one for highlights format and the other is in text format for circos plot
my $generef=$m->genes;
open(OUT,">hg19_genes_MT.highlights.txt") or die $!;
open(OUT2,">hg19_genes_MT.text.txt") or die $!;
my $count=0;
foreach my $g (sort keys %{$generef}){
    $count++;
    my $chr=$count % 24;
    $chr++;
    my $color="chr${chr}";
    $color="chr0" if ($color eq 'chr14');
    my $type=$generef->{$g}->{_type};
#    if($type eq 'mRNA'){
    #       $color="chr1";
    #}elsif($type eq 'tRNA'){
    #    $color='chr2';
    # }else{
    #    $color='chr3';
   #}
    print OUT join "\t",("MT",$generef->{$g}->{_start},$generef->{$g}->{_end},"fill_color=$color\n");
    print OUT2 join "\t",("MT",$generef->{$g}->{_start},$generef->{$g}->{_end},$g,"color=black\n");
}


