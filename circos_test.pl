#!/usr/bin/perl
#############################################
#Author: Jiang Li
#email: riverlee2008@gmail.com
#Creat Time: 
#Vanderbilt Center for Quantitative Sciences
#############################################
use strict;
use warnings;
use Cwd;
my $currentdir=getcwd;
use Circoswrap;

my $c=Circoswrap->new();

$c->cwd("$currentdir/circosplot");
$c->datafile("$currentdir/mt_1/mito1_heteroplasmy.txt");

$c->prepare("heteroplasmy");
$c->plot();


