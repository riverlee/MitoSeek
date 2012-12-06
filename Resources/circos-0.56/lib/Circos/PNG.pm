package Circos::PNG;

=pod

=head1 NAME

Circos::Geometry - utility routines for PNG in Circos

=head1 SYNOPSIS

This module is not meant to be used directly.

=head1 DESCRIPTION

Circos is an application for the generation of publication-quality,
circularly composited renditions of genomic data and related
annotations.

Circos is particularly suited for visualizing alignments, conservation
and intra and inter-chromosomal relationships. However, Circos can be
used to plot any kind of 2D data in a circular layout - its use is not
limited to genomics. Circos' use of lines to relate position pairs
(ribbons add a thickness parameter to each end) is effective to
display relationships between objects or positions on one or more
scales.

All documentation is in the form of tutorials at L<http://www.circos.ca>.

=cut

# -------------------------------------------------------------------

use strict;
use warnings;

use base 'Exporter';
our @EXPORT = qw();

use Carp qw( carp confess croak );
use FindBin;
use GD;
use Params::Validate qw(:all);

use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";

use Circos::Configuration;
use Circos::Colors;
use Circos::Constants;
use Circos::Debug;
use Circos::Error;
use Circos::Image qw(!draw_line);
use Circos::Utils;

use Memoize;

our $default_color = "black";

for my $f ( qw ( ) ) {
  memoize($f);
}

################################################################
# Draw a line

sub draw_line {
  my %params;
  if( fetch_conf("debug_validate") ) {
    %params = validate(@_,{
			   points           => { type    => ARRAYREF },
			   color            => { default => fetch_conf("default_color") || $default_color  },
			   thickness        => { default => 1 },
			   "stroke-linecap" => { default => "round" },
			  });
  } else {
    %params = @_;
    $params{color}            ||= fetch_conf("default_color") || $default_color;
    $params{"stroke-linecap"} ||= "round";
    $params{thickness}        ||= 1;
  }
  
  if(@{$params{points}} != 4) {
    fatal_error("argument","list_size",current_function(),current_package(),4,int(@{$params{points}}));
  }
  
  printdebug_group("png","line",@{$params{points}},$params{color},$params{thickness});
  
  my $line_color_obj;
  # In the current implementation of gd (2.0.35) antialiasing is
  # incompatible with thick lines and transparency. Thus, antialiased lines
  # are available only when thickness=1 and the color has no alpha channel.
  if($params{thickness} == 1 && rgb_color_opacity($params{color}) == 1) {
    # this is a 1-px thick line and the color has no transparency - 
    # go ahead and antialias this line
    $IM->setAntiAliased(fetch_color($params{color}));
    $line_color_obj = gdAntiAliased;
  } else {
    $IM->setThickness($params{thickness}) if $params{thickness} > 1;
    $line_color_obj = fetch_color($params{color});
  }
  $IM->line( @{$params{points}}, $line_color_obj );
  $IM->setThickness(1) if $params{thickness} > 1;
}

1;
