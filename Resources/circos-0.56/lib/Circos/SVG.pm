package Circos::SVG;

=pod

=head1 NAME

Circos::Geometry - utility routines for SVG in Circos

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
our @EXPORT = qw(
style_string
);

use Carp qw( carp confess croak );
use FindBin;
use GD::Image;
use Params::Validate qw(:all);

use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";

use Circos::Configuration; 
use Circos::Colors;
use Circos::Constants;
use Circos::Debug;
use Circos::Error;
use Circos::Utils;
use Circos::Image qw(!draw_line);

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

    my %style = ( "stroke-width"   => $params{thickness},
		  "stroke-linecap" => $params{"stroke-linecap"},
		  "stroke"         => sprintf("rgb(%d,%d,%d)",rgb_color($params{color})) );
    
    my $svg = sprintf(qq{<line x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f' style='%s'/>},
		      @{$params{points}},
		      style_string(%style));

    printdebug_group("svg","line",@{$params{points}},$params{color},$params{thickness},style_string(%style));

    Circos::printsvg($svg);

}

################################################################
# Given a hash, generate a style string
#
# key1=value1; key2=value2; key3=value3; ...
sub style_string {
    my %hash = @_;
    my $style = join(";",map { sprintf("%s:%s", $_, $hash{$_}) } keys %hash);
    return $style;
}

1;
