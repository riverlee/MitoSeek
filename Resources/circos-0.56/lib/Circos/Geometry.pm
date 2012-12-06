package Circos::Geometry;

=pod

=head1 NAME

Circos::Geometry - utility routines for Geometry in Circos

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
getxypos
angle_quadrant
);

use Carp qw( carp confess croak );
use FindBin;
use GD::Image;
use Params::Validate qw(:all);

use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";

use Circos::Configuration; # qw(%CONF $DIMS fetch_conf);
use Circos::Constants;
use Circos::Debug;
use Circos::Error;
use Circos::Utils;

use Memoize;

for my $f ( qw ( ) ) {
memoize($f);
}

################################################################
# given an angle, get the xy position for a certain radius
#
sub getxypos {
    return (
	$DIMS->{image}{radius} + $_[1] * cos( $_[0] * $DEG2RAD ),
	$DIMS->{image}{radius} + $_[1] * sin( $_[0] * $DEG2RAD )
	);
}

sub angle_quadrant {
    my $angle = shift;
    if($angle < -90 || $angle > 270) {
	fatal_error("geometry","out_of_bounds",$angle);
    } else {
	if($angle <= 0) {
	    return 0;
	} elsif ($angle <= 90) {
	    return 1;
	} elsif ($angle <= 180) {
	    return 2;
	} else {
	    return 3;
	}
    }
}

1;
