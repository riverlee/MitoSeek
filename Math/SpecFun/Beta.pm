# -*- mode: perl; coding: utf-8-unix -*-
#
# Author:      Peter John Acklam
# Time-stamp:  2009-03-30 16:22:42 +02:00
# E-mail:      pjacklam@online.no
# URL:         http://home.online.no/~pjacklam

=pod

=head1 NAME

Math::SpecFun::Beta - beta, log(beta) and incomplete beta functions

=head1 SYNOPSIS

    use Math::SpecFun::Beta qw(beta betaln betainc betaincc);

Imports all the routines explicitly.  Use a subset of the list for the
routines you want.

    use Math::SpecFun::Beta qw(:all);

Imports all the routines, as well.

=head1 DESCRIPTION

This module implements the beta function, C<beta>, the natural
logarithm of the beta function, C<betaln>, the incomplete beta
function C<betainc>, and the complement of the incomplete beta
function C<betaincc>.

For references and details about the algorithms,  see below and also
see the comments inside this module.

=head1 FUNCTIONS

=over 8

=item beta A, B

Returns the beta function beta(a,b).  Both A and B must be positive.
The beta function is defined as

   beta(a,b) = integral from 0 to 1 of t^(a-1) * (1-t)^(b-1) dt

The following identity holds

   beta(a,b) = beta(b,a)

=item betaln A, B

Returns the logarithm of the beta function log(beta(a,b)).  Both A and B
must be positive.  The logarithm of the beta function is defined as

   betaln(a, b) = log(beta(a, b))

=item betainc X, A, B

Returns the incomplete beta function betainc(x,a,b).  Both A and B must be
positive; X must be in the closed interval [0,1].  The incomplete beta
function is defined as

    betainc(x,a,b) = 1 / beta(a,b) *
                       integral from 0 to x of t^(a-1) * (1-t)^(b-1) dt

The following identity holds

    betainc(x,a,b) = 1 - betainc(1-x,b,a)
                   = 1 - betaincc(x,a,b)
                   = betaincc(1-x,b,a)

=item betaincc X, A, B

Returns the complement of the incomplete beta function
betaincc(x,a,b).  Both A and B must be positive; X must be in the
closed interval [0,1].  The complement of the incomplete beta function
is defined as

    betaincc(x,a,b) = 1 / beta(a,b) *
                        integral from x to 1 of t^(a-1) * (1-t)^(b-1) dt

The following identity holds

    betaincc(x,a,b) = 1 - betaincc(1-x,b,a)
                    = 1 - betainc(x,a,b)
                    = betainc(1-x,b,a)

=back

=head1 HISTORY

=over 4

=item Version 0.01

First release.

=back

=head1 AUTHORS

Perl translation by Peter John Acklam E<lt>pjacklam@online.noE<gt>

=head1 COPYRIGHT

Copyright 2000 Peter John Acklam.
This program is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

=cut

package Math::SpecFun::Beta;
require 5.000;
require Exporter;

use strict;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
use Carp;

$VERSION = '0.02';

@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(beta betaln betainc);
%EXPORT_TAGS = (all => [ @EXPORT_OK ]);

# the smallest positive floating-point number such that 1+EPS > 1
# Numerical Recipes uses 3.0e-7
use constant EPS   => 2.220446049250313080847263336181640625e-016;

# the smallest positive normalized floating point number
# Numerical Recipes uses 1.0e-30
use constant FPMIN => 2.2250738585072e-308;

# largest number of iterations we allow
use constant MAXIT => 100;

########################################################################
## Internal functions.
########################################################################

########################################################################
# Evaluates continued fraction for incomplete beta function by
# modified Lentz's method (sec 5.2 in Numerical Recipes).
########################################################################

sub _betacf ($$$) {
    # This is an internal function, so we omit error checking of args.
    my ($x, $a, $b) = @_;

    # These q's will be used in factors that occur
    # in the coefficients (6.4.6).
    my $qab = $a+$b;
    my $qap = $a+1;
    my $qam = $a-1;
    my $c = 1;                  # First step of Lentz's method.
    my $d = 1-$qab*$x/$qap;
    if (abs($d) < FPMIN) { $d = FPMIN; }
    $d = 1/$d;
    my $h = $d;
    my $m;
    for ( $m = 1; $m <= MAXIT; $m++ ) {
        my $m2 = 2*$m;
        my $aa = $m*($b-$m)*$x/(($qam+$m2)*($a+$m2));
        $d = 1+$aa*$d;          # One step (the even one) of the recurrence.
        if (abs($d) < FPMIN) { $d = FPMIN; }
        $c = 1+$aa/$c;
        if (abs($c) < FPMIN) { $c = FPMIN; }
        $d = 1/$d;
        $h *= $d*$c;
        $aa = -($a+$m)*($qab+$m)*$x/(($a+$m2)*($qap+$m2));
        $d = 1+$aa*$d;          # Next step of the recurrence (the odd one).
        if (abs($d) < FPMIN) { $d = FPMIN; }
        $c = 1+$aa/$c;
        if (abs($c) < FPMIN) { $c = FPMIN; }
        $d = 1/$d;
        my $del = $d*$c;
        $h *= $del;
        if (abs($del-1) < EPS) { last; }   # Are we done?
    }
    if ($m > MAXIT) {
        die "a or b too big, or MAXIT too small in _betacf";
    }
    return $h;
}

########################################################################
## User functions.
########################################################################

########################################################################
# Beta function
#
# Reference: Abramowitz & Stegun, Handbook of Mathematical Functions,
# sec. 6.2.
########################################################################

sub beta {
    # Get and check input arguments.
    croak "usage: beta(a,b)\n" unless @_ == 2;
    my ($a, $b) = @_;
    croak "Both arguments to beta must be positive"
      unless $a > 0 && $b > 0;

    # Special cases.
    return 1/$b if $a == 1;
    return 1/$a if $b == 1;

    # The general case.
    require Math::SpecFun::Gamma;
    import Math::SpecFun::Gamma qw(gammaln);
    return exp( gammaln($a) + gammaln($b) - gammaln($a+$b) );
}

########################################################################
# Logarithm of beta function.
#
# Reference: Abramowitz & Stegun, Handbook of Mathematical Functions,
# sec. 6.2.
########################################################################

sub betaln {
    # Get and check input arguments.
    croak "usage: betaln(a,b)\n" unless @_ == 2;
    my ($a, $b) = @_;
    croak "Both arguments to betaln must be positive"
      unless $a > 0 && $b > 0;

    # Special cases.
    return -log($b) if $a == 1;
    return -log($a) if $b == 1;

    # The general case.
    require Math::SpecFun::Gamma;
    import Math::SpecFun::Gamma qw(gammaln);
    return gammaln($a) + gammaln($b) - gammaln($a+$b);
}

########################################################################
# Incomplete beta function.
#
# Reference: Abramowitz & Stegun, Handbook of Mathematical Functions,
# sec. 26.5.
########################################################################

sub betainc {
    # Get and check input arguments.
    croak "usage: betainc(x,a,b)\n" unless @_ == 3;
    my ($x, $a, $b) = @_;
    croak "First argument must be in the interval [0,1]"
      if $x < 0 || $x > 1;
    croak "Last two arguments must be positive" unless $a > 0 && $b > 0;

    # Special cases.
    return 0                if $x == 0;
    return 1                if $x == 1;
    return $x**$a           if $b == 1;
    return 1 - (1 - $x)**$b if $b == 1;
    return 0.5              if $x == 0.5 && $a == $b;

    # The general case.
    require Math::SpecFun::Gamma;
    import Math::SpecFun::Gamma qw(gammaln);
    # Factors in front of the continued fraction.
    my $bt = exp( gammaln($a + $b) - gammaln($a) - gammaln($b)
                  + $a*log($x) + $b*log(1 - $x) );
    if ( $x < ($a+1)/($a+$b+2) ) {
        # Use continued fraction directly.
        return $bt * _betacf($x, $a, $b) / $a;
    } else {
        # Use continued fraction after making the symmetry
        # transformation.
        return 1 - $bt * _betacf(1 - $x, $b, $a) / $b;
    }
}

########################################################################
# The complement of the incomplete beta function.
########################################################################

sub betaincc {
    # Get and check input arguments.
    croak "usage: betaincc(x,a,b)\n" unless @_ == 3;
    my ($x, $a, $b) = @_;
    croak "First argument must be in the interval [0,1]"
      if $x < 0 || $x > 1;
    croak "Last two arguments must be positive" unless $a > 0 && $b > 0;

    # Special cases.
    return 1                if $x == 0;
    return 0                if $x == 1;
    return 1 - $x**$a       if $b == 1;
    return (1 - $x)**$b     if $b == 1;
    return 0.5              if $x == 0.5 && $a == $b;

    # The general case.
    require Math::SpecFun::Gamma;
    import Math::SpecFun::Gamma qw(gammaln);
    # Factors in front of the continued fraction.
    my $bt = exp( gammaln($a + $b) - gammaln($a) - gammaln($b)
                  + $a*log($x) + $b*log(1 - $x) );
    if ( $x < ($a+1)/($a+$b+2) ) {
        return 1 - $bt * _betacf($x, $a, $b) / $a;
    } else {
        return $bt * _betacf(1 - $x, $b, $a) / $b;
    }
}

1;              # Modules must return true.
