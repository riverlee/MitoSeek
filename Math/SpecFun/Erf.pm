# -*- mode: perl; coding: utf-8-unix -*-
#
# Author:      Peter John Acklam
# Time-stamp:  2009-03-30 16:22:42 +02:00
# E-mail:      pjacklam@online.no
# URL:         http://home.online.no/~pjacklam

=pod

=head1 NAME

Math::SpecFun::Erf - error and scaled and unscaled complementary error
functions and their inverses

=head1 SYNOPSIS

    use Math::SpecFun::Erf qw(erf erfc erfcx erfinv erfcinv erfcxinv);

Imports all the routines explicitly.  Use a subset of the list for the
routines you want.

    use Math::SpecFun::Erf qw(:all);

Imports all the routines, as well.

=head1 DESCRIPTION

This module implements the error function, C<erf>, and its inverse
C<erfinv>, the complementary error function, C<erfc>, and its inverse
C<erfcinv>, and the scaled complementary error function, C<erfcx>, and its
inverse C<erfcxinv>.

For references and details about the algorithms, see the comments inside
this module.

=head1 FUNCTIONS

=over 8

=item erf EXPR

=item erf

Returns the error function evaluated at EXPR.  If EXPR is omitted, C<$_> is
used.  The error function is

    erf(x) = 2/sqrt(PI) * integral from 0 to x of exp(-t*t) dt

=item erfinv EXPR

=item erfinv

Returns the inverse of the error function evaluated at EXPR.  If EXPR is
omitted, C<$_> is used.

=item erfc EXPR

=item erfc

Returns the complementary error function evaluated at EXPR.  If EXPR is
omitted, C<$_> is used.  The complementary error function is

    erfc(x) = 2/sqrt(PI) * integral from x to infinity of exp(-t*t) dt
            = 1 - erf(x)

Here is a function returning the lower tail probability of the standard
normal distribution function

    use Math::SpecFun::Erf qw(erfc);

    sub ltpnorm ($) {
        erfc( - $_[0] / sqrt(2) )/2;
    }

=item erfcinv EXPR

=item erfcinv

Returns the inverse complementary error function evaluated at EXPR.  If EXPR
is omitted, C<$_> is used.

Here is a function returning the lower tail quantile of the standard normal
distribution function

    use Math::SpecFun::Erf qw(erfcinv);

    sub ltqnorm ($) {
        -sqrt(2) * erfcinv( 2 * $_[0] );
    }

=item erfcx EXPR

=item erfcx

Returns the scaled complementary error function evaluated at EXPR.  If EXPR
is omitted, C<$_> is used.  The scaled complementary error function is

    erfcx(x) = exp(x*x) * erfc(x)

=item erfcxinv EXPR

=item erfcxinv

Returns the inverse scaled complementary error function evaluated at EXPR.
If EXPR is omitted, C<$_> is used.

=back

=head1 HISTORY

=over 4

=item Version 0.03

Added the inverse functions.

=item Version 0.02

Minor code tweaking.

=item Version 0.01

First release.

=back

=head1 AUTHOR

Perl translation by Peter John Acklam E<lt>pjacklam@online.noE<gt>

FORTRAN code by W. J. Cody, Argonne National Laboratory, March 19, 1990.
FORTRAN code can be found at http://www.netlib.org/specfun/erf

=head1 COPYRIGHT

Copyright 1999-2000 Peter John Acklam.
This program is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

=cut

package Math::SpecFun::Erf;
require 5.000;
require Exporter;

use strict;
use vars qw($VERSION @ISA @EXPORT_OK %EXPORT_TAGS);

$VERSION = '0.02';
@ISA = qw(Exporter);
@EXPORT_OK = qw(erf erfc erfcx erfinv erfcinv erfcxinv);
%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

########################################################################
## Internal functions.
########################################################################

sub calerf {
    my ($arg, $result, $jint) = @_;
    local $[ = 1;
#------------------------------------------------------------------
#
# This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
#   for a real argument  x.  It contains three FUNCTION type
#   subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX),
#   and one SUBROUTINE type subprogram, CALERF.  The calling
#   statements for the primary entries are:
#
#                   Y=ERF(X)     (or   Y=DERF(X)),
#
#                   Y=ERFC(X)    (or   Y=DERFC(X)),
#   and
#                   Y=ERFCX(X)   (or   Y=DERFCX(X)).
#
#   The routine  CALERF  is intended for internal packet use only,
#   all computations within the packet being concentrated in this
#   routine.  The function subprograms invoke  CALERF  with the
#   statement
#
#          CALL CALERF(ARG,RESULT,JINT)
#
#   where the parameter usage is as follows
#
#      Function                     Parameters for CALERF
#       call              ARG                  Result          JINT
#
#     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
#     ERFC(ARG)     ABS(ARG) < XBIG           ERFC(ARG)         1
#     ERFCX(ARG)    XNEG < ARG < XMAX         ERFCX(ARG)        2
#
#   The main computation evaluates near-minimax approximations
#   from "Rational Chebyshev approximations for the error function"
#   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
#   transportable program uses rational functions that theoretically
#   approximate  erf(x)  and  erfc(x)  to at least 18 significant
#   decimal digits.  The accuracy achieved depends on the arithmetic
#   system, the compiler, the intrinsic functions, and proper
#   selection of the machine-dependent constants.
#
#*******************************************************************
#*******************************************************************
#
# Explanation of machine-dependent constants
#
#   XMIN   = the smallest positive floating-point number.
#   XINF   = the largest positive finite floating-point number.
#   XNEG   = the largest negative argument acceptable to ERFCX;
#            the negative of the solution to the equation
#            2*exp(x*x) = XINF.
#   XSMALL = argument below which erf(x) may be represented by
#            2*x/sqrt(pi)  and above which  x*x  will not underflow.
#            A conservative value is the largest machine number X
#            such that   1.0 + X = 1.0   to machine precision.
#   XBIG   = largest argument acceptable to ERFC;  solution to
#            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
#            W(x) = exp(-x*x)/[x*sqrt(pi)].
#   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
#            machine precision.  A conservative value is
#            1/[2*sqrt(XSMALL)]
#   XMAX   = largest acceptable argument to ERFCX; the minimum
#            of XINF and 1/[sqrt(pi)*XMIN].
#
#   Approximate values for some important machines are:
#
#                          XMIN       XINF        XNEG     XSMALL
#
#  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
#  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
#  IEEE (IBM/XT,
#    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
#  IEEE (IBM/XT,
#    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
#  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
#  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
#  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
#  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
#
#
#                          XBIG       XHUGE       XMAX
#
#  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
#  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
#  IEEE (IBM/XT,
#    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
#  IEEE (IBM/XT,
#    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
#  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
#  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
#  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
#  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
#
#*******************************************************************
#*******************************************************************
#
# Error returns
#
#  The program returns  ERFC = 0      for  ARG >= XBIG;
#
#                       ERFCX = XINF  for  ARG < XNEG;
#      and
#                       ERFCX = 0     for  ARG >= XMAX.
#
#
# Intrinsic functions required are:
#
#     ABS, AINT, EXP
#
#
#  Author: W. J. Cody
#          Mathematics and Computer Science Division
#          Argonne National Laboratory
#          Argonne, IL 60439
#
#  Latest modification: March 19, 1990
#
# Translation to Perl by Peter John Acklam, December 3, 1999
#
#------------------------------------------------------------------
    my ($i);
    my ($x, $del, $xden, $xnum, $y, $ysq);
#------------------------------------------------------------------
#  Mathematical constants
#------------------------------------------------------------------
    my ($four, $one, $half, $two, $zero) = (4, 1, 0.5, 2, 0);
    my $sqrpi = 5.6418958354775628695e-1;
    my $thresh = 0.46875;
    my $sixten = 16;
#------------------------------------------------------------------
#  Machine-dependent constants
#------------------------------------------------------------------
    my ($xinf, $xneg, $xsmall) = (1.79e308, -26.628, 1.11e-16);
    my ($xbig, $xhuge, $xmax) = (26.543, 6.71e7, 2.53e307);
#------------------------------------------------------------------
#  Coefficients for approximation to  erf  in first interval
#------------------------------------------------------------------
    my @a = (3.16112374387056560e00, 1.13864154151050156e02,
             3.77485237685302021e02, 3.20937758913846947e03,
             1.85777706184603153e-1);
    my @b = (2.36012909523441209e01, 2.44024637934444173e02,
             1.28261652607737228e03, 2.84423683343917062e03);
#------------------------------------------------------------------
#  Coefficients for approximation to  erfc  in second interval
#------------------------------------------------------------------
    my @c = (5.64188496988670089e-1, 8.88314979438837594e0,
             6.61191906371416295e01, 2.98635138197400131e02,
             8.81952221241769090e02, 1.71204761263407058e03,
             2.05107837782607147e03, 1.23033935479799725e03,
             2.15311535474403846e-8);
    my @d = (1.57449261107098347e01, 1.17693950891312499e02,
             5.37181101862009858e02, 1.62138957456669019e03,
             3.29079923573345963e03, 4.36261909014324716e03,
             3.43936767414372164e03, 1.23033935480374942e03);
#------------------------------------------------------------------
#  Coefficients for approximation to  erfc  in third interval
#------------------------------------------------------------------
    my @p = (3.05326634961232344e-1, 3.60344899949804439e-1,
             1.25781726111229246e-1, 1.60837851487422766e-2,
             6.58749161529837803e-4, 1.63153871373020978e-2);
    my @q = (2.56852019228982242e00, 1.87295284992346047e00,
             5.27905102951428412e-1, 6.05183413124413191e-2,
             2.33520497626869185e-3);
#------------------------------------------------------------------
    $x = $arg;
    $y = abs($x);
    if ($y <= $thresh) {
#------------------------------------------------------------------
#  Evaluate  erf  for  |X| <= 0.46875
#------------------------------------------------------------------
        $ysq = $zero;
        if ($y > $xsmall) { $ysq = $y * $y }
        $xnum = $a[5]*$ysq;
        $xden = $ysq;
        for (my $i = 1 ; $i <= 3 ; ++$i) {
            $xnum = ($xnum + $a[$i]) * $ysq;
            $xden = ($xden + $b[$i]) * $ysq;
        }
        $$result = $x * ($xnum + $a[4]) / ($xden + $b[4]);
        if ($jint != 0) { $$result = $one - $$result }
        if ($jint == 2) { $$result = exp($ysq) * $$result }
        goto x800;
#------------------------------------------------------------------
#  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
#------------------------------------------------------------------
    } elsif ($y <= $four) {
        $xnum = $c[9]*$y;
        $xden = $y;
        for (my $i = 1 ; $i <= 7 ; ++$i) {
            $xnum = ($xnum + $c[$i]) * $y;
            $xden = ($xden + $d[$i]) * $y;
        }
        $$result = ($xnum + $c[8]) / ($xden + $d[8]);
        if ($jint != 2) {
            $ysq = int($y*$sixten)/$sixten;
            $del = ($y-$ysq)*($y+$ysq);
            $$result = exp(-$ysq*$ysq) * exp(-$del) * $$result;
        }
#------------------------------------------------------------------
#  Evaluate  erfc  for |X| > 4.0
#------------------------------------------------------------------
    } else {
        $$result = $zero;
        if ($y >= $xbig) {
            if (($jint != 2) || ($y >= $xmax)) { goto x300 }
            if ($y >= $xhuge) {
                $$result = $sqrpi / $y;
                goto x300;
            }
        }
        $ysq = $one / ($y * $y);
        $xnum = $p[6]*$ysq;
        $xden = $ysq;
        for (my $i = 1 ; $i <= 4 ; ++$i) {
            $xnum = ($xnum + $p[$i]) * $ysq;
            $xden = ($xden + $q[$i]) * $ysq;
        }
        $$result = $ysq *($xnum + $p[5]) / ($xden + $q[5]);
        $$result = ($sqrpi -  $$result) / $y;
        if ($jint != 2) {
            $ysq = int($y*$sixten)/$sixten;
            $del = ($y-$ysq)*($y+$ysq);
            $$result = exp(-$ysq*$ysq) * exp(-$del) * $$result;
        }
    }
#------------------------------------------------------------------
#  Fix up for negative argument, erf, etc.
#------------------------------------------------------------------
  x300:
    if ($jint == 0) {
        $$result = ($half - $$result) + $half;
        if ($x < $zero) { $$result = -$$result }
    } elsif ($jint == 1) {
        if ($x < $zero) { $$result = $two - $$result }
    } else {
        if ($x < $zero) {
            if ($x < $xneg) {
                $$result = $xinf;
            } else {
                $ysq = int($x*$sixten)/$sixten;
                $del = ($x-$ysq)*($x+$ysq);
                $y = exp($ysq*$ysq) * exp($del);
                $$result = ($y+$y) - $$result;
            }
        }
    }
  x800:
    return 1;
#---------- Last card of CALERF ----------
}

sub erf {
    my $x = @_ ? $_[0] : $_;
#--------------------------------------------------------------------
#
# This subprogram computes approximate values for erf(x).
#   (see comments heading CALERF).
#
#   Author/date: W. J. Cody, January 8, 1985
#
# Translation to Perl by Peter John Acklam, December 3, 1999
#
#--------------------------------------------------------------------
    my $result;
    my $jint = 0;
    calerf($x, \$result, $jint);
    my $erf = $result;
    return $erf;
#---------- Last card of ERF ----------
}

########################################################################
## User functions.
########################################################################

sub erfc {
    my $x = @_ ? $_[0] : $_;
#--------------------------------------------------------------------
#
# This subprogram computes approximate values for erfc(x).
#   (see comments heading CALERF).
#
#   Author/date: W. J. Cody, January 8, 1985
#
# Translation to Perl by Peter John Acklam, December 3, 1999
#
#--------------------------------------------------------------------
    my ($result);
    my $jint = 1;
    calerf($x, \$result, $jint);
    my $erfc = $result;
    return $erfc;
#---------- Last card of ERFC ----------
}

sub erfcx {
    my $x = @_ ? $_[0] : $_;
#------------------------------------------------------------------
#
# This subprogram computes approximate values for exp(x*x) * erfc(x).
#   (see comments heading CALERF).
#
#   Author/date: W. J. Cody, March 30, 1987
#
# Translation to Perl by Peter John Acklam, December 3, 1999
#
#------------------------------------------------------------------
    my ($result);
    my $jint = 2;
    calerf($x, \$result, $jint);
    my $erfcx = $result;
    return $erfcx;
#---------- Last card of ERFCX ----------
}

sub erfinv {
    my $y = @_ ? $_[0] : $_;

    return  0             if $y == 0;
    return  erfcinv(1-$y) if $y >  0.5;
    return -erfcinv(1+$y) if $y < -0.5;

    #
    # Halley's rational 3rd order method:
    #   u <- f(x)/f'(x)
    #   v <- f''(x)/f'(x)
    #   x <- x - u/(1-u*v/2)
    #
    # Here:
    #   f(x) = erf(x) - y
    #   f'(x) = 2/sqrt(pi)*exp(-x*x)
    #   f''(x) = -4/sqrt(pi)*x*exp(-x*x)
    #
    my $x = 0;
    my $dx;
    my $c = .88622692545275801364908374167055;  # sqrt(pi)/2
    my $eps = 5e-15;
    do {
        my $f = erf($x) - $y;
        my $u = $c*$f*exp($x*$x);
        $dx = -$u/(1+$u*$x);
        $x += $dx;
    } until abs($dx/$x) <= $eps;
    return $x;
}

sub erfcinv {
    my $y = @_ ? $_[0] : $_;

    return 0 if $y == 1;

    my $flag = $y > 1;
    $y = 2 - $y if $flag;

    #
    # Halley's rational 3rd order method:
    #   u <- f(x)/f'(x)
    #   v <- f''(x)/f'(x)
    #   x <- x - u/(1-u*v/2)
    #
    # Here:
    #   f(x) = erfc(x) - y
    #   f'(x) = -2/sqrt(pi)*exp(-x*x)
    #   f''(x) = 4/sqrt(pi)*x*exp(-x*x)
    #
    my $x = 0;
    my $dx;
    my $c = -.88622692545275801364908374167055; # sqrt(pi)/2
    my $eps = 5e-15;
    do {
        my $u = $c*(erfcx($x) - $y*exp($x*$x));
        $dx = -$u/(1+$u*$x);
        $x += $dx;
    } until abs($dx/$x) <= $eps;

    return $flag ? -$x : $x;
}

sub erfcxinv {
    my $y = @_ ? $_[0] : $_;

    return 0 if $y == 1;

    #
    # Halley's 3rd order method:
    #   u <- f(x)/f'(x)
    #   v <- f''(x)/f'(x)
    #   x <- x - u/(1-u*v/2)
    #
    # Here:
    #   f(x) = erfcx(x) - y
    #   f'(x) = 2*(x*erfcx(x)-1/sqrt(pi));
    #   f''(x) = (2+4*x*x)*erfcx(x) - 4*x/sqrt(pi);
    #
    my $x = 0;
    my $dx;
    my $c = .56418958354775628694807945156079;  # 1/sqrt(pi)
    my $d = 2.2567583341910251477923178062432;  # 4/sqrt(pi)
    my $eps = 5e-15;
    do {
        my $f = erfcx($x) - $y;
        my $df = 2*($x*erfcx($x)-$c);
        my $ddf = (2+4*$x*$x)*erfcx($x) - $x*$d;
        my $u = $f/$df;
        my $v = $ddf/$df;
        $dx = -$u/(1-$u*$v/2);
        $x += $dx;
    } until abs($dx/$x) <= $eps;
    return $x;
}
