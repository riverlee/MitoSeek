# -*- mode: perl; coding: utf-8-unix -*-
#
# Author:      Peter John Acklam
# Time-stamp:  2009-03-30 16:22:42 +02:00
# E-mail:      pjacklam@online.no
# URL:         http://home.online.no/~pjacklam

=pod

=head1 NAME

Math::SpecFun::Gamma - gamma, log(gamma), incomplete gamma etc.

=head1 SYNOPSIS

    use Math::SpecFun::Gamma qw(gamma gammaln gammainc gammaincc);

Imports all the routines explicitly.  Use a subset of the list for the
routines you want.

    use Math::SpecFun::Gamma qw(:all);

Imports all the routines, as well.

=head1 DESCRIPTION

This module implements the gamma function, C<gamma>, the natural
logarithm of the gamma function, C<gammaln>, the incomplete gamma
function, C<gammainc>, and its complement C<gammaincc>.

For references and details about the algorithms, see the comments inside
this module.

=head1 FUNCTIONS

=over 8

=item gamma EXPR

=item gamma

Returns the gamma function evaluated at EXPR.  If EXPR is omitted, C<$_>
is used.  The gamma function is defined as

    gamma(x) = integral from 0 to infinity of t^(x-1) exp(-t) dt

The following reflection formula holds

    gamma(x) = PI / ( sin(PI*x) * gamma(1-x) )

Note that the gamma function has poles at all non-positive integer
values of x.

Here is a factorial function implemented in terms of the gamma
function

    use Math::SpecFun::Gamma qw(gamma);

    sub factorial ($) {
        gamma 1 + shift;
    }

=item gammaln EXPR

=item gammaln

Returns the natural logarithm of the gamma function evaluated at EXPR.
If EXPR is omitted, C<$_> is used.

Here is a function that implements the binomial function in terms of
C<gammaln>

    use Math::SpecFun::Gamma qw(gammaln);

    sub binomial ($$) {
        my ($n, $k) = @_;
        int 0.5 + exp gammaln($n+1) - gammaln($k+1) - gammaln($n-$k+1);
    }

which is, for large arguments, somewhat faster than the more
"traditional" approach

    sub binomial ($$) {
        my ($n, $k) = @_;
        $k = $n - $k if $k > $n / 2;
        return 1 unless $k > 0;
        my $r = 1;
        do { $r *= $n--/$k } while --$k;
        int 0.5 + sprintf "%.f", $r;
    }

=item gammainc X, A

Returns the incomplete gamma function.  The incomplete gamma function
is defined as:

    gammainc(x, a) = 1 / gamma(a) *
        integral from 0 to x of t^(a-1) exp(-t) dt

For any a >= 0, as x approaches infinity, gammainc(x,a) approaches 1.
For small x and a, gammainc(x,a) ~= x^a, so gammainc(0,0) = 1.

=item gammaincc X, A

Returns the complement of the incomplete gamma function.  The
incomplete gamma function is defined as:

    gammaincc(x, a) = 1 / gamma(a) *
        integral from x to infinity of t^(a-1) exp(-t) dt

    gammaincc(x, a) = 1 - gammainc(x, a)

=item factorial EXPR

=item factorial

Returns the factorial function evaluated at EXPR.  If EXPR is omitted, C<$_>
is used.  The factorial function is defined as

    factorial(x) = x! = 1*2*...*x

=back

=head1 HISTORY

=over 4

=item Version 0.03

Added the incomplete gamma function and its complement.

=item Version 0.02

Added the factorial function.  Extended the user functions.

=item Version 0.01

First release.

=back

=head1 AUTHORS

Perl translation by Peter John Acklam <pjacklam@online.no>

FORTRAN code for C<gamma> by W. J. Cody and L. Stoltz, Argonne
National Laboratory, October 12, 1989.  FORTRAN code can be found at
http://www.netlib.org/specfun/algama

FORTRAN code for C<gammaln> by W. J. Cody and L. Stoltz, Argonne
National Laboratory, June 16, 1988.  FORTRAN code can be found at
http://www.netlib.org/specfun/gamma

C code for C<gammainc> and C<gammaincc> found in Numerical Recipes in C.

=head1 COPYRIGHT

Copyright 1999-2000 Peter John Acklam.
This program is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

=cut

package Math::SpecFun::Gamma;
require 5.000;
require Exporter;

use strict;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
use Carp;

$VERSION = '0.03';

@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(gamma gammaln gammainc gammaincc factorial);
%EXPORT_TAGS = (all => [ @EXPORT_OK ]);

use constant PI => 4 * atan2 1, 1;

########################################################################
## Internal functions.
########################################################################

##
## This is a straightforward Perl translation of the FORTRAN routine
## NETLIB/SPECFUN/GAMMA.F.  It should not be called by users.  Changes
## to the interface and handling of exceptional cases etc. should
## preferably be done in the user functions rather than in this
## routine.
##

sub _gamma ($) {
#----------------------------------------------------------------------
#
# This routine calculates the GAMMA function for a real argument X.
#   Computation is based on an algorithm outlined in reference 1.
#   The program uses rational functions that approximate the GAMMA
#   function to at least 20 significant decimal digits.  Coefficients
#   for the approximation over the interval (1,2) are unpublished.
#   Those for the approximation for X >= 12 are from reference 2.
#   The accuracy achieved depends on the arithmetic system, the
#   compiler, the intrinsic functions, and proper selection of the
#   machine-dependent constants.
#
#
#*******************************************************************
#*******************************************************************
#
# Explanation of machine-dependent constants
#
# beta   - radix for the floating-point representation
# maxexp - the smallest positive power of beta that overflows
# XBIG   - the largest argument for which GAMMA(X) is representable
#          in the machine, i.e., the solution to the equation
#                  GAMMA(XBIG) = beta**maxexp
# XINF   - the largest machine representable floating-point number;
#          approximately beta**maxexp
# EPS    - the smallest positive floating-point number such that
#          1.0+EPS .GT. 1.0
# XMININ - the smallest positive floating-point number such that
#          1/XMININ is machine representable
#
#     Approximate values for some important machines are:
#
#                            beta       maxexp        XBIG
#
# CRAY-1         (S.P.)        2         8191        966.961
# Cyber 180/855
#   under NOS    (S.P.)        2         1070        177.803
# IEEE (IBM/XT,
#   SUN, etc.)   (S.P.)        2          128        35.040
# IEEE (IBM/XT,
#   SUN, etc.)   (D.P.)        2         1024        171.624
# IBM 3033       (D.P.)       16           63        57.574
# VAX D-Format   (D.P.)        2          127        34.844
# VAX G-Format   (D.P.)        2         1023        171.489
#
#                            XINF         EPS        XMININ
#
# CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
# Cyber 180/855
#   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
# IEEE (IBM/XT,
#   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
# IEEE (IBM/XT,
#   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
# IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
# VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
# VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
#
#*******************************************************************
#*******************************************************************
#
# Error returns
#
#  The program returns the value XINF for singularities or
#     when overflow would occur.  The computation is believed
#     to be free of underflow and overflow.
#
# References: "An Overview of Software Development for Special
#              Functions", W. J. Cody, Lecture Notes in Mathematics,
#              506, Numerical Analysis Dundee, 1975, G. A. Watson
#              (ed.), Springer Verlag, Berlin, 1976.
#
#              Computer Approximations, Hart, Et. Al., Wiley and
#              sons, New York, 1968.
#
#  Latest modification: October 12, 1989
#
#  Authors: W. J. Cody and L. Stoltz
#           Applied Mathematics Division
#           Argonne National Laboratory
#           Argonne, IL 60439
#
#----------------------------------------------------------------------
    my $x = $_[0];
    my ($y1, $res, $z, $i);
#----------------------------------------------------------------------
#  Mathematical constants
#----------------------------------------------------------------------
    my $sqrtpi = 0.9189385332046727417803297e0;     # log(sqrt(2*pi))
    my $pi     = 3.1415926535897932384626434e0;
#----------------------------------------------------------------------
#  Machine dependent parameters
#----------------------------------------------------------------------
    my $xbig   = 171.624e0;
    my $xminin = 2.23e-308;
    my $eps    = 2.22e-16;
    my $xinf   = 1.79e308;
#----------------------------------------------------------------------
#  Numerator and denominator coefficients for rational minimax
#     approximation over (1,2).
#----------------------------------------------------------------------
    my @P = (-1.71618513886549492533811e+0, 2.47656508055759199108314e+1,
             -3.79804256470945635097577e+2, 6.29331155312818442661052e+2,
             8.66966202790413211295064e+2, -3.14512729688483675254357e+4,
             -3.61444134186911729807069e+4, 6.64561438202405440627855e+4);
    my @Q = (-3.08402300119738975254353e+1, 3.15350626979604161529144e+2,
             -1.01515636749021914166146e+3, -3.10777167157231109440444e+3,
             2.25381184209801510330112e+4, 4.75584627752788110767815e+3,
             -1.34659959864969306392456e+5, -1.15132259675553483497211e+5);
#----------------------------------------------------------------------
#  Coefficients for minimax approximation over (12, INF).
#----------------------------------------------------------------------
    my @C = (-1.910444077728e-03, 8.4171387781295e-04,
             -5.952379913043012e-04, 7.93650793500350248e-04,
             -2.777777777777681622553e-03, 8.333333333333333331554247e-02,
             5.7083835261e-03);
#----------------------------------------------------------------------
#  Statement functions for conversion between integer and float
#----------------------------------------------------------------------
    my $parity = 0;                      # false
    my $fact = 1;
    my $n = 0;
    my $y = $x;
    if ($y <= 0) {
#----------------------------------------------------------------------
#  Argument is negative
#----------------------------------------------------------------------
        $y = -$x;
        $y1 = int($y);
        $res = $y - $y1;
        if ($res != 0) {
            if ($y1 != int($y1*0.5)*2) { $parity = 1 }
            $fact = -$pi / sin($pi*$res);
            $y = $y + 1;
        } else {
            $res = $xinf;
            return $res;
        }
    }
#----------------------------------------------------------------------
#  Argument is positive
#----------------------------------------------------------------------
    if ($y < $eps) {
#----------------------------------------------------------------------
#  Argument <= EPS
#----------------------------------------------------------------------
        if ($y >= $xminin) {
            $res = 1 / $y;
        } else {
            $res = $xinf;
            return $res;
        }
    } elsif ($y < 12) {
        $y1 = $y;
        if ($y < 1) {
#----------------------------------------------------------------------
#  0.0 <= argument <= 1.0
#----------------------------------------------------------------------
            $z = $y;
            $y = $y + 1;
        } else {
#----------------------------------------------------------------------
#  1.0 <= argument <= 12.0, reduce argument if necessary
#----------------------------------------------------------------------
            $n = int($y) - 1;
            $y = $y - $n;
            $z = $y - 1;
        }
#----------------------------------------------------------------------
#  Evaluate approximation for 1.0 <= argument <= 2.0
#----------------------------------------------------------------------
        my $xnum = 0;
        my $xden = 1;
        foreach my $i (0 .. 7) {
            $xnum = ($xnum + $P[$i]) * $z;
            $xden = $xden * $z + $Q[$i];
        }
        $res = $xnum / $xden + 1;
        if ($y1 < $y) {
#----------------------------------------------------------------------
#  Adjust result for case  0.0 <= argument <= 1.0
#----------------------------------------------------------------------
            $res = $res / $y1;
        } elsif ($y1 > $y) {
#----------------------------------------------------------------------
#  Adjust result for case  2.0 <= argument <= 12.0
#----------------------------------------------------------------------
            foreach $i ( 1 .. $n ) {
                $res = $res * $y;
                $y = $y + 1;
            }
        }
    } else {
#----------------------------------------------------------------------
#  Evaluate for argument >= 12.0,
#----------------------------------------------------------------------
        if ($y <= $xbig) {
            my $ysq = $y * $y;
            my $sum = $C[6];
            foreach my $i (0 .. 5) {
                $sum = $sum / $ysq + $C[$i];
            }
            $sum = $sum/$y - $y + $sqrtpi;
            $sum = $sum + ($y-0.5)*log($y);
            $res = exp($sum);
        } else {
            $res = $xinf;
            return $res;
        }
    }
#----------------------------------------------------------------------
#  Final adjustments and return
#----------------------------------------------------------------------
    if ($parity) { $res = -$res }
    if ($fact != 1) { $res = $fact / $res }
    return $res;
}

##
## This is a straightforward Perl translation of the FORTRAN routine
## NETLIB/SPECFUN/ALGAMA.F.  It should not be called by users.
## Changes to the interface and handling of exceptional cases
## etc. should preferably be done in the user functions rather than in
## this routine.
##

sub _gammaln ($) {
#----------------------------------------------------------------------
#
# This routine calculates the LOG(GAMMA) function for a positive real
#   argument X.  Computation is based on an algorithm outlined in
#   references 1 and 2.  The program uses rational functions that
#   theoretically approximate LOG(GAMMA) to at least 18 significant
#   decimal digits.  The approximation for X > 12 is from reference
#   3, while approximations for X < 12.0 are similar to those in
#   reference 1, but are unpublished.  The accuracy achieved depends
#   on the arithmetic system, the compiler, the intrinsic functions,
#   and proper selection of the machine-dependent constants.
#
#
#*********************************************************************
#*********************************************************************
#
# Explanation of machine-dependent constants
#
# beta   - radix for the floating-point representation
# maxexp - the smallest positive power of beta that overflows
# XBIG   - largest argument for which LN(GAMMA(X)) is representable
#          in the machine, i.e., the solution to the equation
#                  LN(GAMMA(XBIG)) = beta**maxexp
# XINF   - largest machine representable floating-point number;
#          approximately beta**maxexp.
# EPS    - The smallest positive floating-point number such that
#          1.0+EPS .GT. 1.0
# FRTBIG - Rough estimate of the fourth root of XBIG
#
#
#     Approximate values for some important machines are:
#
#                           beta      maxexp         XBIG
#
# CRAY-1        (S.P.)        2        8191       9.62E+2461
# Cyber 180/855
#   under NOS   (S.P.)        2        1070       1.72E+319
# IEEE (IBM/XT,
#   SUN, etc.)  (S.P.)        2         128       4.08E+36
# IEEE (IBM/XT,
#   SUN, etc.)  (D.P.)        2        1024       2.55D+305
# IBM 3033      (D.P.)       16          63       4.29D+73
# VAX D-Format  (D.P.)        2         127       2.05D+36
# VAX G-Format  (D.P.)        2        1023       1.28D+305
#
#                           XINF        EPS        FRTBIG
#
# CRAY-1        (S.P.)   5.45E+2465   7.11E-15    3.13E+615
# Cyber 180/855
#   under NOS   (S.P.)   1.26E+322    3.55E-15    6.44E+79
# IEEE (IBM/XT,
#   SUN, etc.)  (S.P.)   3.40E+38     1.19E-7     1.42E+9
# IEEE (IBM/XT,
#   SUN, etc.)  (D.P.)   1.79D+308    2.22D-16    2.25D+76
# IBM 3033      (D.P.)   7.23D+75     2.22D-16    2.56D+18
# VAX D-Format  (D.P.)   1.70D+38     1.39D-17    1.20D+9
# VAX G-Format  (D.P.)   8.98D+307    1.11D-16    1.89D+76
#
#**************************************************************
#**************************************************************
#
# Error returns
#
#  The program returns the value XINF for X <= 0.0 or when
#     overflow would occur.  The computation is believed to
#     be free of underflow and overflow.
#
# References:
#
#  1) W. J. Cody and K. E. Hillstrom, 'Chebyshev Approximations for
#     the Natural Logarithm of the Gamma Function,' Math. Comp. 21,
#     1967, pp. 198-203.
#
#  2) K. E. Hillstrom, ANL/AMD Program ANLC366S, DGAMMA/DLGAMA, May,
#     1969.
#
#  3) Hart, Et. Al., Computer Approximations, Wiley and sons, New
#     York, 1968.
#
#  Authors: W. J. Cody and L. Stoltz
#           Argonne National Laboratory
#
#  Latest modification: June 16, 1988
#
#----------------------------------------------------------------------
    my $x = $_[0];
    my ($res, $corr, $xnum, $xden, $xm1, $xm2, $xm4);
#----------------------------------------------------------------------
#  Mathematical constants
#----------------------------------------------------------------------
    my $pnt68  = 0.6796875e0;
    my $sqrtpi = 0.9189385332046727417803297e0;
#----------------------------------------------------------------------
#  Machine dependent parameters
#----------------------------------------------------------------------
    my $xbig   = 2.55e305;
    my $xinf   = 1.79e308;
    my $eps    = 2.22e-16;
    my $frtbig = 2.25e76;
#----------------------------------------------------------------------
#  Numerator and denominator coefficients for rational minimax
#     approximation over (0.5,1.5).
#----------------------------------------------------------------------
    my $D1 = -5.772156649015328605195174e-1;
    my @P1 = (4.945235359296727046734888e0, 2.018112620856775083915565e2,
              2.290838373831346393026739e3, 1.131967205903380828685045e4,
              2.855724635671635335736389e4, 3.848496228443793359990269e4,
              2.637748787624195437963534e4, 7.225813979700288197698961e3);
    my @Q1 = (6.748212550303777196073036e1, 1.113332393857199323513008e3,
              7.738757056935398733233834e3, 2.763987074403340708898585e4,
              5.499310206226157329794414e4, 6.161122180066002127833352e4,
              3.635127591501940507276287e4, 8.785536302431013170870835e3);
#----------------------------------------------------------------------
#  Numerator and denominator coefficients for rational minimax
#     Approximation over (1.5,4.0).
#----------------------------------------------------------------------
    my $D2 = 4.227843350984671393993777e-1;
    my @P2 = (4.974607845568932035012064e0, 5.424138599891070494101986e2,
              1.550693864978364947665077e4, 1.847932904445632425417223e5,
              1.088204769468828767498470e6, 3.338152967987029735917223e6,
              5.106661678927352456275255e6, 3.074109054850539556250927e6);
    my @Q2 = (1.830328399370592604055942e2, 7.765049321445005871323047e3,
              1.331903827966074194402448e5, 1.136705821321969608938755e6,
              5.267964117437946917577538e6, 1.346701454311101692290052e7,
              1.782736530353274213975932e7, 9.533095591844353613395747e6);
#----------------------------------------------------------------------
#  Numerator and denominator coefficients for rational minimax
#     Approximation over (4.0,12.0).
#----------------------------------------------------------------------
    my $D4 = 1.791759469228055000094023e0;
    my @P4 = (1.474502166059939948905062e4, 2.426813369486704502836312e6,
              1.214755574045093227939592e8, 2.663432449630976949898078e9,
              2.940378956634553899906876e10, 1.702665737765398868392998e11,
              4.926125793377430887588120e11, 5.606251856223951465078242e11);
    my @Q4 = (2.690530175870899333379843e3, 6.393885654300092398984238e5,
              4.135599930241388052042842e7, 1.120872109616147941376570e9,
              1.488613728678813811542398e10, 1.016803586272438228077304e11,
              3.417476345507377132798597e11, 4.463158187419713286462081e11);
#----------------------------------------------------------------------
#  Coefficients for minimax approximation over (12, INF).
#----------------------------------------------------------------------
    my @C = (-1.910444077728e-03, 8.4171387781295e-04,
             -5.952379913043012e-04, 7.93650793500350248e-04,
             -2.777777777777681622553e-03, 8.333333333333333331554247e-02,
             5.7083835261e-03);
#----------------------------------------------------------------------
    my $y = $x;
    if (($y > 0) && ($y <= $xbig)) {
        if ($y <= $eps) {
            $res = -log($y);
        } elsif ($y <= 1.5) {
#----------------------------------------------------------------------
#  EPS < X <= 1.5
#----------------------------------------------------------------------
            if ($y < $pnt68) {
                $corr = -log($y);
                $xm1 = $y
            } else {
                $corr = 0;
                $xm1 = ($y - 0.5) - 0.5;
            }
            if (($y <= 0.5) || ($y >= $pnt68)) {
                $xden = 1;
                $xnum = 0;
                foreach my $i (0 .. 7) {
                    $xnum = $xnum*$xm1 + $P1[$i];
                    $xden = $xden*$xm1 + $Q1[$i];
                }
                $res = $corr + ($xm1 * ($D1 + $xm1*($xnum/$xden)));
            } else {
                $xm2 = ($y - 0.5) - 0.5;
                $xden = 1;
                $xnum = 0;
                foreach my $i (0 .. 7) {
                    $xnum = $xnum*$xm2 + $P2[$i];
                    $xden = $xden*$xm2 + $Q2[$i];
                }
                $res = $corr + $xm2 * ($D2 + $xm2*($xnum/$xden));
            }
        } elsif ($y <= 4) {
#----------------------------------------------------------------------
#  1.5 < X <= 4.0
#----------------------------------------------------------------------
            $xm2 = $y - 2;
            $xden = 1;
            $xnum = 0;
            foreach my $i (0 .. 7) {
                $xnum = $xnum*$xm2 + $P2[$i];
                $xden = $xden*$xm2 + $Q2[$i];
            }
            $res = $xm2 * ($D2 + $xm2*($xnum/$xden));
        } elsif ($y <= 12) {
#----------------------------------------------------------------------
#  4.0 < X <= 12.0
#----------------------------------------------------------------------
            $xm4 = $y - 4;
            $xden = -1;
            $xnum = 0;
            foreach my $i (0 .. 7) {
                $xnum = $xnum*$xm4 + $P4[$i];
                $xden = $xden*$xm4 + $Q4[$i];
            }
            $res = $D4 + $xm4*($xnum/$xden);
        } else {
#----------------------------------------------------------------------
#  Evaluate for argument .GE. 12.0,
#----------------------------------------------------------------------
            $res = 0;
            if ($y <= $frtbig) {
                $res = $C[6];
                my $ysq = $y * $y;
                foreach my $i (0 .. 5) {
                    $res = $res / $ysq + $C[$i];
                }
            }
            $res = $res/$y;
            $corr = log($y);
            $res = $res + $sqrtpi - 0.5*$corr;
            $res = $res + $y*($corr-1);
        }
    } else {
#----------------------------------------------------------------------
#  Return for bad arguments
#----------------------------------------------------------------------
        $res = $xinf;
    }
#----------------------------------------------------------------------
#  Final adjustments and return
#----------------------------------------------------------------------
    return $res;
}

##
## Returns the incomplete gamma function P(a; x) evaluated by its
## series representation.
##

sub _gser ($$) {
    my ($x, $a) = @_;

    my $itmax = 100;
    my $eps = 3e-7;
    my $gln = _gammaln $a;

    if ($x <= 0) {
        croak "x less than 0 in routine _gser" if $x < 0;
        return 0;
    } else {
        my $ap = $a;
        my $del = my $sum = 1/$a;
        for ( my $n = 1 ; $n <= $itmax ; $n++ ) {
            ++$ap;
            $del *= $x/$ap;
            $sum += $del;
            if (abs($del) < abs($sum)*$eps) {
                return $sum * exp( -$x + $a*log($x) - $gln );
            }
        }
        croak "a too large, ITMAX too small in routine _gser";
    }
}

##
## Returns the incomplete gamma function Q(a; x) evaluated by its
## continued fraction representation.
##

sub _gcf ($$) {
    my ($x, $a) = @_;

    my $itmax = 100;        # Maximum allowed number of iterations.
    my $eps = 3e-7;         # Relative accuracy.
    my $fpmin = 1e-30;      # Number near the smallest representable floating-point number.
    my $gln = gammaln $a;

    my $b = $x+1-$a;              # Set up for evaluating continued fraction
    my $c = 1/$fpmin;             # by modified Lentz's method ( x 5.2)
    my $d = 1/$b;                 # with b0 = 0.
    my $h = $d;
    my $i;
    for ( $i = 1 ; $i <= $itmax ; $i++ ) {      # Iterate to convergence.
        my $an = -$i*($i-$a);
        $b += 2;
        $d = $an*$d + $b;
        if (abs($d) < $fpmin) { $d = $fpmin };
        $c = $b + $an/$c;
        if (abs($c) < $fpmin) { $c = $fpmin };
        $d = 1/$d;
        my $del = $d*$c;
        $h *= $del;
        if (abs($del-1) < $eps) { last };
    }
    croak "a too large, ITMAX too small in _gcf" if $i > $itmax;
    return exp( -$x + $a*log($x) - $gln )*$h;   # Put factors in front.
}

########################################################################
## User functions.
########################################################################

sub gamma {
    # Get and check input arguments.
    croak "usage: gamma() or gamma(x)\n" if @_ > 1;
    my $x = @_ ? shift : $_;

    croak sprintf "gamma of %g is undefined\n", $x
      if $x == int $x && $x <= 0;

    return _gamma $x;
}

sub gammaln {
    # Get and check input arguments.
    croak "usage: gammaln() or gammaln(x)\n" if @_ > 1;
    my $x = @_ ? shift : $_;

    if ( $x <= 0 ) {
        croak sprintf "gammaln of %g is undefined\n", $x if $x == int $x;
        return log(PI) - log(sin(PI*$x)) - gammaln(1-$x);
    }
    return _gammaln $x;
}

sub gammainc ($$) {
    # Returns the incomplete gamma function P(a; x).
    my ($x, $a) = @_;
    croak "First argument must be non-negative in gammainc" if $x < 0;
    croak "Second argument must be positive in gammainc"    if $a <= 0;
    if ($x < ($a+1)) {              # Use the series representation.
        return _gser($x, $a);
    } else {                        # Use the continued fraction representation
        return 1 - _gcf($x, $a);    # and take its complement.
    }
}

sub gammaincc ($$) {
    # Returns the incomplete gamma function Q(a; x) = 1 - P(a; x).
    my ($x, $a) = @_;
    croak "First argument must be non-negative in gammaincc" if $x < 0;
    croak "Second argument must be positive in gammaincc"    if $a <= 0;
    if ($x < ($a+1)) {              # Use the series representation
        return 1 - _gser($x, $a);   #   and take its complement.
    } else {                        # Use the continued fraction representation.
        return _gcf($x, $a);
    }
}

##
## This function seemed to fit into the company of the functions above.
##

sub factorial {
    # Get and check input arguments.
    croak "usage: factorial() or factorial(x)\n" if @_ > 1;
    my $x = @_ ? shift : $_;
    croak sprintf "factorial of %g is undefined\n", $x
      if $x != int $x || $x < 0;

    return 1 unless $x;
    my $y = $x;
    $y *= $x while --$x;
    return $y;
}

1;              # Modules must return true.
