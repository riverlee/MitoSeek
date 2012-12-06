package Circos::Debug;

=pod

=head1 NAME

Circos::Debug - debugging routines for Circos

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

All documentation is in the form of tutorials at L<http://www.circos.ca/tutorials>.

=cut

# -------------------------------------------------------------------

use strict;
use warnings;

use base 'Exporter';
our @EXPORT = qw(get_timer start_timer stop_timer format_timer report_timer printerror printdebug printdebug_group printdumper printwarning printinfo printout debug_or_group debug_group_add debug_group_delete);

use Carp qw( carp confess croak );
use Data::Dumper;
use FindBin;
use Memoize;
use Time::HiRes qw(gettimeofday tv_interval);

use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";

use Circos::Constants;

# -------------------------------------------------------------------
sub list_as_string {
    my $empty_text = Circos::Configuration::fetch_conf("debug_empty_text")     || "_emptylist_";
    return $empty_text if ! @_;
    my $sep        = Circos::Configuration::fetch_conf("debug_word_separator") || " ";
    my $undef_text = Circos::Configuration::fetch_conf("debug_undef_text")     || "_undef_";
    # speedup option - remove field remapping
    my @fields     = Circos::Configuration::fetch_conf("debug_output_tidy") ? remap_fields(@_) : @_;
    return join($sep, map { defined $_ ? $_ : $undef_text } @fields);
}

# -------------------------------------------------------------------
# Opportunity to change the way text is displayed. 
# - reduce precision of floats to 3 decimals
sub remap_fields {
    my @fields = @_;
    my @remapped = ();
    for my $item (map { defined $_ ? split(" ",$_) : $_ } @fields) {
	my $value = $item;
	if(defined $item ) {
	    if($item =~ /^[-+]?\d+\.\d{5,}$/) {
		$value = sprintf("%.3f",$item);
	    }
	}
	push @remapped, $value;
    }
    return @remapped;
}

# -------------------------------------------------------------------
sub errorheader {
    my %args   = @_;
    my $width  = $args{width} || 50;
    my $margin_width = $args{margin} || 2;
    my $delim  = $args{delim} || "*";
    my $text   = $args{text}  || "error";
    my $hdr    = $delim x $width;
    my $margin = " " x $margin_width;
    substr($hdr,(length($hdr) - length($text) - $margin_width*2)/2,length($text)+2*$margin_width) = $margin.$text.$margin;
    return $hdr;
}

# -------------------------------------------------------------------
sub printerror {
    printinfo();
    printinfo(errorheader());
    printinfo();
    printinfo(@_);
    printinfo(errorheader(text=>"debugging"));
    printinfo();
}

# -------------------------------------------------------------------
sub printdebug {
    if(Circos::Configuration::fetch_conf("debug")) {
	printinfo('debug', @_);
    }
}

# -------------------------------------------------------------------
{
    my $t = [gettimeofday];
    sub printdebug_group {
	my ($group,@msg) = @_;
	my $group_label;
	# if the group name is preceeded by !, always print the message
	if($group =~ /^!/) {
	    $group_label = $group;
	} elsif (defined $group && debug_or_group($group)) {
	    $group_label = ref $group eq "ARRAY" ? join($COMMA,@$group) : $group;
	}
	if($group_label) {
	    printinfo('debuggroup',$group_label,sprintf("%.2fs",tv_interval($t)),@msg);
	} 
    }
}

sub format_timer {
    my $t = shift;
    return sprintf("%.3f s",tv_interval($t));
}

{
    my $timers    = {};
    my @lasttimer = ();
    sub get_timer {
	my $timer = shift;
	return $timers->{$timer}{elapsed};
    }
    sub start_timer {
	my $timer = shift;
	return if $timers->{$timer}{start};
	$timers->{$timer}{start}   = [gettimeofday];
	$timers->{$timer}{elapsed} ||= 0;
	push @lasttimer, $timer;
    }
    sub stop_timer {
	my $timer = shift;
	if(! defined $timer) {
	    $timer = pop @lasttimer;
	}
	return unless defined $timers->{$timer}{start};
	@lasttimer = grep($_ ne $timer, @lasttimer);
	$timers->{$timer}{elapsed} += tv_interval( $timers->{$timer}{start} );
	delete $timers->{$timer}{start};
    }
    sub report_timer {
	my $timer = shift;
	my @timers;
	if(! defined $timer) {
	    @timers = sort keys %$timers;
	} else {
	    @timers = ($timer);
	}
	for my $t (@timers) {
	    stop_timer($t);
	    if(defined $timers->{$t}{elapsed}) {
	      printdebug_group("timer","report",$t,sprintf("%.3f s",$timers->{$t}{elapsed}));
	    } else {
	      # no such timer
	    }
	  }
      }
  }

# -------------------------------------------------------------------
{
    my $lookup;
    my $prevkey;
    sub debug_group_add {
	my $group = shift;
	$lookup->{$group} = 1;
    }
    sub debug_group_delete {
	my $group = shift;
	delete $lookup->{$group};
    }
    sub debug_or_group {
	my $group  = shift;
	confess "No debug group defined." if ! defined $group;
	# initialize debug group lookup table
	my $key = $Circos::Configuration::CONF{"debug_group"};
	if($key) {
	    if(! $prevkey || $key ne $prevkey) {
		$lookup = {};		
		for my $g (split($COMMA,$key)) {
		    $lookup->{$g} = 1;
		}
	    }
	    $prevkey = $key;
	}
	return if ! $lookup;
	# groups for which reporting is asked for
	my @groups = ref $group eq "ARRAY" ? @$group : split(/,/,$group);
	# groups for which reporting will be done
	return grep($lookup->{ $_ }, @groups);
    }
}

# -------------------------------------------------------------------
sub printdumper {
    $Data::Dumper::Sortkeys  = 1;
    $Data::Dumper::Indent    = 2;
    $Data::Dumper::Quotekeys = 0;
    $Data::Dumper::Terse     = 0;
    print Dumper(@_) unless Circos::Configuration::fetch_conf("silent") || $Circos::OPT{silent};
}

# -------------------------------------------------------------------
sub printwarning {
    if(Circos::Configuration::fetch_conf("warnings") ||
       Circos::Configuration::fetch_conf("paranoid")) {
	Circos::Error::error("warning","general",join(" ",@_));
    }
    if(Circos::Configuration::fetch_conf("paranoid")) {
	Circos::Error::fatal_error("warning","paranoid");
    }
}

# -------------------------------------------------------------------
sub printinfo {
    if(! @_) {
	printout();
    } else {
	printout( list_as_string(@_) );
    }
}

# -------------------------------------------------------------------
sub printout {
    print "@_\n" unless Circos::Configuration::fetch_conf("silent") || $Circos::OPT{silent};
}

1;
