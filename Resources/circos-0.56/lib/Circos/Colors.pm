package Circos::Colors;

=pod

=head1 NAME

Circos::Colors - Color handling for Circos

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
allocate_colors
allocate_color
rgb_color
rgb_color_opacity
rgb_color_transparency
rgb_to_color
fetch_color
);

use Carp qw( carp confess croak );
use Digest::MD5 qw(md5_hex);
use FindBin;
use File::Basename;
use File::Spec::Functions;
use File::Temp qw(tempdir);
use List::MoreUtils qw( uniq );
use Memoize;
use Math::Round;
use Params::Validate qw(:all);
use Regexp::Common;
use Storable;
use Sys::Hostname;

#use Time::HiRes qw(gettimeofday tv_interval);
#use List::Util qw( max min );

use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";

use Circos::Configuration;
use Circos::Constants;
use Circos::Debug;
use Circos::Error;
use Circos::Image;
use Circos::Utils;

memoize("validate_rgb_list");

# -------------------------------------------------------------------
sub allocate_colors {

    # return undef if ! $CONF{image}{pngmake};
    
    my $image            = shift;
    my $allocated_colors = 0;
    my $colors           = {};
    
    # scan the <colors> block and first allocate all colors
    # specified as r,g,b or r,g,b,a.
    #
    # resolution of name lookups or lists is avoided at this point

    start_timer("colordefinitions");
    for my $color_name ( sort keys %{ $CONF{colors} } ) {
	if(ref $CONF{colors}{$color_name} eq "ARRAY") {
	    my @unique_definitions = uniq @{$CONF{colors}{$color_name}};
	    if(@unique_definitions == 1) {
		printwarning("The color [$color_name] has multiple identical definitions: ".join(" ",@unique_definitions));
		$CONF{colors}{$color_name} = $unique_definitions[0];
	    } else {
		fatal_error("color","multiple_defn",$color_name,join($NEW_LINE, map { " $_" } @unique_definitions));
	    }
	} elsif( my $cref = ref $CONF{colors}{$color_name}) {
	    fatal_error("color","malformed_structure",$color_name,$cref);
	}
	my $color_definition = $CONF{colors}{$color_name};
	if(validate_rgb($color_definition)) {
	    printdebug_group("color","parsing_color RGB",$color_definition);
	    allocate_color($color_name,$color_definition,$colors,$image);
	} elsif (validate_hsv($color_definition)) {
	    my @hsv      = validate_hsv($color_definition);
	    my @rgb255   = hsv_to_rgb(@hsv);
	    my @rgb      = rgb_to_rgb255(@rgb255);
	    my $rgb_text = join(",",@rgb);
	    printdebug_group("color","parsing_color HSV",$color_definition,"RGB",$rgb_text);
	    allocate_color($color_name,$rgb_text,$colors,$image);
	}
    }
    stop_timer("colordefinitions");
    
# now resolve name lookups
    start_timer("colorlookups");
    for my $color_name ( sort keys %{ $CONF{colors} } ) {
	my $color_definition = $CONF{colors}{$color_name};
	# if this color has already been allocated, skip it
	next if exists $colors->{$color_name};
	my %lookup_seen;
	while( exists $CONF{colors}{$color_definition} ) {
	    printdebug_group("color","colorlookup",$color_definition);
	    if($lookup_seen{$color_definition}++) {
		fatal_error("color","circular_defn",$color_definition,$CONF{color}{$color_definition});
	    }
	    $colors->{$color_name} = $colors->{$color_definition};
	    printdebug_group("color","colorlookupassign",$color_name,$color_definition,$CONF{colors}{$color_definition});
	    $color_definition = $CONF{colors}{$color_definition};
	}
    }
    stop_timer("colorlookups");

    # automatic transparent colors
    start_timer("colortransparency");
    create_transparent_colors($colors,$image);
    stop_timer("colortransparency");

    # now resolve lists - employ caching since this can be slow (2-5 seconds);
    start_timer("colorlists");
    my $hostname = hostname;
    my $user     = $ENV{USERNAME} ? $ENV{USERNAME} . $PERIOD : $EMPTY_STRING;
    my $cache_file;
    my $cache_file_root = sprintf("%s.%s.%sdat",
				  Circos::Configuration::fetch_configuration("color_cache_file") || "circos.colorlist",
				  $hostname,
				  $user);				  
    if(my $cache_file_dir = Circos::Configuration::fetch_configuration("color_cache_dir")) {
	$cache_file = catfile($cache_file_dir,$cache_file_root);
    } else {
	# use File::Temp to temporarily create a directory and use this
	# to figure out the system's temporary directory root (e.g. /tmp)
	my $cache_dir      = tempdir();
	rmdir($cache_dir);
	if(! $cache_dir) {
	    fatal_error("io","temp_dir_not_created","color_cache_dir");
	}
	my $cache_dir_root = dirname($cache_dir);
	printdebug_group("cache","temporary file dir",$cache_dir_root);
	$cache_file = catfile($cache_dir_root,$cache_file_root);
    }
    my $allocated_color_list = [keys %$colors];
    my $list_cache;
    my $cache_ok;
    my $is_cache_static = Circos::Configuration::fetch_configuration("color_cache_static");
    my $rebuild_cache   = Circos::Configuration::fetch_configuration("rebuild_color_cache");
    if($rebuild_cache) {
	printdebug_group("cache","colorlist cache rebuild forced");
    } elsif(-e $cache_file) {
	start_timer("colorcache");
	printdebug_group("cache","colorlist cache",$cache_file,"found");
	if ($is_cache_static || -M $cache_file < -M $CONF{configfile}) {
	    printdebug_group("cache","colorlist cache",$cache_file,"useable - static or more recent than configfile");
	    # cache file younger than config file, read cache
	    eval {
		$list_cache = retrieve($cache_file);
	    };
	    if($@) {
		printwarning("Problem reading color cache file $cache_file");
		$cache_ok = 0;
	    } else {
		printdebug_group("cache","colorlist cache",$cache_file,"read in");
		my $target_hash = Digest::MD5::md5_hex(join("", sort keys %{$CONF{colors}}));
		if($list_cache->{colorhash} eq $target_hash) {
		    printdebug_group("cache","color list hash",$target_hash,"matches that of cache file - using cache file");
		    $cache_ok = 1;
		} elsif ($is_cache_static) {
		    printdebug_group("cache","color list hash",$target_hash,"doesn't match that of cache file - using cache anyway because it is static");
		    $cache_ok = 1;
		} else {
		    printdebug_group("cache","color list hash",$target_hash,"does not match - colors changed? - recomputing file");
		}
	    }
	} else {
	    printdebug_group("cache","colorlist cache",$cache_file,"older than configfile - recreating cache");
	}
	stop_timer("colorcache");
    } else {
	printdebug_group("cache","colorlist cache",$cache_file,"not found");
    }
    if(! $cache_ok) {
	# create cache
	$list_cache->{colorhash} = Digest::MD5::md5_hex(join("", sort keys %{$CONF{colors}}));
	printdebug_group("cache","creating colorlist cache, hash",$list_cache->{colorhash});
	for my $color_name ( sort keys %{ $CONF{colors} } ) {
	    # skip if this color has already been allocated
	    next if exists $colors->{$color_name};
	    my @color_definitions = str_to_list($CONF{colors}{$color_name});
	    my @match_set;
	    for my $color_definition (@color_definitions) {
		# do a very quick match to narrow down the colors with fast grep()
		my $rx = $color_definition;
		if($rx =~ /rev\((.+)\)/) {
		    $rx  = $1;
		}
		my @early_matches = grep($_ =~ /$rx/, @$allocated_color_list);
		my @matches;
		# now do a full match, including sorting results
		if(@early_matches) {
		    @matches = sample_list($color_definition,\@early_matches); #$allocated_color_list);
		}
		if(! @matches) {
		    fatal_error("color","bad_pointer",$color_name,$color_definition);
		}
		push @match_set, @matches;
	    }
	    $list_cache->{list2color}{$color_name} = \@match_set;
	    printdebug_group("color","colorlist",$color_name,@match_set);
	}
	# store cache
	eval { 
	    printdebug_group("cache","writing to colorlist cache file [$cache_file]");
	    store($list_cache,$cache_file);
	};
	if($@) {
	    printwarning("Could not write to color list cache file $cache_file - store() gave error");
	    printinfo($@);
	} else {
	    if(-e $cache_file) {
		printdebug_group("cache","wrote to colorlist cache file [$cache_file]");
	    } else {
		printwarning("Could not find the cache file we supposedly just created $cache_file");
	    }
	}
    }
    for my $color (keys %{$list_cache->{list2color}}) {
	$colors->{$color} = $list_cache->{list2color}{$color};
	push @$allocated_color_list, $color;
    }
    stop_timer("colorlists");
    return $colors;
}
 

# -------------------------------------------------------------------
sub rgb_color_opacity {
  # Returns the opacity of a color, based on its name. Colors with a
  # trailing _aNNN have a transparency level in the range
  # 0..auto_alpha_steps. 
  my $color = shift;
  return 1 if ! defined $color;
  if ( $color =~ /(.+)_a(\d+)/ ) {
    unless ( $CONF{image}{auto_alpha_colors}
	     && $CONF{image}{auto_alpha_steps}
	   ) {
      die "you are trying to process a transparent color ($color) ",
	"but do not have auto_alpha_colors or auto_alpha_steps defined";
    }
    my $color_root = $1;
    my $opacity    = 1 - $2 / (1+$CONF{image}{auto_alpha_steps});
  } else {
    return 1;
  }
}


# -------------------------------------------------------------------
sub allocate_color {
    my ($name,$definition,$colors,$image) = @_;
    my @rgb = validate_rgb($definition);
    my $idx;
    printdebug_group("color","allocate_color 0",@rgb);
    if ( @rgb == 3 ) {
	if($name =~ /.+_a\d+$/) {
	    fatal_error("color","reserved_name_a",$name,$definition);
	}
	eval {
	    my $color_index = $image->colorExact(@rgb);
	    if ( $color_index == -1 ) {
		$colors->{$name} = $image->colorAllocate(@rgb);
	    } else {
		$colors->{$name} = $color_index;
	    }
	};
	printdebug_group("color","allocate_color 1",@rgb,$image->colorExact(@rgb));
	if ($@) {
	    fatal_error("color","cannot_allocate",$name,$definition,$@);
	}
    } elsif ( @rgb == 4 ) {
	if($rgb[3] < 0 || $rgb[3] > 127) {
	    fatal_error("color","bad_alpha",$rgb[3]);
	}
	$rgb[3] *= 127 if $rgb[3] < 1;
	eval {
	    printdebug_group("color","allocate_color 2",@rgb);
	    $colors->{$name} = $image->colorAllocateAlpha(@rgb);
	};
	if ($@) {
	    fatal_error("color","cannot_allocate",$name,$definition,$@);
	}
    }
    printdebug_group("color","allocate_color","idx",$colors->{$name},$name,@rgb,"now have",int(keys %$colors),"colors");
}

# -------------------------------------------------------------------
sub create_transparent_colors {
    # Automatically allocate colors with alpha values, if asked for.
    # The number of steps is determined by auto_alpha_steps in the
    # <image> block
    # Colors with alpha values have names COLOR_aN for N=1..num_steps
    # The alpha value (out of max 127) for step i is 127*i/(num_steps+1)
    #
    # For example, if the number of steps is 5, then for the color
    # chr19=153,0,204, the follow additional 5 colors will be
    # allocated (see full list in lines with 'auto_alpha_color' with -debug).
    #
    # Now add automatic transparenc levels to all the defined colors
    # using _aN suffix
    my ($colors,$image) = @_;
    return unless $CONF{image}{auto_alpha_colors};
    my @c = keys %$colors;
    for my $color_name (@c) {
	# if this color is already transparent, skip it
	next if $color_name =~ /.*_a\d+$/;
	my @rgb = $image->rgb( $colors->{$color_name} );
	# provide _a0 synonym
	$colors->{ sprintf("%s_a0",$color_name) } = $colors->{ $color_name };
	for my $i ( 1 .. $CONF{image}{auto_alpha_steps} ) {
	    my $alpha = round( 127 * $i / ( $CONF{image}{auto_alpha_steps} + 1 ) );
	    my $color_name_alpha = $color_name . "_a$i";
	    printdebug_group("color","allocate","auto_alpha_color",$color_name_alpha,@rgb,$alpha);
	    allocate_color($color_name_alpha,[@rgb,$alpha],$colors,$image);
	}
    }
}

sub validate_rgb_list {
  my @rgb = @_;
  my $n = @rgb;
  if($n == 3 || $n == 4) {
    if(grep($_ =~ /$RE{num}{int}/ 
	    &&
	    $_ >= 0 
	    &&
	    $_ <= 255, @rgb) == $n) {
      return @rgb;
    }
  }
  return undef;
}

sub validate_hsv {
    my ($definition,$strict) = @_;
    $strict = 1 if ! defined $strict;
    if($definition =~ /^hsv/i) {
	if($definition =~ /hsv\s*\(\s*(.+)\s*\)/i) {
	    my $hsv_list = $1;
	    my @hsv = str_to_list($hsv_list);
	    return @hsv if @hsv;
	}
	fatal_error("color","malformed_hsv",$definition) if $strict;
	return;
    }
    return;
}

sub validate_rgb {
    my ($definition,$strict) = @_;
    if($definition =~ /^\s*\d+\s*,\s*\d+\s*,\s*\d+\s*(,\s*\d+)?$/ || ref $definition eq "ARRAY") {
	my @rgb = ref $definition eq "ARRAY" ? @$definition : str_to_list($definition);
	return @rgb if @rgb;
    }
    fatal_error("color","malformed_rgb",$definition) if $strict;
    return;
}

# -------------------------------------------------------------------
sub rgb_color_transparency {
  my $color = shift;
  $color = lc $color;
  return 1 - rgb_color_opacity($color);
}

# -------------------------------------------------------------------
sub rgb_color {
  my $color = shift;
  $color = lc $color;
  #confess if ! defined $color;
  return undef if ! defined $color;
  if ( $color =~ /(.+)_a(\d+)/ ) {
      my $color_root = $1;
      return rgb_color($color_root);
  } else {
      return undef unless defined $color && exists $COLORS->{$color};
      #my $colordef = $CONF{colors}{$color};
      my @rgb = $IM->rgb( $COLORS->{$color} );
      #printinfo($color,@rgb);
      return @rgb;
      my $colordef  = $COLORS->{$color};
      if($COLORS->{$colordef}) {
	  return rgb_color($colordef);
      }
      @rgb = split( $COMMA, $colordef );
      return @rgb;
  }
}

sub rgb_to_color {
    my @rgb = @_;
    for my $color (keys %$COLORS) {
	next if $color =~ /_a\d+$/;
	my @crgb = $IM->rgb( fetch_color($color) );
	if($rgb[0] == $crgb[0] &&
	   $rgb[1] == $crgb[1] &&
	   $rgb[2] == $crgb[2]) {
	    return $color;
	}
    }
    fatal_error("color","bad_rgb_lookup",@rgb);
}

sub hsv_to_rgb {
    my ($h, $s, $v) = @_;
    my @rgb;

    $h = $h % 360 if $h < 0 || $h > 360;
    # hue segment 
    $h /= 60;

    my $i = POSIX::floor( $h );
    my $f = $h - $i; 
    my $p = $v * ( 1 - $s );
    my $q = $v * ( 1 - $s * $f );
    my $t = $v * ( 1 - $s * ( 1 - $f ) );

    if($i == 0) {
	@rgb = ($v,$t,$p);
    } elsif ($i == 1) {
	@rgb = ($q,$v,$p);
    } elsif ($i == 2) {
	@rgb = ($p,$v,$t);
    } elsif ($i == 3) {
	@rgb = ($p,$q,$v);
    } elsif ($i == 4) {
	@rgb = ($t,$p,$v);
    } else {
	@rgb = ($v,$p,$q);
    }
    return @rgb;
}

sub rgb_to_rgb255 {
    my @rgb = @_;
    # make sure values are in [0,1]
    @rgb = map { put_between($_,0,1) } @rgb;
    my @rgb255 = map { round( 255 * $_) } @rgb;
    return @rgb255;
}

sub fetch_color {
    my ($color_name,$color_table) = shift;
    $color_table ||= $COLORS;
    if(exists $COLORS->{$color_name}) {
	return $COLORS->{$color_name};
    } elsif ($COLORS->{lc $color_name}) {
	printwarning("Circos colors should be lowercase. You have asked for color [$color_name] and it was interpreted as [".lc $color_name."]");
	return $COLORS->{lc $color_name};
    } else {
	fatal_error("color","undefined",$color_name);
    }
}

1;
