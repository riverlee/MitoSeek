package Circos::Configuration;

=pod

=head1 NAME

Circos::Configuration - Configuration handling for Circos

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
our @EXPORT_OK = qw(%CONF $DIMS);
our @EXPORT    = qw(fetch_configuration fetch_conf %CONF $DIMS);

use Carp qw( carp confess croak );
use Config::General 2.50;
use File::Basename;
use File::Spec::Functions;
use Data::Dumper;
use Math::VecStat qw(sum min max average);
use Math::Round qw(round);
use FindBin;
use IO::File;
use Params::Validate qw(:all);
use List::MoreUtils qw(uniq);

use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";

use Circos::Constants;
use Circos::Debug;
use Circos::Utils;
use Circos::Error;

our %CONF;
our $DIMS;

# -------------------------------------------------------------------
sub increment_counter {
  my ($counter,$value) = @_;
  $CONF{counter}{$counter} ||=0;
  $CONF{counter}{$counter} += $value;
  printdebug_group("counter","incrementing counter",$counter,$value,"now",$CONF{counter}{$counter});
}
sub set_counter {
  my ($counter,$value) = @_;
  $CONF{counter}{$counter} = $value;
  printdebug_group("counter","set counter",$counter,$value,"now",$CONF{counter}{$counter});
}
sub init_counter {
  my ($counter,$value) = @_;
  $CONF{counter}{$counter} = $value if ! defined $CONF{counter}{$counter};
  printdebug_group("counter","init counter",$counter,$value,"now",$CONF{counter}{$counter});
}

# -------------------------------------------------------------------
# Return the configuration hash leaf for a parameter path.
#
# fetch_configuration("ideogram","spacing")
#
# returns
#
# $CONF{ideogram}{spacing}
#
# If the leaf, or any of its parents, do not exist, undef is returned.
sub fetch_configuration {
    my @config_path = @_;
    my $conf_item   = \%CONF;
    for my $path_element (@config_path) {
	if(! exists $conf_item->{$path_element}) {
	    return undef;
	} else {
	    $conf_item = $conf_item->{$path_element};
	}
    }
    return $conf_item;
}
sub fetch_conf {
    return fetch_configuration(@_);
}

# -------------------------------------------------------------------
# 
#
#
#
sub fetch_parameter_list_item {
    my ($list,$item,$delim) = @_;
    my $parameter_hash = make_parameter_list_hash($list,$delim);
    return $parameter_hash->{$item};
}

# -------------------------------------------------------------------
# Given a string that contains a list, like
#
# hs1:0.5;hs2:0.25;hs3:0.10;...
#
# returns a hash keyed by the first field before the delimiter with
# the second field as value.
#
# { hs1=>0.5, hs2=>0.25, hs3=>0.10, ... }
#
# The delimiter can be set as an optional second field. By default,
# the delimiter is \s*[;,]\s*
#
sub make_parameter_list_hash {
    my ($list_str,$record_delim,$field_delim) = @_;
    $record_delim ||= fetch_configuration("list_record_delim") || qr/\s*[;,]\s*/;
    $field_delim  ||= fetch_configuration("list_field_delim")  || qr/\s*[:=]\s*/;
    my $parameter_hash;
    for my $pair_str (split($record_delim,$list_str)) {
	my ($parameter,$value) = split($field_delim,$pair_str);
	if(exists $parameter_hash->{$parameter}) {
	    fatal_error("configuration","multiple_defn_in_list",$list_str,$parameter);
	} else {
	    $parameter_hash->{$parameter} = $value;
	}
    }
    return $parameter_hash;
}

# -------------------------------------------------------------------
# Given a string that contains a list, like
#
# file1,file2,...
#
# returns an array of these strings.
#
# [ "file1", "file2", ... ]
#
# The delimiter can be set as an optional second field. By default,
# the delimiter is \s*[;,]\s*
#
sub make_parameter_list_array {
    my ($list_str,$record_delim) = @_;
    $record_delim ||= fetch_configuration("list_record_delim") || qr/\s*[;,]\s*/;
    my $parameter_array;
    for my $str (split($record_delim,$list_str)) {
	push @$parameter_array, $str;
    }
    return $parameter_array;
}

# -------------------------------------------------------------------
sub populateconfiguration {

  my %OPT = @_;

  for my $key ( keys %OPT ) {
      if(defined $OPT{$key}) {
	  if($key eq "debug_group") {
	      $CONF{$key} .= ",".$OPT{$key};
	  } else {
	      $CONF{$key} = $OPT{$key};
	  }
      }
  }

  # any configuration fields of the form __XXX__ are parsed and replaced
  # wiht eval(XXX).
  # The configuration can therefore depend on itself.
  #
  # flag = 10
  # note = __2*$CONF{flag}__ # would become 2*10 = 20

  resolve_synonyms( \%CONF, [] );
   
  repopulateconfiguration( \%CONF );
  
  override_values( \%CONF );

  check_multivalues( \%CONF );

  # populate some defaults
  $CONF{'anglestep'}    ||= 1;
  $CONF{'minslicestep'} ||= 5;

}

sub resolve_synonyms {
    my $root = shift;
    my $tree = shift;
    if(ref $root eq "HASH") {
	for my $key (keys %$root) {
	    my $value = $root->{$key};
	    if(ref $value eq "HASH" ) {
		resolve_synonyms($value, [ @$tree, $key ]);
	    } elsif (ref $value eq "ARRAY") {
		map { resolve_synonyms($_, [ @$tree, $key ]) } @$value;
	    } else {
		my ($new_key,$action) = apply_synonym($value,$key,$tree);
		if(defined $new_key) {
		    if($action eq "copy") {
			$root->{$new_key} = $root->{$key};
		    } else {
			$root->{$new_key} = $root->{$key};
			delete $root->{$key};
		    }
		}
	    }
	}
    }
}

sub apply_synonym {
    my ($value,$key,$tree) = @_;

    my @synonyms = (
	{ key_rx => ".*::label_tangential", new_key => "label_parallel", action => "copy" },
	);

    my $key_name  = join(":",@$tree)."::".$key;
    my ($new_key,$action);
    for my $s (@synonyms) {
	printdebug_group("conf","testing synonym",$s->{key_rx},$key_name);
	if($key_name =~ /$s->{key_rx}/) {
	    $new_key = $s->{new_key};
	    $action  = $s->{action};
	    printdebug_group("conf","applying synonym",$action,$key_name,$new_key);
	    return ($new_key,$action);
	}
    }
    return;
}

# -------------------------------------------------------------------
# Parameters with *, **, ***, etc suffixes override those with
# fewer "*";
sub override_values {
    my $root = shift;
    my $parameter_missing_ok = 1;
    if(ref $root eq "HASH") {
	my @keys = keys %$root;
	# do we have any parameters to override?
	my @lengths = uniq ( map { length($_) } ( map { $_ =~ /([*]+)$/ } @keys ) );
	for my $len (sort {$a <=> $b} @lengths) {
	    for my $key (@keys) {
		my $rx = qr/^(.+)[*]{$len}$/;
		#printinfo("rx",$rx,"key",$key);
		if($key =~ $rx) {
		    my $key_name = $1;
		    # do not require that the parameter be present to override it
		    if($parameter_missing_ok || grep($_ eq $key_name, @keys)) {
			#printinfo("overriding",$key_name,$root->{$key_name},$key,$root->{$key});
			$root->{$key_name} = $root->{$key};
		    }
		}
	    }
	}
	for my $key (keys %$root) {
	    my $value = $root->{$key};
	    if(ref $value eq "HASH" ) {
		#printinfo("iter",$key);
		override_values($value);
	    }
	}
    }
}

# -------------------------------------------------------------------
# Parameters with *, **, ***, etc suffixes override those with
# fewer "*";
sub check_multivalues {
    my $root = shift;

    return unless ref $root eq "HASH";

    my @keys = keys %$root;
    for my $key (@keys) {
	my $value = $root->{$key};
	if(ref $value eq "ARRAY") {
	    my @list_ok = qw(rule tick plot radius zoom highlight);
	    if(! grep($key eq $_, @list_ok)) {
		printdumper($root);
		fatal_error("configuration","multivalue",$key);
	    }
	} elsif (ref $value eq "HASH") {
	    # this is a block
	}
    }

    for my $key (keys %$root) {
	my $value = $root->{$key};
	if(ref $value eq "HASH" ) {
		check_multivalues($value);
	} elsif (ref $value eq "ARRAY") {
	    map { check_multivalues($_) } @$value;
	}
    }
}

# -------------------------------------------------------------------
sub repopulateconfiguration {
  my ($root,$parent) = @_;
  if (my $value = $root->{init_counter}) {
      for my $counter_txt (split(",",$value)) {
	  my ($counter,$value) = split(":",$counter_txt);
	  init_counter($counter,$value);
      }
  }
  if (my $value = $root->{pre_increment_counter}) {
      for my $counter_txt (split(",",$value)) {
	  my ($counter,$incr) = split(":",$counter_txt);
	  increment_counter($counter,$incr);
      }
  }
  if (my $value = $root->{pre_set_counter}) {
      for my $counter_txt (split(",",$value)) {
	  my ($counter,$incr) = split(":",$counter_txt);
	  set_counter($counter,$incr);
      }
  }

  for my $key ( keys %$root ) {
    my $value = $root->{$key};
    if ( ref($value) eq 'HASH' ) {
      $CONF{counter}{$key} ||= 0;
      repopulateconfiguration($value,$key);
    } elsif ( ref($value) eq 'ARRAY' ) {
      for my $item (@$value) {
	if (! defined $CONF{counter}{$key}) {
	  printdebug_group("counter","zeroing counter",$key);
	  $CONF{counter}{$key} = 0;
	}
	repopulateconfiguration($item,$key) if ref($item);
	increment_counter($key, 
			  exists $root->{increment_counter} ? $root->{increment_counter} : 
			  exists $item->{increment_counter} ? $item->{increment_counter} : 1) if ref($root) && ref($item);
      }
    } else {
      if ($key =~ /\s+/) {
	  fatal_error("configuration","multi_word_key",$key); 
      }
      my $delim = "__";
      while ( $value =~ /$delim([^_].+?)$delim/g ) {
	my $source = $delim . $1 . $delim;
	my $target = eval $1;
	printdebug_group("conf","repopulate","key",$key,"value",$value,"var",$1,"target",$target);
	$value =~ s/\Q$source\E/$target/g;
	printdebug_group("conf","repopulate",$key,$value,$target);
      }
      if ($value =~ /\s*eval\s*\((.+)\)/ && $parent ne "rule") {
	my $fn = $1;
	$value = eval $fn;
	if ($@) {
	    fatal_error("rules","parse_error",$fn,$@);
	}
	printdebug_group("conf","repopulateeval",$fn,$value);
      }
      $root->{$key} = $value;
    }
  }
  if (my $value = $root->{post_increment_counter}) {
      for my $counter_text (split(",",$value)) {
	  my ($counter,$incr) = split(":",$counter_text);
	  increment_counter($counter,$incr);
      }
  }
  if (my $value = $root->{post_set_counter}) {
      for my $counter_text (split(",",$value)) {
	  my ($counter,$incr) = split(":",$counter_text);
	  set_counter($counter,$incr);
      }
  }
}

# -------------------------------------------------------------------
sub loadconfiguration {
  my $arg = shift;
  printdebug_group("conf","looking for conf file",$arg);
  my @possibilities = (
		       $arg,
		       catfile( $FindBin::RealBin, $arg ),
		       catfile( $FindBin::RealBin, '..', $arg ),
		       catfile( $FindBin::RealBin, 'etc', $arg ),
		       catfile( $FindBin::RealBin, '..', 'etc', $arg ),
		       catfile( '/home', $ENV{'LOGNAME'}, ".${APP_NAME}.conf" ),
		       catfile( $FindBin::RealBin, "${APP_NAME}.conf" ),
		       catfile( $FindBin::RealBin, 'etc', "${APP_NAME}.conf"),
		       catfile( $FindBin::RealBin, '..', 'etc', "${APP_NAME}.conf"),
		       
		      );
  
  my $file;
  for my $f ( @possibilities ) { 
    if ( -e $f && -r _ ) {
      printdebug_group("conf","found conf file",$f);
      $file = $f;
      last;
    }
  }
  
  if ( !$file ) {
    fatal_error("configuration","missing");
  }
  
  my @configpath = (
		    dirname($file),
		    dirname($file)."/etc",
		    "$FindBin::RealBin/etc", 
		    "$FindBin::RealBin/../etc",
		    "$FindBin::RealBin/..",  
		    $FindBin::RealBin,
		   );
    
  my $conf = Config::General->new(
				  -SplitPolicy       => 'equalsign',
				  -ConfigFile        => $file,
				  -AllowMultiOptions => 1,
				  -LowerCaseNames    => 1,
				  -IncludeAgain      => 1,
				  -ConfigPath        => \@configpath,
				  -AutoTrue => 1
				 );
  
  %CONF = $conf->getall;
}

# -------------------------------------------------------------------
sub validateconfiguration {

  for my $parsekey ( keys %CONF ) {
    if ( $parsekey =~ /^(__(.+)__)$/ ) {
      if ( !defined $CONF{$1} ) {
	  fatal_error("configuration","bad_pointer",$1,$2);
      }

      my ( $token, $parsevalue ) = ( $1, $CONF{$1} );
      for my $key ( keys %CONF ) {
	$CONF{$key} =~ s/$token/$parsevalue/g;
      }
    }
  }

  if($CONF{debug} && ! $CONF{debug_group}) {
      $CONF{debug_group} = "summary,io,timer";
  }

  $CONF{chromosomes_units} ||= 1;
  $CONF{svg_font_scale} ||= 1;

  if ( ! $CONF{karyotype} ) {
      fatal_error("configuration","no_karyotype");
  }

  $CONF{image}{image_map_name} ||= $CONF{image_map_name};
  $CONF{image}{image_map_use}  ||= $CONF{image_map_use};
  $CONF{image}{image_map_file} ||= $CONF{image_map_file};
  $CONF{image}{image_map_missing_parameter} ||= $CONF{image_map_missing_parameter};
  $CONF{image}{"24bit"} = 1;
  $CONF{image}{png}  = $CONF{png} if exists $CONF{png};
  $CONF{image}{svg}  = $CONF{svg} if exists $CONF{svg};
  $CONF{image}{file} = $CONF{outputfile} if $CONF{outputfile};
  $CONF{image}{dir}  = $CONF{outputdir}  if $CONF{outputdir};

  if ( $CONF{image}{angle_offset} > 0 ) {
    $CONF{image}{angle_offset} -= 360;
  }

  #
  # Make sure these fields are initialized
  #

  for my $fld ( qw(chromosomes chromosomes_breaks chromosomes_radius) ) {
    $CONF{ $fld } = $EMPTY_STR if !defined $CONF{ $fld };
  }
  
}

1;
