package Circos::Error;

=pod

=head1 NAME

Circos::Error - error handling for Circos

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
error
fatal_error
);

use Carp qw( carp confess croak );
use Params::Validate;
use Text::Format;

use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";

use Circos::Constants;
use Circos::Debug;
use Circos::Utils;

our %GROUPERROR = (configuration => "configuration file error");
our %ERROR;

$ERROR{color} =
{
    bad_name_in_list => "The color list [%s] included a color definition [%s] that does not match any previously defined color.",
    clear_redefined  => "You defined a color named 'clear' but this name is reserved as synonym for 'transparent'. Plese remove this definition and call your color something else.",
    undefined => "You've asked for color named [%s] but this color has not been defined.\nPlease verify that you've included all the color files you wanted in the <color> block.\nIf you've asked for a transparent color (e.g. blue_a3), make sure that in the <image> block you have auto_alpha_colors=yes and an appropriate value for auto_alpha_steps.",
    multiple_defn => "The color [%s] has multiple distinct definitions - please use only one of these.\n%s",
    malformed_structure => "The color [%s] is not defined correctly. Saw a data structure of type [%s] instead of a simple color assignment.",
    circular_defn => "You have a circular color definition in your <color> block. Color [%s] and [%s] refer to each other.\nYou can define one color in terms of another, such as\n red=255,0,0\n favourite=red\nbut you must avoid loops, such as\n red=favourite\n favourite=red",
    reserved_name_a => "You are trying to allocate color [%s] with definition [%s], but names ending in _aN are reserved for colors with transparency.",
    cannot_allocate => "Could not allocate color [%s] with definition [%s] (error: %s).",
    bad_alpha => "Alpha value of [%s] cannot be used. Please use a range 0-127.",
    malformed_hsv => "HSV definition [%s] is not in the correct format. You must use\n h,s,v = 60,1,0.5\nor\n h,s,v,a = 0,1,0.5,100\nwhere a is the alpha channel (0-127), h is the hue (0-360), s is saturation (0-1) and v is the value (0-1).",
    malformed_rgb => "RGB definition [%s] is not in the correct format. You must use\n r,g,b = 60,120,180\nor\n r,g,b,a = 60,120,180,100\nwhere a is the alpha channel (0-127), and r,g,b is in the range 0-255.",
    bad_rgb_lookup => "Could not find a color with RGB value %d,%d,%d.",
};

$ERROR{ideogram} =
{
    use_undefined => "Entry in 'chromosomes' parameter [%s] mentions chromosome [%s] which is not defined the karyotype file. Make sure that list_record_delim and list_field_delim are defined (see etc/housekeeping.conf) in order for the 'chromosomes' parameter to be parsed correctly.",
    multiple_start_anchors => "Only one chromosome order group can have a start '^' anchor",
    multiple_end_anchors => "Only one chromosome order group can have an end '\$' anchor",
    multiple_tag => "Incorrectly formatted chromosomes_order field (or content of chromosomes_order_file). Tag [%s] appears multiple times, but must be unique.",
    orphan_tag => "Incorrectly formatted chromosomes_order field (or content of chromosomes_order_file). Tag [%s] is not associated with any chromosome.",
    reserved_tag => "You have an ideogram with the tag [%s] which is not allowed, as this is a reserved keyword",
    reserved_chr => "You have an ideogram with the name [%s] which is not allowed, as this is a reserved keyword",
    start_and_end_anchors => "You have a chromosome order group with both start '^' and end '\$' anchors.\n %s\nThis is not supported.\nIf you want to limit which ideograms are drawn, use '-' in front of their names in the chromosomes field. For example,\n chromosomes = -hs1,-hs2",
    cannot_place => "Chromosomes_order string cannot be processed because group\n %s\ncannot be placed in the figure. This may be due to more tags in the chromosomes_order field than ideograms.",
    unparsable_def => "Chromosome definition\n %s\ncould not be parsed. It must be in the format CHR_NAME or CHR_NAME:RUN_LIST. For example\n hs1\n hs1:10-20\n hs1:10-20,30-50\n hs1:(-20,30-)",
    regex_tag => "You have used a regular expression in the 'chromosomes' parameter in the string\n %s\ntogether with a tag [%s]. This combination is not supported.",
    no_ideograms_to_draw => "No ideograms to draw. Either define some in 'chromosomes' parameter or set\n chromosomes_display_default = yes",
    no_such_idx => "Tried to fetch an ideogram with index [%d], but no such ideogram exists.",
    bad_scaled_position => "Could not correctly apply scaling to find pixel position for ideogram [%s] position [%s].",

};

$ERROR{pattern} = 
{
    no_file_def => "You asked for pattern [%s] but it is not associated with a file definition in the <pattern> block.",
    no_file     => "You asked for pattern [%s] but its associated image file [%s] does not exist.",
    cannot_create => "There was a problem creating the pattern [%s] from image file [%s].",
    
};

$ERROR{configuration} = 
  {
   multi_word_key => "Your parameter [%s] contains a white space. This is not allowed. You either forgot a '=' in assignment (e.g. 'red 255,0,0' vs 'red = 255,0,0') or used a multi-word parameter name\n (e.g. 'my red = 255,0,0' vs 'my_red = 255,0,0'",
   missing => "file(error/configuration.missing.txt)",
   bad_parameter_type => "You attempted to reference a configuration parameter group with name [%s], but it is not defined",
   defined_twice => "Parameter [%s] of type [%s] is defined twice. Some parameters can have multiple values, but not this one.",
   multiple_defn_in_list => "Configuration value [%s] defines parameter [%s] more than once. This is not allowed.",
   multivalue => "Configuration parameter [%s] has been defined more than once in the block shown above, and has been interpreted as a list. This is not allowed. Did you forget to comment out an old value of the parameter?",
   unsupported_parameter => "Parameter [%s] of type [%s] is not supported.",
   bad_pointer => "Problem with variable lookup in configuration file. You referenced variable [%s] (seen as %s) in another parameter, but this variable is not defined.",
   no_karyotype => "You did not define a karyotype file. Are you sure your configuration file is well formed?",
  };

$ERROR{font} = 
{
    no_def  => "Non-existent font definition for font [%s] requested for [%s].",
    no_file => "Could not find file for font [%s] which has definition [%s] requested for [%s]",
    no_name => "Non-existent font definition for font [%s] requested for [%s].",
    no_ttf  => "There was a problem with True Type font support. Circos could not render text from the font file\n  %s\nPlease check that gd (system graphics library) and GD (Perl's interface to gd) are compiled with True Type support.\nOn UNIX systems, try\n  gdlib-config --all\nand look for GD_FREETYPE in the 'features' line and -lfreetype in the 'libs' line. If these are there, it's likely that your Perl GD module needs recompiling.\nFor help in installing libgd and/or GD, see\n  http://www.perlmonks.org/?node_id=621579",
};

$ERROR{argument} =
{
    list_size => "Function [%s] in package [%s] expected a list of size [%d] but only saw [%d] elements.",
};

$ERROR{geometry} =
{
    angle_out_of_bounds => "The angle [%f] is out of bounds. Expected value in [-90,270].",
};

$ERROR{graphics} =
{
    brush_zero_size => "Cannot create a brush with zero width",
};

$ERROR{rules} =
{
    no_such_field => "You set up a rule [%s] that uses the parsable field [%s] but the data point you are testing does not have the field [%s].",
    no_field_value => "You set up a rule [%s] that uses the parsable field [%s], but this field has no associated value.",
    wrong_num_elements => "You set up a rule [%s] that uses the parsable field [%s] but the data point you are testing does not have [%d] elements.",
    parse_error => "There was a problem evaluating the string [%s] as code (error: %s)",
};

$ERROR{io} =
{
    temp_dir_not_created => "Attempted to automatically determine a directory for temporary files, but failed. Please set the %s parameter to define it.",
    cannot_write => "Cannot open file [%s] for writing %s. (error %s)",
    cannot_read => "Cannot read file [%s] for reading %s. (error %s)",
    cannot_find => "Cannot guess the location of file [%s]. Tried to look in the following directories\n%s",
};

$ERROR{warning} =
{
    general => "%s",
    paranoid => "Circos produced a warning and quit because you are running it in paranoid mode. The warning does not nececessarily mean that something is wrong - Circos was likely trying to guess a parameter or massage data to fit the figure.\nTo change this setting, see etc/housekeeping.conf.",
};

$ERROR{track} =
{
    min_larger_than_max => "Plot min value [%f] is larger than max [%f]",
    start_larger_than_end => "Input data line in file [%s] for track type [%s] has start position [%s] greater than end position [%s].",

};

$ERROR{links} = 
{
    duplicate_names => "Multiple link data sets with name [%s] are defined. This is not supported.",
    single_entry => "Link data file for [%s] has a single positional entry for link [%s]. A link must have two entries - a start and an end.",
    too_thick =>"You are attempting to draw a bezier curve of thickness greater than 100 [%d]. This would take a very long time and you don't want to do this.",
    too_thin =>"You are attempting to draw a bezier curve of thickness less than 1 [%d]. This would produce nothing. Is this what you want? If so, hide the link. If not, set the thickness to be at least 1.",

};

$ERROR{map} = 
{
    url_param_not_set => "You have tried to use the URL [%s] for an image map, but the parameter in the url [%s] has no value defined for this data point or data set.\nTo make this error go away, either (a) define the parameter, (b) set\n image_map_missing_parameter = blank\nto remove the undefined parameter from the image element, or (c) set\n image_map_missing_parameter = removeurl\nto remove the URL from the image element.",
};

$ERROR{data} =
{
    malformed_span => "There was a problem initializing a span. Saw start [%s] > end [%s].",
    repeated_chr_in_karyotype => "Chromosome [%s] defined more than once in karyotype file.",
    malformed_karyotype_coordinates => "Start [%s] and/or end [%s] coordinate in karyotype file are not digits.",
    unknown_karyotype_line => "You have a line named [%s] in the karyotype file but currently only 'chr' or 'band' lines are supported.",
    inconsistent_karyotype_coordinates => "Start [%s] must be smaller than end [%s] coordinate in karyotype file.",
    band_on_missing_chr => "Bands for chromosome [%s] are defined but the chromosome itself has no definition.\nIs there a 'chr' line for this chromosome in the karyotype file?",
    band_sticks_out => "Band [%s] on chromosome [%s] has coordinates that extend outside chromosome.",
    band_overlaps => "Band [%s] on chromosome [%s] overlaps with another band by more than [%s].",
};

$ERROR{function} =
{
    remap_wrong_num_args => "You must pass five parameters to remap(). The syntax is\n remap(value,min,max,remap_min,remap_max)",
    remap_min_max =>        "You must pass five parameters to remap(). The syntax is\n remap(value,min,max,remap_min,remap_max)\nYou have set min [%s] = max [%s], but this only makes sense if remap_min [%s] = remap_max [%s].",
    pairwise => "You have tried to use the %s(a,b,x,y) function but one of the arguments a=[%s] b=[%s] c=[%s] d=[%s] is undefined.",
    sample_list_bad_arg => "Argument to sample_list must be a list but saw [%s] of type [%s].",
};

$ERROR{unit} =
{
    conversion_fail => "Unable to convert a value [%s] from one unit [%s] to another [%s]. The following from->to combinations were expected: %s",
};

$ERROR{system} =
{
    bad_error_name => "What do you know - a fatal error caused by bad error handling.\nThe error category [%s] and name [%s] is invalid.",
    missing_units_ok => "The parameter 'units_ok' is not defined (usually found in etc/housekeeping.conf). This parameter defines allowable units and is required. Set it to\n units_ok = %s",
    missing_units_nounit => "The parameter 'units_nounit' is not defined (usually found in etc/housekeeping.conf). This parameter defines the explicit suffix for unitless quantities. Set it to\n units_nounit = %s",
    wrong_unit => "The parameter [%s] value [%s] does not have the correct unit. Saw [%s] but expected one of [%s].",
    undef_parameter => "The parameter [%s] was not defined. It needs to be defined and have one of these units [%s].",
    unit_format_fail => "The unit [%s] failed format check. The list of allowable units is set by 'units_ok' (usually in etc/housekeeping.conf).",
    bad_dimension => "Dimension [%s] is not defined in expression [%s]",
};

$ERROR{support} =
{
    googlegroup=>"If you are having trouble debugging this error, use this tutorial to learn how to use the debugging facility\n http://www.circos.ca/tutorials/lessons/configuration/debugging\nIf you're still stumped, get support in the Circos Google Group\n  http://groups.google.com/group/circos-data-visualization",
};

sub error {
  my ($cat,$errorid,@args) = @_;
  my $error_text = $ERROR{$cat}{$errorid};
  if(! defined $error_text) {
      fatal_error("system","bad_error_name",$cat,$errorid);
  }
  if($error_text =~ /file\((.*)\)/) {
      my $file = Circos::Utils::locate_file(file=>$1,name=>"error file",return_undef=>1);
      if($file && open(F,$file)) {
	  $error_text = join("",<F>);
	  close(F);
      } else {
	  $error_text = "...error text from [$1] could not be read...";
      }
  }
  my (@text,$format);
  my $undef_text = Circos::Configuration::fetch_conf("debug_undef_text") || "_undef_";
  @args = map { defined $_ ? $_ : $undef_text } @args;
  if($cat eq "warning") {
      if(Circos::Configuration::fetch_conf("paranoid")) {
	  @text = ("*** CIRCOS WARNING ***",
		   uc $GROUPERROR{$cat},
		   sprintf($error_text,@args));
	  $format = 1;
      } else {
	  printdebug_group("!circoswarning",uc $GROUPERROR{$cat},sprintf($error_text,@args));
	  return;
      }
  } else {
      @text = ("*** CIRCOS ERROR ***",
	       uc $GROUPERROR{$cat},
	       sprintf($error_text,@args),
	       $ERROR{support}{googlegroup},
	       "",
	       "Stack trace:",
	  );
      $format = 1;
  }
  if($format) {
      my @text_fmt;
      for my $t (@text) {
	  for my $line (split(/\n/, $t)) {
	      if($line =~ /^\s/) {
		  push @text_fmt, Text::Format->new({leftMargin=>6,columns=>80,firstIndent=>0})->paragraphs($line);
	      } else {
	      push @text_fmt, Text::Format->new({leftMargin=>2,columns=>80,firstIndent=>0})->paragraphs($line);
	      }
	  }
      }
      print "\n" . join("\n", @text_fmt);
  } else {
      print join(" ",@text)."\n";
  }
}

sub fatal_error {
    error(@_);
    confess;
}

1;
