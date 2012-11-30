#!/usr/bin/perl
#############################################
#Author: Jiang Li
#email: riverlee2008@gmail.com
#Creat Time: 
#Vanderbilt Center for Quantitative Sciences
#############################################
#use strict;
#use warnings;
#use GD::Graph::histogram;
#use GD::Graph::boxplot;


use Getopt::Long;

my $savebam=1;

GetOptions("b!"=>\$savebam);

print "savebam=$savebam\n\n";

_samtools_flagstat("Examples/brca_tumor.bam");
sub _samtools_flagstat{
	my ($inbam) = @_;
	my ($total,$mapped)=(0,0);
	open(IN,"samtools flagstat $inbam |") or die $!;
	my @array=<IN>;close IN;
	print join "\n",@array;
	print "\n\n=======\n\n";
	@array=map {/(^\d+)/;$1} @array;
	print join "\n",@array;
	$total=$array[0];
	$mapped=$array[2];
	return ($total,$mapped);
}


my $a=10000000000000000000000000000;
my $b=1234;
print $a,"/",$b,"=",$a/$b,"\n\n";

=head a
#Test function seek
open(IN,"samtools view Examples/TCGA-BH-A0BM-01A-11W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome_MT.bam|") or die $!;
my $five=5;
my $count=1;
while(<IN>){
	my ($r,@t)=split "\t";
	print $count++,"$r\n";
	last if ($count>$five);
}

print "\n\nRedo from begin\n";

$count=1;
seek(IN,0,0);
while(<IN>){
	my ($r,@t)=split "\t";
	print $count++,"$r\n";
	last if ($count>$five);
}


print "\nTest warn\n";
warn "yes, it's a warning\n";

my         $data = [1,5,7,8,9,10,11,3,3,5,5,5,7,2,2];
 my $graph = new GD::Graph::histogram(400,600);
$graph->set( 
                x_label         => 'X Label',
                y_label         => 'Count',
                title           => 'A Simple Count Histogram Chart',
                x_labels_vertical => 1,
                bar_spacing     => 0,
                shadow_depth    => 1,
               # shadowclr       => qw(green pink blue cyan)  ,
               dclrs => [ qw(lgreen pink blue cyan) ],
                transparent     => 0,
        ) 
        or warn $graph->error;
  my $gd = $graph->plot($data) or die $graph->error;
          open(IMG, '>histogram.png') or die $!;
        binmode IMG;
        print IMG $gd->png;
          
    $one = [210..275];
        $two = [180, 190, 200, 220, 235, 245];
        $three = [40, 140..150, 160..180, 250];
        $four = [100..125, 136..140];
        $five = [10..50, 100, 180];

        @data = ( 
                ["1st", "2nd", "3rd", "4th", "5th"],
                [$one, $two, $three, $four, $five ],
                [ [-25, 1..15], [-45, 25..45, 100], [70, 42..125], [undef], [180..250] ],
                # as many sets of data sets as you like         
                );
                
          $my_graph = new GD::Graph::boxplot( );
   $my_graph->set( 
                x_label           => 'X Label',
                y_label           => 'Y label',
                title             => 'Some simple graph',
                upper_percent     => 70,
                lower_percent     => 35,
                step_const        => 1.8
                );
   $gd = $my_graph->plot( \@data );

    open(IMG, '>box.png') or die $!;
    binmode IMG;
    print IMG $gd->png;                           
    
    
use GD::Graph::bars3d;
use GD::Text::Wrap;
        my @data = (
           [qw (Sun Mon Tue Wed Thu Fri Sat) ],
           [ 1100, 1500,  1800,  1950,  2324,  2162,  140],
           [ 103, 2000,  1173,  908,  232,  1200,  100],
        );

my $graph = new GD::Graph::bars3d( 495, 211 );
        $graph->set(
                x_label           => 'Day of the week',
                dclrs             => [ qw(red dgreen) ],
                bar_spacing       => 10,
                cumulate          => 1,
        );
        my ($gd) =  $graph->plot( \@data )->png;
my $text = <<EOSTR;
Lorem ipsum dolor sit amet, consectetuer adipiscing elit,
sed diam nonummy nibh euismod tincidunt ut laoreet dolore
magna aliquam erat volutpat.
EOSTR


  my $wrapbox = GD::Text::Wrap->new( $gd,
      line_space  => 4,
      color       => 'black',
      text        => $text,
  );

binmode STDOUT;
print $gd;

=cut
