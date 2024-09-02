#! /usr/bin/perl

# This script was copied from /opt/phenix/core/bin/preliminary.pl 
# and modified to work on a local machine

# This depends on the PDF::API2 from CPAN.

use PDF::API2;
use Getopt::Long;

sub PrintHelp
  {
    $help = <<END;

Use this command to put a PHENIX preliminary "watermark" on a PDF file 
(usually a plot). Call it like this:

 preliminary.pl --input=in.pdf --output=out.pdf --x=100 --y=100 --scale=0.5

where x and y and scale will need to be adjusted to put the watermark
in a suitable place on your plot.

END
    print $help;
    exit;
  }

$input = '';
$output = '';
$x = 0;
$y = 0;
$scale = 1;
GetOptions('input=s' => \$input,
	   'output=s' => \$output,
	   'x=f' => \$x,
	   'y=f' => \$y,
	   'scale=f' => \$scale,
	   'help' => \&PrintHelp);

die 'You must provide an output file' if $output eq '';
die 'The input and output files must be different' if $input eq $output;

$pdf = PDF::API2->open("$input");
$pdf->default('nounrotate',1);
$watermark = PDF::API2->open("../assets/phenix_preliminary_logo.pdf");
$page = $pdf->openpage(1);
$gfx = $page->gfx();
$xo = $pdf->importPageIntoForm($watermark, 1);
$gfx->formimage($xo, $x, $y, $scale);

$pdf->saveas("$output");
