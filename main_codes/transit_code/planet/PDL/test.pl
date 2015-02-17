use blib;

use strict;
use PDL;
use PDL::Planet;
use constant PI => 4 * atan2(1, 1);
use PDL::Graphics::PGPLOT::Window;

#for (my $i = 0; $i<= 20; $i++){
#    my $phase = ($i/20) * 37.18;
#    my $res = PDL::Planet::RV($phase, 0.98e0, 37.18e0, 0.755 + 5. * PI / 180., 0.361, 0.0e0, 37.18e0);
#    printf("%10.7f %10.7f\n", $res, $phase/37.18);
#}

my $t = 2 * 30.433 * sequence(1000)/999;

my $y = PDL::Planet::RV($t, 0, 968.6e0, (189.1 )* PI / 180, 0.5170, 0.0e0, 5.633e0);

my $win = PDL::Graphics::PGPLOT::Window->new({Device=>"Test.eps/vcps",
                                               AspectRatio =>1,
                                               YTitle => 'v (m sec\u-1\d)',
                                               XTitle => 'Phase'});

$win->env(-0.1, 2.1, -1800, 750);
$win->lines($t / 5.633, $y);

my $y2 = PDL::Planet::RV($t, 0, 968.6e0, (189.1 + 0.5 )* PI / 180, 0.5170, 0.0e0, 5.633e0);
#$win->lines($t / 5.633, $y2, {Color => 'red'});

$win->close();

$win = PDL::Graphics::PGPLOT::Window->new({Device=>"Test2.eps/vcps",
                                               AspectRatio =>1,
                                               YTitle => '\gD v (m sec\u-1\d)',
                                               XTitle => 'Phase'});

$win->env(-0.1, 2.1, -20, 20);
$win->lines($t / 5.633, $y-$y2);


$win->close();

$t = 2 * 460.0 * sequence(1000)/999;
$y = PDL::Planet::RV($t, 0, 51, (271 )* PI / 180, 0.53, 0.0e0, 460.0);
$win = PDL::Graphics::PGPLOT::Window->new({Device=>"Test3.eps/vcps",
                                               AspectRatio =>1,
                                               YTitle => 'v (m sec\u-1\d)',
                                               XTitle => 'Phase'});

$win->env(-0.1, 2.1, -60, 60);
$win->lines($t / 460.0, $y);


$win->close();

$t = sequence(100)/99;
my $e = 0.5;
my $w = 0.0;
my $Omega = 0.0;
my $t0 = 0;
my $Period = 1;
my $inc = PI/2 ;
my $a = 1;

my ($X,$Y,$Z) = PDL::Planet::Orbit($t,$t0,$Period,$a,$e,$w,$inc, $Omega);

$win = PDL::Graphics::PGPLOT::Window->new({Device=>"Test4.eps/vcps",
					   AspectRatio =>1,
					   YTitle => 'Y',
					   XTitle => 'X'});

$win->env(-1.5,1.5,-1.5,1.5);
$win->points($X,$Y);
$win->close();
