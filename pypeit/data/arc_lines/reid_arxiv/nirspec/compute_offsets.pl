use PDL;
use PDL::Graphics::PGPLOT::Window;
use PDL::Image2D;
use PDL::Fit::Gaussian;

sub fitgauss {

  my $xs = shift;
  my $ys = shift;

  my($xcen,$pk,$fwhm2,$back,$err,$fit) = fitgauss1d($xs,$ys);
  
  return($xcen,$fwhm2/2.35,$pk);
}

sub cent {

  my $stamp = shift;
  my $offset = shift;
  my $win;
  if (@_) {
    $win = shift;
  }

  $stamp = $stamp->medover;
  $stamp->where($stamp < 10 ) .=0;
  my $xs = $stamp->xvals;
  $xs += $offset;
  my ($cenx,$sigx,$pkflux) = fitgauss($xs,$stamp);

  if ($win) {
    $win->line($xs,$stamp);
    $win->hold;
    $win->line(pdl([$cenx,$cenx]),pdl([-1e6,1e6]));
    $win->release;
  }

 
  return($cenx,$sigx,$pkflux);
}

my $win = PDL::Graphics::PGPLOT::Window->new(device=>'/xs');

opendir(D,".");
my @sfos = grep { /tdfor\.fits/ } readdir(D);
closedir(D);

my $slyb = 70;
my $slye = 290;

my $sl = "484:579,$slyb:$slye";
my $first = 1;
my $sfoname;
foreach $sfoname (@sfos) {
  my $sfo = rfits "$sfoname";
  print "$sfoname ";
  $win->erase;
  $win->env($slyb,$slye,-10,100,{PlotPosition=>[0.1,.99,.2,.4]});
  my($cenx,$sigx,$pkflux) = cent($sfo->slice($sl),$slyb,$win);
  if ($first) {
    $slx = $cenx;
  } else {
    $cenx -= $slx;
  }
  printf "%.2f %.2f %.2f\n",$cenx, $sigx,$pkflux;
  $first = 0;
#  $_ = <STDIN>;
}
