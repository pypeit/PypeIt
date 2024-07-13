use PDL;
use PDL::Graphics::PGPLOT::Window;
use PDL::Image2D;
use PDL::Fit::Gaussian;
use lib qw(/Users/holden/Documents/my_perl);
use Statistics;
use strict;

sub gaussmod {

  my $mean = shift;
  my $sig = shift;
  my $norm = shift;
  my $xs = shift;

  my $modvals = $xs->zeroes;
  $modvals += ($xs-$mean)**2;
  $modvals /= 2*$sig**2;
  $modvals = exp(-$modvals);
  $modvals /= $modvals->sum;
  $modvals *= $norm;

  return($modvals);

}

sub like {

  my $modvals = shift;
  my $ys = shift;

  my $like = ($modvals - $ys)**2;
#  $like /= $modvals;
  $like = $like->sum;

  return($like);

}

sub compare {

  my $guesses = shift;
  my $xs = shift;
  my $ys = shift;

  my $norm = $ys->sum;

  my $modvals = gaussmod($$guesses{mean},$$guesses{sig},$norm,$xs);
  my $chisq = like($modvals,$ys);

  return($chisq);

}

sub makeguesses {

  my $xcenguess = shift;
  my $xsigguess = shift;

  my %simplex;

  $simplex{mean} =zeroes 3;
  $simplex{sig} =zeroes 3;

  $simplex{mean} .= $xcenguess;
  $simplex{sig} .= $xsigguess;

  $simplex{mean} += 5*(random(3) - 0.5);
  $simplex{sig} += 5*(random(3) - 0.5);

  return(\%simplex);

}



sub amotry {

  my $simplex = shift;
  my $chis = shift;
  my $hiindex = shift;
  my $fac = shift;

  my @params = @_;

  my $fac1 = 1.0-$fac;
  $fac1 /= 2;
  my $fac2 = $fac1 - $fac;

  my %point = makepoint($simplex,$hiindex);
  my $param;
  my %psum;
  foreach $param (@params) {
    $psum{$param} = sum($$simplex{$param});
  }
  foreach $param (@params) {
    $point{$param} = $psum{$param}*$fac1 -
      $fac2*$$simplex{$param}->slice("($hiindex)");
  }

  my $chi = compare(\%point,$$simplex{xs},$$simplex{ys});
  if ($chi < $chis->slice("($hiindex)")) {
    $chis->slice("($hiindex)") .= $chi;
    foreach $param (@params) {
      $$simplex{$param}->slice("($hiindex)") .= $point{$param};
    }
  }
  return($chi);
}


sub makepoint {

  my $simplex = shift;
  my $index = shift;

  my %point;
  my $param;
  $point{mean} = $$simplex{mean}->slice("($index)");
  $point{sig} = $$simplex{sig}->slice("($index)");
  return(%point);
}


sub amoeba {

  my $simplex = shift;

  my @params = ("mean","sig");
  my $chis = zeroes (3);
  my $n;
  # I need to setup the chis PDL for the input parameters
  # the compare subroutine needs a point structure, so I make
  # one and copy the data for each corner of the simplex into it

  my $param;
  my %point = makepoint($simplex,0);
  $chis->slice("0") += compare(\%point,$$simplex{xs},$$simplex{ys});
  foreach $n (1..2) {
    %point = makepoint($simplex,$n);
    $chis->slice("$n") += compare(\%point,$$simplex{xs},$$simplex{ys});
  }

  my $niter  =0;
  my $tiny = 1e-10;
  my $ftol = 1e-8;
  my $lowindex;
  while ($niter < 10000) {

    # find best and worst values.
    my $hiindex = maximum_ind($chis);
    my $hi = maximum($chis);
    $lowindex = minimum_ind($chis);
    my $low = minimum($chis);

    my $rtol = 2.0*($hi - $low)/($hi+$low + $tiny);
    # tolerance value;
    if ($rtol < $ftol) {
      # hey! we are done
      return($lowindex,$chis->slice("($lowindex)"));
    }
    # first stab
    my $trychi = amotry($simplex,$chis,$hiindex,-1.0,@params);
    if ($trychi < $low) {
      amotry($simplex,$chis,$hiindex,2.0,@params);
      # keep going a good direction
    } elsif ($trychi > $hi) {
      # really bad guess
      foreach $n (0..2) {
	next if $n == $lowindex;
	foreach $param (@params) {
	  $$simplex{$param}->slice("($n)") .=
	    0.5*($$simplex{$param}->slice("($n)") +
		 $$simplex{$param}->slice("($lowindex)"));
	  $point{$param} = $$simplex{$param}->slice("($n)");
	}
	$chis->slice("$n") .= compare(\%point,$$simplex{xs},$$simplex{ys});
      }
    }
    $lowindex = minimum_ind($chis);
#    print "$niter ",$chis->slice("($lowindex)"),"\n";
    foreach $param (@params) {
#      print $param," ",$$simplex{$param}->slice("($lowindex)"),"\n";
    }
    $niter++;
  }
 warn "amoeba did not converge\n";
  return($lowindex,$chis->slice("($lowindex)"));
}


sub myfitgauss1d {


  my $xs = shift;
  my $ys = shift;

  my $yscopy = $ys->copy;
  $yscopy->where($yscopy < 10) .= 0;
  $ys->where($ys < -15) .= 0;

  my($xcen,$pk,$fwhm2,$back,$err,$fit) = fitgauss1d($xs,$ys);
  my $sig = $fwhm2 / 2.35;

  my $simplex = makeguesses($xcen,$sig);
  $$simplex{xs} = $xs;
  $$simplex{ys} = $ys;
  my($lowindex,$lowchi) = amoeba($simplex);
  my $cen = $$simplex{mean}->slice("($lowindex)");
  my $sig = $$simplex{sig}->slice("($lowindex)");
  return($cen,$sig);

}

sub fitgauss {

  my $xs = shift;
  my $ys = shift;

#  my($xcen,$pk,$fwhm,$back,$err,$fit) = fitgauss1d($xs,$ys);
  my($xcen,$fwhm) = myfitgauss1d($xs,$ys);
  
  return($xcen,$fwhm/2.35);
}

sub cent {

  my $stamp = shift;
  my $offset = shift;
  my $win;
  if (@_) {
    $win = shift;
  }

  $stamp = $stamp->medover;
#  $stamp -= $stamp->slice("0:10")->median;

  my $xs = $stamp->xvals;
  $xs += $offset;
  my ($cenx,$sigx) = fitgauss($xs,$stamp);

  if ($win) {
    $win->line($xs,$stamp);
    $win->hold;
    $win->line(pdl([$cenx,$cenx]),pdl([-1e6,1e6]));
    $win->release;
  }

  return($cenx,$sigx);
}

my $win = PDL::Graphics::PGPLOT::Window->new(device=>'/xs');

opendir(D,".");
#my @sfos = grep { /todfcr\.fits/ } readdir(D);
my @sfos = grep { /tr\.fits/ } readdir(D);
closedir(D);
my $slyb = 70;
my $slye = 320;

my $sl = "484:579,$slyb:$slye";

my $sfoname;
my $slx;

foreach $sfoname (@sfos) {
  my $sfo = rfits "$sfoname";
  print "$sfoname ";
  my $first = 0;
  if ($sfoname =~ m|apr22s0030|) {
    $first = 1;
  }

  $win->erase;
  $win->env($slyb,$slye,-10,100,{PlotPosition=>[0.1,.99,.2,.4]});
  my($cenx,$sigx) = cent($sfo->slice($sl),$slyb,$win);
  if ($first) {
    $slx = $cenx;
  } else {
#    $cenx -= $slx;
  }
  printf "%.2f %.2f \n",$cenx, $sigx;
#  $_ = <STDIN>;
}
