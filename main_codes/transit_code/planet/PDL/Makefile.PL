
# Use this as a template for the Makefile.PL for
# any external PDL module.

use ExtUtils::MakeMaker;

use PDL::Core::Dev;
@pack = (["pdl_planet.pd",Planet,PDL::Planet]);

%hash = pdlpp_stdargs(@::pack);
# $hash{'OPTIMIZE'}  = '-g'; # If you want to debug, uncomment this.
# $hash{INC} .= " -I/usr/local/include";  # uncomment as required
#$hash{INC} .= " -I/scisoft/i386/Packages/gsl-1.9/include ";
#$hash{MYEXTLIB} .= " /Users/ajordan/devel/planet/libplanet.a ";
$hash{LIBS}[0] .= " -L/home/ajordan/devel/planet/ -lplanet -lgsl -lgslcblas -lm "; # uncomment as required

WriteMakefile(%hash);

# Add genpp rule
# add other makefile additions as required (see also ExtUtils::MakeMaker)
sub MY::postamble {
	pdlpp_postamble(@::pack);
}

