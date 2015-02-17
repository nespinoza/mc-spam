
#
# GENERATED WITH PDL::PP! Don't modify!
#
package PDL::Planet;

@EXPORT_OK  = qw( PDL::PP RV PDL::PP true_anomaly PDL::PP Orbit PDL::PP transit_times );
%EXPORT_TAGS = (Func=>[@EXPORT_OK]);

use PDL::Core;
use PDL::Exporter;
use DynaLoader;



   
   @ISA    = ( 'PDL::Exporter','DynaLoader' );
   push @PDL::Core::PP, __PACKAGE__;
   bootstrap PDL::Planet ;








=head1 FUNCTIONS



=cut






=head2 RV

=for sig

  Signature: ( double t(); double v0(); double K(); double w(); double e(); 
double t0(); double period(); double [o] rv();)


=for ref

info not available




=cut






*RV = \&PDL::RV;




=head2 true_anomaly

=for sig

  Signature: (double t(); double t0(); double e(); double period(); 
double [o] ta();)


=for ref

info not available




=cut






*true_anomaly = \&PDL::true_anomaly;




=head2 Orbit

=for sig

  Signature: (double t(); double t0(); double period(); double a(); double e(); double w(); double I(); double Omega(); double [o] xx(); double [o] yy(); double [o] zz())


=for ref

info not available




=cut






*Orbit = \&PDL::Orbit;




=head2 transit_times

=for sig

  Signature: (double e(); double w(); double t0(); double period(); double [o] t1(); double [o] t2())


=for ref

info not available




=cut






*transit_times = \&PDL::transit_times;


;



# Exit with OK status

1;

		   