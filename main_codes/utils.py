import numpy as np
import sys
import interpolator

def read_and_interpolate_lds(methods = 'espinoza14'):
    """
    Description
    -----------

    Function made to read limb-darkening coefficients tables and create a linear interpolator. 
    The idea is that, given a file with the limb-darkening coefficients with differents temperatures,  
    logg, M/H and vturb, this function returns an interpolation function, that is, a function f(Teff,logg,M/H,vturb,law) 
    that returns the linearly interpolated coefficients of a given law at a given Teff, logg, M/H, vturb.

    Inputs
    ------   

	methods		Method from which limb darkening coefficients will be obtained. These values will then be interpolated.


    Outputs
    -------

	interpolator	A dictionary containing interpolator objects as well as upper and lower limits of those interpolations.

    """

    interpolators = {}
    interpolators['atlas'] = {}
    interpolators['phoenix'] = {}
    stellar_models = ['atlas','phoenix']
    for stellar_model in stellar_models:
        if methods == 'espinoza14':
           Teff,logg,mh,vturb,a,u1,u2,b1,b2,b3,c1,c2,c3,c4,l1,l2,e1,e2,s1,s2 = \
                np.loadtxt('data/tabulated_lds/'+stellar_model+'_kepler_lds.dat',unpack=True,\
                            usecols=(3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22))
        elif methods == 'claret':
           if stellar_model == 'phoenix':
              logg1,teff1,mh1,vturb1,u1_1,u2_1 = np.loadtxt('data/tabulated_lds/claret/phoenix_kepler_lowteff_claret_quadratic.dat',unpack=True,usecols=(0,1,2,3,9,10))
              logg2,teff2,mh2,vturb2,u1_2,u2_2 = np.loadtxt('data/tabulated_lds/claret/phoenix_kepler_highteff_claret_quadratic.dat',unpack=True,usecols=(0,1,2,3,9,10))
              logg = np.append(logg1,logg2)
              Teff = np.append(teff1,teff2)
              mh = np.append(mh1,mh2)
              vturb = np.append(vturb1,vturb2)
              u1 = np.append(u1_1,u1_2)
              u2 = np.append(u2_1,u2_2)
              logg1,teff1,mh1,vturb1,c1_1,c2_1,c3_1,c4_1 = np.loadtxt('data/tabulated_lds/claret/phoenix_kepler_lowteff_claret_nonlinear.dat',unpack=True,usecols=(0,1,2,3,4,5,6,7))
              logg2,teff2,mh2,vturb2,c1_2,c2_2,c3_2,c4_2 = np.loadtxt('data/tabulated_lds/claret/phoenix_kepler_highteff_claret_nonlinear.dat',unpack=True,usecols=(0,1,2,3,4,5,6,7))
              c1 = np.append(c1_1,c1_2)
              c2 = np.append(c2_1,c2_2)
              c3 = np.append(c3_1,c3_2)
              c4 = np.append(c4_1,c4_2)
           else:
              logg,Teff,mh,vturb,u1,u2 = np.loadtxt('data/tabulated_lds/claret/atlas_kepler_claret_quadratic.dat',unpack=True,usecols=(0,1,2,3,4,5))
              logg,teff,mh,vturb,c1,c2,c3,c4 = np.loadtxt('data/tabulated_lds/claret/atlas_kepler_claret_nonlinear.dat',unpack=True,usecols=(0,1,2,3,4,5,6,7))

        # Define lower and upper limits of the grid:
        lo = np.array([np.min(Teff), np.min(logg), np.min(mh), np.min(vturb)])
        hi = np.array([np.max(Teff), np.max(logg), np.max(mh), np.max(vturb)])

        # Generate the grid of data:
        teff_vals, teff_idx = np.unique(Teff, return_inverse = True)
        logg_vals, logg_idx = np.unique(logg, return_inverse = True)
        mh_vals, mh_idx = np.unique(mh, return_inverse=True)
        vturb_vals, vturb_idx = np.unique(vturb, return_inverse=True)
        u1_array = np.empty(teff_vals.shape + logg_vals.shape + mh_vals.shape + vturb_vals.shape)
        u2_array = np.empty(teff_vals.shape + logg_vals.shape + mh_vals.shape + vturb_vals.shape)
        c1_array = np.empty(teff_vals.shape + logg_vals.shape + mh_vals.shape + vturb_vals.shape)
        c2_array = np.empty(teff_vals.shape + logg_vals.shape + mh_vals.shape + vturb_vals.shape)
        c3_array = np.empty(teff_vals.shape + logg_vals.shape + mh_vals.shape + vturb_vals.shape)
        c4_array = np.empty(teff_vals.shape + logg_vals.shape + mh_vals.shape + vturb_vals.shape)

        u1_array[teff_idx, logg_idx, mh_idx, vturb_idx] = u1
        u2_array[teff_idx, logg_idx, mh_idx, vturb_idx] = u2
        c1_array[teff_idx, logg_idx, mh_idx, vturb_idx] = c1
        c2_array[teff_idx, logg_idx, mh_idx, vturb_idx] = c2
        c3_array[teff_idx, logg_idx, mh_idx, vturb_idx] = c3
        c4_array[teff_idx, logg_idx, mh_idx, vturb_idx] = c4

        # Now, use the Intergrid function to interpolate values:
        maps = [teff_vals, logg_vals, mh_vals, vturb_vals]
        interfunc_u1 = interpolator.Intergrid( u1_array, lo=lo, hi=hi, maps=maps, verbose=0 )
        interfunc_u2 = interpolator.Intergrid( u2_array, lo=lo, hi=hi, maps=maps, verbose=0 )
        interfunc_c1 = interpolator.Intergrid( c1_array, lo=lo, hi=hi, maps=maps, verbose=0 )
        interfunc_c2 = interpolator.Intergrid( c2_array, lo=lo, hi=hi, maps=maps, verbose=0 )
        interfunc_c3 = interpolator.Intergrid( c3_array, lo=lo, hi=hi, maps=maps, verbose=0 )
        interfunc_c4 = interpolator.Intergrid( c4_array, lo=lo, hi=hi, maps=maps, verbose=0 )


        # Save the interpolators:
        interpolators[stellar_model]['u1'] = interpolator.Intergrid( u1_array, lo=lo, hi=hi, maps=maps, verbose=0 )
        interpolators[stellar_model]['u2'] = interpolator.Intergrid( u2_array, lo=lo, hi=hi, maps=maps, verbose=0 )
        interpolators[stellar_model]['c1'] = interpolator.Intergrid( c1_array, lo=lo, hi=hi, maps=maps, verbose=0 )
        interpolators[stellar_model]['c2'] = interpolator.Intergrid( c2_array, lo=lo, hi=hi, maps=maps, verbose=0 )
        interpolators[stellar_model]['c3'] = interpolator.Intergrid( c3_array, lo=lo, hi=hi, maps=maps, verbose=0 )
        interpolators[stellar_model]['c4'] = interpolator.Intergrid( c4_array, lo=lo, hi=hi, maps=maps, verbose=0 )

        # Save limits of the interpolators:
        interpolators[stellar_model]['lower limits'] = np.copy(lo)
        interpolators[stellar_model]['upper limits'] = np.copy(hi)

    return interpolators

def sample_from_errors(loc, sigma1, sigma2, n, low_lim = None, up_lim = None):
    """
    Description
    -----------

    Function made to sample points given the 0.16, 0.5 and 0.84 quantiles of a parameter

    In the case of unequal variances, this algorithm assumes a skew-normal distribution and samples from it. 
    If the variances are equal, it samples from a normal distribution.

    Inputs
    ------

          loc:           Location parameter (we hope it is the 0.5 quantile, i.e., the median).

       sigma1:           Upper error of the parameter (we hope loc+sigma1 is the 0.84 quantile, i.e., the "upper 1-sigma bound").
 
       sigma2:           Lower error of the parameter (we hope loc-sigma2 is the 0.16 quantile, i.e., the "lower 1-sigma bound").

            n:           Number of samples you want to generate from this.

      low_lim:           (Optional) Lower limits on the values of the samples*.

       up_lim:           (Optional) Upper limits on the values of the samples*.

    Outputs
    -------

        The output are n samples from the distribution that best-matches the quantiles.

    *The optional inputs (low_lim and up_lim) are lower and upper limits that the samples have to have; if any of the samples 
    surpasses those limits, new samples are drawn until no samples do. Note that this changes the actual variances of the samples.

    """

    if (sigma1 != sigma2):
            """
             If errors are assymetric, sample from a skew-normal distribution given
             the location parameter (assumed to be the median), sigma1 and sigma2.
            """

            # First, find the parameters mu, sigma and alpha of the skew-normal distribution that 
            # best matches the observed quantiles:
            sknorm = skew_normal()
            sknorm.fit(loc,sigma1,sigma2)

            # And now sample n values from the distribution:
            samples = sknorm.sample(n)

            # If a lower limit or an upper limit is given, then search if any of the samples surpass 
            # those limits, and sample again until no sample surpasses those limits:
            if low_lim is not None:
               while True:
                     idx = np.where(samples<low_lim)[0]
                     l_idx = len(idx)
                     if l_idx > 0:
                        samples[idx] = sknorm.sample(l_idx)
                     else:
                        break
            if up_lim is not None:
               while True:
                     idx = np.where(samples>up_lim)[0]
                     l_idx = len(idx)
                     if l_idx > 0:
                        samples[idx] = sknorm.sample(l_idx)
                     else:
                        break
            return samples       

    else:
            """
             If errors are symmetric, sample from a gaussian
            """
            samples = np.random.normal(loc,sigma1,n)
            # If a lower limit or an upper limit is given, then search if any of the samples surpass 
            # those limits, and sample again until no sample surpasses those limits:
            if low_lim is not None:
               while True:
                     idx = np.where(samples<low_lim)[0]
                     l_idx = len(idx)
                     if l_idx > 0:
                        samples[idx] = np.random.normal(loc,sigma1,l_idx)
                     else:
                        break
            if up_lim is not None:
               while True:
                     idx = np.where(samples>up_lim)[0]
                     l_idx = len(idx)
                     if l_idx > 0:
                        samples[idx] = np.random.normal(loc,sigma1,l_idx)
                     else:
                        break
            return samples
           
from scipy.integrate import quad
from scipy.optimize import leastsq

class skew_normal: 
      """
      Description
      -----------

      This class defines a skew_normal object, which generates a skew_normal distribution given the quantiles 
      from which you can then sample datapoints from.

      """
      def __init__(self):
         self.mu = 0.0
         self.sigma = 0.0
         self.alpha = 0.0

      def fit(self, median, sigma1, sigma2):
            """
            This function fits a Skew Normal distribution given
            the median, upper error bars (sigma1) and lower error bar (sigma2).
            """

            # First, define the sign of alpha, which should be positive if right skewed
            # and negative if left skewed:
            alpha_sign = np.sign(sigma1-sigma2)

            # Now define the residuals of the least-squares problem:
            def residuals(p, data, x):
                mu, sqrt_sigma, sqrt_alpha = p
                return data - model(x, mu, sqrt_sigma, sqrt_alpha)

            # Define the model used in the residuals:
            def model(x, mu, sqrt_sigma, sqrt_alpha):
                """
                Note that we pass the square-root of the scale (sigma) and shape (alpha) parameters, 
                in order to define the sign of the former to be positive and of the latter to be fixed given
                the values of sigma1 and sigma2:
                """
                return self.cdf (x, mu, sqrt_sigma**2, alpha_sign * sqrt_alpha**2)

            # Define the quantiles:
            y = np.array([0.15866, 0.5, 0.84134])

            # Define the values at which we observe the quantiles:
            x = np.array([median - sigma2, median, median + sigma1])

            # Start assuming that mu = median, sigma = mean of the observed sigmas, and alpha = 0 (i.e., start from a gaussian):
            guess = (median, np.sqrt ( 0.5 * (sigma1 + sigma2) ), 0)
            # Perform the non-linear least-squares optimization:
            plsq = leastsq(residuals, guess, args=(y, x))[0]

            self.mu, self.sigma, self.alpha = plsq[0], plsq[1]**2, alpha_sign*plsq[2]**2

      def sample(self, n):
            """
            This function samples n points from a skew normal distribution using the 
            method outlined by Azzalini here: http://azzalini.stat.unipd.it/SN/faq-r.html.
            """
            # Define delta:
            delta = self.alpha/np.sqrt(1+self.alpha**2)
            
            # Now sample u0,u1 having marginal distribution ~N(0,1) with correlation delta:
            u0 = np.random.normal(0,1,n)
            v = np.random.normal(0,1,n)
            u1 = delta*u0 + np.sqrt(1-delta**2)*v

            # Now, u1 will be random numbers sampled from skew-normal if the corresponding values 
            # for which u0 are shifted in sign. To do this, we check the values for which u0 is negative:
            idx_negative = np.where(u0<0)[0]
            u1[idx_negative] = -u1[idx_negative]

            # Finally, we change the location and scale of the generated random-numbers and return the samples:
            return self.mu + self.sigma*u1

      @staticmethod
      def cdf(x, mu, sigma, alpha):
            """
            This function simply calculates the CDF at x given the parameters 
            mu, sigma and alpha of a Skew-Normal distribution. It takes values or 
            arrays as inputs.
            """
            if type(x) is np.ndarray:
               out = np.zeros(len(x))
               for i in range(len(x)):
                   out[i] = quad(lambda x: skew_normal.pdf(x,mu,sigma,alpha), -np.inf, x[i])[0]
               return out

            else:
               return quad(lambda x: skew_normal.pdf(x,mu,sigma,alpha), -np.inf, x)[0]

      @staticmethod
      def pdf(x, mu, sigma, alpha):
            """
            This function returns the value of the Skew Normal PDF at x, given
            mu, sigma and alpha
            """
            def erf(x):
                # save the sign of x
                sign = np.sign(x)
                x = abs(x)

                # constants
                a1 =  0.254829592
                a2 = -0.284496736
                a3 =  1.421413741
                a4 = -1.453152027
                a5 =  1.061405429
                p  =  0.3275911

                # A&S formula 7.1.26
                t = 1.0/(1.0 + p*x)
                y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*np.exp(-x*x)
                return sign*y

            def palpha(y,alpha):
                phi = np.exp(-y**2./2.0)/np.sqrt(2.0*np.pi)
                PHI = ( erf(y*alpha/np.sqrt(2)) + 1.0 )*0.5
                return 2*phi*PHI

            return palpha((x-mu)/sigma,alpha)*(1./sigma)

#############################################################################################

import sys
sys.path.append('transit_code')
from Transit import Transit, Transit_nl

def get_transit_duration(P, p, r_a, inclination, eccentricity):
    """
    Description
    -----------

    This function returns the maximum transit duration given the period (P), planet-to-star radius ratio (p = Rp/R_*), 
    the stellar-to-semi-major axis ratio (r_a = R_*/a), the eccentricity e, and the inclination in radians.

    """


    b = np.cos(inclination)/r_a
    num1 = ( 1. + p )**2 - b**2
    den1 = 1. - np.cos(inclination)**2

    return (P/np.pi)*np.arcsin(r_a * np.sqrt(num1/den1))*np.sqrt((1.+eccentricity)/(1.-eccentricity))

def getTransit(inclination, r_a, p, ld_coeffs, times = None, e=0.0, omega=0.0, ld_law= 'quadratic', npoints = 1000):
    """
    Description
    -----------
    This function returns npoints values of the Mandel & Agol (2002) lightcurve of a transiting planet assuming unitary 
    Period (i.e., sampling in phase space). It samples 400 points out-of-transit and 1000 points in-transit. 


    Inputs
    ------

    inclination:  Inclination in radians.

    r_a:          Stellar radius over semi-major axis (R_*/a)

    p:            Planet-to-star radius ratio.

    ld_coeffs:    Array with limb darkening coefficients (assumed quadratic, but can be changed with the ld_law parameter).

    times:        (optional) array of times.

    e:            (optional) Eccentricity.

    omega:        (optional) Argument of perihastron.

    ld_law:       (optional) String which defines the limb-darkening law to be used.

    npoints:      (optional) Number of in-transit points.

    Outputs
    -------

    If the times are given, this code returns the relative flux at those times. If, not, it returns both the times and 
    the relative fluxes.

    """

    Period = 1.0
    t0 = 0.0
    given_times = True
    if times is None:
       # Calculate duration of the transit:
       transit_time = get_transit_duration(Period, p, r_a, inclination,e)
       # Generate times based on this duration:
       times = np.linspace(-(transit_time)/2.0,(transit_time)/2.0,npoints)
       # Add two hundred points before and after transit, just to have some points off-transit:
       delta_times = np.diff(times)[0]
       time_points_before = times[0]-(np.arange(1,201,1)*delta_times)
       time_points_after = times[-1]+(np.arange(1,201,1)*delta_times)
       times = np.append( time_points_before ,times )
       times = np.append( times, time_points_after )
       given_times = False

    if ld_law == 'quadratic':
       trans,trans00=Transit(times, Period, inclination, r_a, p, t0, ld_coeffs[0], ld_coeffs[1], e, omega)
    elif ld_law == 'non-linear':
       trans,trans00=Transit_nl(times, Period, inclination, r_a, p, t0, ld_coeffs[0], ld_coeffs[1], ld_coeffs[2], ld_coeffs[3], e, omega)
    if given_times:
       return trans
    else:
       return times,trans

def getParams(dist,alpha = 0.68):
    """
    Description
    -----------
    This function returns, given samples from a distribution, the median, the upper "sigma" error 
    and lower sigma error (i.e., the distance to the 0.16 and 0.84 quantiles).

    Inputs
    ------
 
    dist		Samples from a distribution.

    alpha		% of the confidence bands around the median that you want to sample.    

    Outputs
    -------

    The value at the median, lower "error" and upper "error"

    """
    ordered_dist = dist[np.argsort(dist)]
    param = 0.0
    nsamples = len(dist) # Number of samples from posterior
    nsamples_at_each_side = int(nsamples*(alpha/2.)+1)
    med_idx = 0
    if(nsamples%2 == 0.0): # Number of points is even
       med_idx_up = int(nsamples/2.)
       med_idx_down = med_idx_up-1
       param = (ordered_dist[med_idx_up]+ordered_dist[med_idx_down])/2.
       return param,np.abs(ordered_dist[med_idx_up+nsamples_at_each_side]-param),np.abs(param-ordered_dist[med_idx_down-nsamples_at_each_side])
    else:
       med_idx = int(nsamples/2.)
       param = ordered_dist[med_idx]
       return param,np.abs(ordered_dist[med_idx+nsamples_at_each_side]-param),np.abs(param-ordered_dist[med_idx-nsamples_at_each_side])

def extract_data(pdata_fname = "data/estimated_parameters/planet_data.dat", sdata_fname = "data/estimated_parameters/star_data.dat"):

    # First, extract data. Stellar data first, from the files:
    f = open(sdata_fname,'r')
    stellar_data = np.array([])
    counter = 0
    while True:
       line = f.readline()
       if line == '':
          break
       elif line[0] != '#':
          line = line.split('\n')[0]
          if counter == 0:
             stellar_data = np.array(line.split('\t'))
             counter = 1
          else:
             stellar_data = np.vstack((stellar_data,np.array(line.split('\t'))))
    f.close()

    # Now get planet data:
    f = open(pdata_fname,'r')
    planet_data = np.array([])
    counter = 0
    while True:
       line = f.readline()
       if line == '':
	  break
       elif line[0] != '#':
	  line = line.split('\n')[0]
	  if counter == 0:
	     planet_data = np.array(line.split('\t'))
	     counter = 1
	  else:
	     planet_data = np.vstack((planet_data,np.array(line.split('\t'))))
    f.close()

    # Now save all planetary and stellar data in dictionaries:
    pdata = {}
    sdata = {}

    # Get names of all the planets:
    planet_names = planet_data[:,0]

    # Get names of the parameters and their indexes in the planet_data.dat file:
    pparams_names = ['p', 'i', 'aR','e','omega']
    pparams_idx = [1, 4, 7, 10, 13]

    # Same for star_data:
    sparams_names = ['Teff', 'logg', 'MH','vturb']
    sparams_idx = [1, 4, 7, 10]

    # Iterate through all the planet names:
    for i in range(len(planet_names)):
        pname = planet_names[i].strip()
        params = {}
        for j in range(len(pparams_idx)):
            c_idx = pparams_idx[j]
            c_pparam_name = pparams_names[j]
            params[ c_pparam_name ] = [ np.double(planet_data[i,c_idx]), np.double(planet_data[i,c_idx+1]), np.double(planet_data[i,c_idx+2])  ]

        # Search which star is the host of the current planet:
        stellar_host_data = []
        for j in range(len(stellar_data[:,0])):
            if stellar_data[j,0].strip() in pname:   
               stellar_host_data = stellar_data[j,:]
               break

        if len(stellar_host_data) == 0:
           print '\t ERROR: No stellar host data in the star_data.dat file for planet '+pname+'.'
           sys.exit()

        for j in range(len(sparams_idx)):
            c_idx = sparams_idx[j]
            c_sparam_name = sparams_names[j]
            params[ c_sparam_name ] = [ np.double(stellar_host_data[c_idx]), np.double(stellar_host_data[c_idx+1]), np.double(stellar_host_data[c_idx+2])  ]

        pdata[pname] = params
        
    return pdata
