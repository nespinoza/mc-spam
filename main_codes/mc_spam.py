import numpy as np
import utils

from scipy.optimize import leastsq

# Some constants:
degree_to_radian = np.pi/180.0

class system:

      def __init__(self, param_dict):
          self.p = 0.0
          self.i = 0.0
          self.aR = 0.0
          self.e = 0.0
          self.omega = 0.0
          self.Teff = 0.0
          self.logg = 0.0
          self.MH = 0.0
          self.vturb = 0.0
          self.chain_length = 1e2
          self.methods = 'espinoza14'
          self.pname = 'unnamed_planet-1b'
          self.verbose = False
          self.ok = False
          self.results = None
          self.input_dict = param_dict

      def draw(self):
          param_dict = self.input_dict
          # Define the needed parameters for the MC-SPAM:

          planetary_parameter_names = ['p', 'i', 'aR','e','omega']

          planetary_parameter_limits = [ [0,None], \
                                         [0,90.0], \
                                         [0,None], \
                                         [0,1],\
                                         [0,180.0]]

          planetary_parameter_description = ['Planetary-to-star radius ratio (R_p/R_*)',\
                                             'Inclination of the orbit (in degrees)',\
                                             'Semi-major axis to stellar radius ratio (a/R_*)',\
                                             'Eccentricity of the orbit',\
                                             'Argument of periastron (degrees)']

          stellar_parameter_names = ['Teff','logg','MH','vturb']

          stellar_parameter_limits = [[0,None],\
                                      [0,None],\
                                      [None, None],\
                                      [0,None]]

          stellar_parameter_description = ['Effective temperature of the host-star (in Kelvins)',\
                                           'Log-gravity (cgs)',\
                                           'Metal content of the star (can be approximated by Fe/H)',\
                                           'Turbulent velocity (in km/s)']

          # First, check parameters that have been determined by MCMC:

          mcmc_params = param_dict['mcmc']
          chain_length = -1
          detected_parameter = False
          for parameter in mcmc_params.keys():
              if parameter in planetary_parameter_names:
                 idx = planetary_parameter_names.index(parameter)
                 planetary_parameter_description.pop(idx)
                 planetary_parameter_names.pop(idx)
                 limits = planetary_parameter_limits.pop(idx)
                 detected_parameter = True
              if parameter in stellar_parameter_names:
                 idx = stellar_parameter_names.index(parameter)
                 stellar_parameter_description.pop(idx)
                 stellar_parameter_names.pop(idx)
                 limits = stellar_parameter_limits.pop(idx)
                 detected_parameter = True
              # If a chain has been detected, save it to the current object:
              if detected_parameter:
                 exec "self."+parameter+"=mcmc_params['"+parameter+"']"
                 if chain_length == -1:
                    chain_length = len(mcmc_params[parameter])
                 detected_parameter = False

          # Now check parameters that have been determined by estimates, and simulate 
          # values by assuming the given values are the median, upper 1-sigma error and 
          # lower 1-sigma errors:

          estimated_params = param_dict['estimated']

          # If no mcmc chains supplied, set the number of samples to the default (or user defined) length:
          if chain_length == -1:
             chain_length = self.chain_length
          for parameter in estimated_params.keys():
              if parameter in planetary_parameter_names:
                 idx = planetary_parameter_names.index(parameter)
                 planetary_parameter_description.pop(idx)
                 planetary_parameter_names.pop(idx)
                 limits = planetary_parameter_limits.pop(idx)
                 detected_parameter = True
              if parameter in stellar_parameter_names:
                 idx = stellar_parameter_names.index(parameter)
                 stellar_parameter_description.pop(idx)
                 stellar_parameter_names.pop(idx)
                 limits = stellar_parameter_limits.pop(idx)
                 detected_parameter = True
              # If a parameter has been detected, simulate values for it: 
              if detected_parameter:
                 samples = utils.sample_from_errors(estimated_params[parameter][0], estimated_params[parameter][1], estimated_params[parameter][2], \
                                                    chain_length, low_lim = limits[0], up_lim = limits[1])
                 exec "self."+parameter+"= samples"
                 detected_parameter = False 

          # Now check parameters that have been fixed:
          fixed_params = param_dict['fixed']
          for parameter in fixed_params.keys():
              if parameter in planetary_parameter_names:
                 idx = planetary_parameter_names.index(parameter)
                 planetary_parameter_description.pop(idx)
                 planetary_parameter_names.pop(idx)
                 detected_parameter = True
              if parameter in stellar_parameter_names:
                 idx = stellar_parameter_names.index(parameter)
                 stellar_parameter_description.pop(idx)
                 stellar_parameter_names.pop(idx)
                 detected_parameter = True
              # If a parameter has been detected, generate the same number of copies at the 
              # chain length for it:
              if detected_parameter:
                 exec "self."+parameter+"= np.ones(chain_length)*fixed_params['"+parameter+"']"

          # Make sure all the necessary parameters have been sampled:
          if len(planetary_parameter_names) != 0 or len(stellar_parameter_names) !=0:
             print '\n'
             print '\t ERROR: MC-SPAM will not work because not all input parameters have been defined.'
             print '\t        The parameters missing from the input dictionary are: \n'
             if len(planetary_parameter_names) != 0:
                print '\t Planetary parameters:'
                print '\t --------------------\n'
                for i in range(len(planetary_parameter_names)):
                    print '\t + '+planetary_parameter_names[i]+' ('+planetary_parameter_description[i]+').'
 
             if len(stellar_parameter_names) != 0:
                print '\n'
                print '\t Stellar parameters:'
                print '\t ------------------\n'
                for i in range(len(stellar_parameter_names)):
                    print '\t + '+stellar_parameter_names[i]+' ('+stellar_parameter_description[i]+').'
          else:
                self.ok = True

          # And finally, check that all the sampled values have 0 < b < 1, where b = cos(i) a/R_* is the 
          # impact parameter (that b > 0 is guaranteed from our sampling scheme; b < 1 is not):

          # First, compute impact parameter and the number of samples that match b>=1:
          i = np.array(self.i)
          aR = np.array(self.aR)

          b = np.cos(i * degree_to_radian)*aR
          idx = np.where(b>=1)[0]
          n_idx = len(idx)
          while True:
                    # For the n_idx samples that have b>=1, recompute n_idx samples and save those to the indexes that match b>=1:
                    if n_idx > 0:
                       if 'i' in estimated_params.keys():
                          i_samps = utils.sample_from_errors(estimated_params['i'][0], estimated_params['i'][1], estimated_params['i'][2], \
                                                             n_idx, low_lim = 0, up_lim = 90.0)

                          i[idx] = i_samps

                       if 'aR' in estimated_params.keys():
                          aR_samps = utils.sample_from_errors(estimated_params['aR'][0], estimated_params['aR'][1], estimated_params['aR'][2], \
                                                              n_idx, low_lim = 0, up_lim = None)

                          aR[idx] = aR_samps
                    else:
                       break

                    # Recompute the impact parameters for the drawn samples, identify samples that have
                    # b>=1, repeat:
                    b = np.cos(i[idx] * degree_to_radian)*aR[idx]
                    idx_idx = np.where(b>=1)[0]
                    if len(idx_idx) == 0:
                       break

                    idx = idx[idx_idx]
                    n_idx = len(idx) 

          self.i, self.aR = list(i), list(aR)

      def get_mcspam_lds(self):

          def fit_quadratic_lightcurve(times, y, inc, r_a, p, e, omega, u1_guess, u2_guess):

              def residuals(parameters, data):
                  u1, u2 = parameters 
                  return data - utils.getTransit(inc, r_a, p, [u1, u2], times = times, e = e, omega = omega)

              
              guess = (u1_guess, u2_guess)
 
              plsq = leastsq(residuals, guess, args= ( y ))[0]
              return plsq

          def check_limit(value, upper_lim, lower_lim):
              if value > upper_lim:
                 value = upper_lim
              elif value < lower_lim:
                 value = lower_lim
              return value   
         
          if not self.ok:
                 print '\t ERROR: The MC-SPAM limb-darkening coefficients cannot be obtained without samples from'
                 print '\t        the parameters. Run the draw() function in order to draw them.'
          else:
                 # Define variables that will save results. First, get the functions that interpolate the 
                 # limb-darkening coefficients and, using them, define the stellar models for which
                 # limb-darkening coefficients will be obtained:
                    
                 interpolators = utils.read_and_interpolate_lds(methods = self.methods)
                 stellar_models = interpolators.keys()
 
                 # Now define the variables that will save the MC-SPAM limb 
                 # darkening coefficients for the given stellar models, and the 
                 # associated dictionary that will be outputed:
 
                 out_lds = {}
                 
                 for stellar_model in stellar_models:
                     out_lds[stellar_model] = {}
                     out_lds[stellar_model]['u1*'] = []
                     out_lds[stellar_model]['u2*'] = []
                     out_lds[stellar_model]['u1'] = []
                     out_lds[stellar_model]['u2'] = []

                 # Get how many times we will repeat the algorithm:
                 length = len(self.p)

                 # Start the algorithm:
                 for i in range(length):
                     # Get stellar parameters of the given draw:
                     teff_d = self.Teff[i]
                     logg_d = self.logg[i]
                     mh_d = self.MH[i]
                     vturb_d = self.vturb[i]
                     for stellar_model in stellar_models:
                         # Check that none of the stellar parameters is out of the limits of the interpolators of the given
                         # stellar model. If someone is, set it to the closest value:
                         low_teff, low_logg, low_mh, low_vturb = interpolators[stellar_model]['lower limits']
                         up_teff, up_logg, up_mh, up_vturb = interpolators[stellar_model]['upper limits']
                         c_teff = check_limit(np.copy(teff_d), up_teff, low_teff)
                         c_logg = check_limit(np.copy(logg_d), up_logg, low_logg)
                         c_mh = check_limit(np.copy(mh_d), up_mh, low_mh)
                         c_vturb = check_limit(np.copy(vturb_d), up_vturb, low_vturb)

                         # Now, get non-linear and quadratic (as input guess for the least squares problem) limb-darkening 
                         # coefficients given those stellar parameters:
                         point = np.array([c_teff, c_logg, c_mh, c_vturb])
                         c1, c2, c3, c4 = interpolators[stellar_model]['c1'](point), interpolators[stellar_model]['c2'](point), \
                                          interpolators[stellar_model]['c3'](point), interpolators[stellar_model]['c4'](point)

                         # The guesses are the "limiting coefficients":
                         u1_guess, u2_guess = (12./35.)*c1 + c2 + (164./105.)*c3 + 2.*c4, (10./21.)*c1 - (34./63.)*c3 - c4


                         if self.verbose:
                            print '\t ---------- \n \t Generating transit for:'
                            print '\t i     : ',self.i[i]
                            print '\t r_a   : ',1./self.aR[i]
                            print '\t (b    : ',np.cos(self.i[i]*degree_to_radian)*self.aR[i],')'
                            print '\t p     : ',self.p[i]
                            print '\t e     : ',self.e[i]
                            print '\t omega : ',self.omega[i]

                         # Generate a transit lightcurve given the current parameters:
                         times,lightcurve = utils.getTransit(self.i[i]*degree_to_radian, 1./self.aR[i], self.p[i], [c1, c2, c3, c4], \
                                                       e = self.e[i], omega = self.omega[i], ld_law = 'non-linear')


                         # Obtain the SPAM coefficients of the given draw:
                         u1_d, u2_d = fit_quadratic_lightcurve(times, lightcurve, self.i[i]*degree_to_radian, 1./self.aR[i], self.p[i], self.e[i], \
                                                               self.omega[i], u1_guess, u2_guess)

                         #if self.verbose:
                         #   from matplotlib.pyplot import plot,legend,show
                         #   plot(times, lightcurve,label = 'Original lightcurve')
                         #   plot(times, utils.getTransit(self.i[i]*degree_to_radian, 1./self.aR[i], self.p[i], [u1_d, u2_d], times = times, \
                         #        e = self.e[i], omega = self.omega[i]),label='Fit')
                         #   legend()
                         #   show()
                         # Save them:
                         out_lds[stellar_model]['u1*'].append(np.copy(u1_d))
                         out_lds[stellar_model]['u2*'].append(np.copy(u2_d))
                         # Save the drawn model LDs:
                         out_lds[stellar_model]['u1'].append(np.copy(u1_guess))
                         out_lds[stellar_model]['u2'].append(np.copy(u2_guess))

                 # Finally, convert the lists to numpy arrays and save the output:
                 for stellar_model in stellar_models:
                     out_lds[stellar_model]['u1'] = np.array(out_lds[stellar_model]['u1'])
                     out_lds[stellar_model]['u2'] = np.array(out_lds[stellar_model]['u2'])
                     out_lds[stellar_model]['u1*'] = np.array(out_lds[stellar_model]['u1*'])
                     out_lds[stellar_model]['u2*'] = np.array(out_lds[stellar_model]['u2*'])

                 self.results = out_lds

      def save(self, out_folder = 'results'):

                 # If folder out_folder does not exist, create it:
                 if not os.path.isdir(out_folder):
                    os.mkdir( out_folder )

                 

                
                 out_lds = self.results
		 # Get and save parameters:
		 # First of ATLAS results:
		 u1_atlas,u1_atlas_up_error,u1_atlas_low_err = utils.getParams(out_lds['atlas']['u1'])
		 u2_atlas,u2_atlas_up_error,u2_atlas_low_err = utils.getParams(out_lds['atlas']['u2'])
		 u1_mcspam_atlas,u1_mcspam_atlas_up_error,u1_mcspam_atlas_low_err = utils.getParams(out_lds['atlas']['u1*'])
		 u2_mcspam_atlas,u2_mcspam_atlas_up_error,u2_mcspam_atlas_low_err = utils.getParams(out_lds['atlas']['u2*'])
		 # Now for PHOENIX results:
		 u1_phoenix,u1_phoenix_up_error,u1_phoenix_low_err = utils.getParams(out_lds['phoenix']['u1'])
		 u2_phoenix,u2_phoenix_up_error,u2_phoenix_low_err = utils.getParams(out_lds['phoenix']['u2'])
		 u1_mcspam_phoenix,u1_mcspam_phoenix_up_error,u1_mcspam_phoenix_low_err = utils.getParams(out_lds['phoenix']['u1*'])
		 u2_mcspam_phoenix,u2_mcspam_phoenix_up_error,u2_mcspam_phoenix_low_err = utils.getParams(out_lds['phoenix']['u2*'])

		 # First, save the coefficients and the correspondig errors:
		 f.write(planet_name+'\t'+str(u1_atlas)+'\t'+str(u1_atlas_up_error)+'\t'+str(u1_atlas_low_err)+'\t'+\
			 str(u2_atlas)+'\t'+str(u2_atlas_up_error)+'\t'+str(u2_atlas_low_err)+'\t'+\
			 str(u1_mcspam_atlas)+'\t'+str(u1_mcspam_atlas_up_error)+'\t'+str(u1_mcspam_atlas_low_err)+'\t'+\
			 str(u2_mcspam_atlas)+'\t'+str(u2_mcspam_atlas_up_error)+'\t'+str(u2_mcspam_atlas_low_err)+'\t'+\
			 str(u1_phoenix)+'\t'+str(u1_phoenix_up_error)+'\t'+str(u1_phoenix_low_err)+'\t'+\
			 str(u2_phoenix)+'\t'+str(u2_phoenix_up_error)+'\t'+str(u2_phoenix_low_err)+'\t'+\
			 str(u1_mcspam_phoenix)+'\t'+str(u1_mcspam_phoenix_up_error)+'\t'+str(u1_mcspam_phoenix_low_err)+'\t'+\
			 str(u2_mcspam_phoenix)+'\t'+str(u2_mcspam_phoenix_up_error)+'\t'+str(u2_mcspam_phoenix_low_err)+'\n')

		 # Now save the distributions of the limb-darkening coefficients:
		 if not os.path.isdir(out_folder+'/'+planet_name):
		    os.mkdir(out_folder+'/'+planet_name)
		 pyfits.PrimaryHDU(out_lds['atlas']['u1']).writeto(out_folder+'/'+planet_name+'/u1_atlas.fits')
		 pyfits.PrimaryHDU(out_lds['atlas']['u2']).writeto(out_folder+'/'+planet_name+'/u2_atlas.fits')
		 pyfits.PrimaryHDU(out_lds['atlas']['u1*']).writeto(out_folder+'/'+planet_name+'/u1_mcspam_atlas.fits')
		 pyfits.PrimaryHDU(out_lds['atlas']['u2*']).writeto(out_folder+'/'+planet_name+'/u2_mcspam_atlas.fits')
		 pyfits.PrimaryHDU(out_lds['phoenix']['u1']).writeto(out_folder+'/'+planet_name+'/u1_phoenix.fits')
		 pyfits.PrimaryHDU(out_lds['phoenix']['u2']).writeto(out_folder+'/'+planet_name+'/u2_phoenix.fits')
		 pyfits.PrimaryHDU(out_lds['phoenix']['u1*']).writeto(out_folder+'/'+planet_name+'/u1_mcspam_phoenix.fits')
		 pyfits.PrimaryHDU(out_lds['phoenix']['u2*']).writeto(out_folder+'/'+planet_name+'/u2_mcspam_phoenix.fits')
