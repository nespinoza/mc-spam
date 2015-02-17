import numpy as np
import os,sys,pyfits
sys.path.append('main_codes')
sys.path.append('main_codes/transit_code')
import utils
import mc_spam

########################  OPTIONS #######################################################

# Define the output folder:
out_folder = 'results'

# Number of samples to generate with MC-SPAM (gets overriden if MCMC chains are used):
links = 100

#########################################################################################
if not os.path.exists(out_folder):
   os.mkdir(out_folder)

version = 'v.1.0'

f = open(out_folder+'/mc_spam_results.dat','w')
f.write('#             MC-SPAM '+version+'\n')
f.write('# \n')
f.write('#    Authors: Nestor Espinoza (nespino@astro.puc.cl) \n')
f.write('#             Andres Jordan   (ajordan@astro.puc.cl) \n')
f.write('# \n')
f.write('# If you use this results or parts of the code, please cite Espinoza & Jordan (2015) \n')
f.write('# \n')
f.write('# \n')
f.write('# Nomenclature: \n')
f.write('# ------------- \n')
f.write('# \n')
f.write('# (A): \t Coefficients generated using the ATLAS models.\n')
f.write('# \n')
f.write('# \n')
f.write('# (P): \t Coefficients generated using the PHOENIX models.\n')
f.write('# \n')
f.write('# ui : \t Coefficients sampled only due to error on stellar parameters (limiting coefficients).\n')
f.write('# \n')
f.write('# ui*: \t MC-SPAM coefficients, generated through fits to synthetic lightcurves.\n')
f.write('# \n')
f.write('# planet name \t u1 coeff (A) \t upper error (A) \t lower error (A) \t u2 coeff (A) \t upper error (A) \t lower error (A) \t u1* coeff (A) \t '+\
          'upper error (A) \t lower error (A) \t u2* coeff (A) \t upper error (A) \t lower error (A) \t'+
          'u1 coeff (P) \t upper error (P) \t lower error (P) \t u2 coeff (P) \t upper error (P) \t lower error (P) \t u1* coeff (P) \t '+\
          'upper error (P) \t lower error (P) \t u2* coeff (P) \t upper error (P) \t lower error (P)''\n')

print '\n'
print '\t ###############################################################'
print ''
print '\t                   MC-SPAM '+version+'\n'
print ''
print '\t         Authors: Nestor Espinoza (nespino@astro.puc.cl)'
print '\t                  Andres Jordan   (ajordan@astro.puc.cl)'
print ''
print '\t DISCLAIMER: If you use the results or parts of the code, please'
print '\t             cite Espinoza & Jordan (2015).'
print ''
print '\t > Starting MC-SPAM...'

# Extract all the stored data from the files:
pdata = utils.extract_data()
planet_names = pdata.keys()

for i in range(len(planet_names)):
	planet_name = planet_names[i]
	print '\t > Working with ',planet_name
	input_dict = {}

        # Here we create the dictionaries that will contain the MCMC samples for each parameter,
        # the estimations or the fixed coefficients:
	input_dict['mcmc'] = {}
	input_dict['estimated'] = {}
	input_dict['fixed'] = {}

        # If you have an MCMC chain, put it here as a numpy array. For example, suppose I generate an MCMC chain
        # for the planet-to-star radius ratio, p, as:
        #
        # p = np.random.normal(0.1,0.01,1e4)
        #
        # Then, to save it and use it in the MC-SPAM algorithm, you simply do:
        #
        # input_dict['mcmc']['p'] = p

        # Now, for each parameter that is not already in the 'mcmc' part of the 
        # dictionary, we save the parameters on the 'estimated' part if the errors
        # are different from zero or on the 'fixed' part if they equal zero:
        for parameter in ['p', 'i', 'aR', 'e', 'omega', 'Teff', 'logg', 'MH', 'vturb']:
            if parameter not in input_dict['mcmc'].keys():
               theta, theta_up_err, theta_down_err = pdata[planet_name][parameter]
               if theta_up_err == 0. and theta_down_err == 0.:
                  # Values of -999 are identified as values of vturb for which no measurement
                  # exist and, therefore, 2 km/s are defined as the fixed value:
                  if theta == -999:
                     theta = 2.0
                  input_dict['fixed'][parameter] = theta
               else:
                  input_dict['estimated'][parameter] = [ theta, theta_up_err, theta_down_err ]

	system = mc_spam.system(input_dict)
        system.chain_length = links
        system.verbose = False
	system.draw()
	system.get_mcspam_lds()
	print '\t > ...done, saving...'
	out_lds = system.results
        # Finally, get and save parameters:
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

f.close()
print '\t ###############################################################'
print ''
print '\t > Finished without problems! The results have been saved in the '+out_folder
print '\t   folder.'
print ''
