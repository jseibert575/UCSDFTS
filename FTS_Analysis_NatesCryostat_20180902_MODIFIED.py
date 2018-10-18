####################################################################################################
#
# FTS Analysis Code for UCSD FTS with Nates Cryostat
#
# Lindsay Ng Lowry
# llowry@ucsd.edu
# 20180904
#
# This program was written by Lindsay to analyze data taken with the UCSD FTS using Nate's Cryostat.
# The program uses functions from the FTS_Analysis_KEK.py code used in KEK during Run38 testing of
# PB2a (written by Fred Matsuda), but in a very simplified way because this data is taken using only
# a single bolometer in the cryostat and without rotating the output polarizer.
# 
# The input data is simply a pickle file keyed by linear stage position whose items are arrays of
# bolometer data from the lock-in amplifier (with the chopper chopping at 18 Hz).
#
####################################################################################################



# Import necessary modules/functions ===============================================================
import os
import sys
import numpy as np
import matplotlib
import pylab as pl
import cPickle as pickle
import argparse
from scipy.optimize import leastsq, fmin_bfgs, minimize, curve_fit
import time

import lib_FTSAnalysis_NatesCryostat_20180902 as lFA

def getargs():
	"""Parse the input arguments"""
	parser = argparse.ArgumentParser(description='Argument Parser for analysis of FTS Data with Nates Cryostat')

	parser.add_argument('-i', dest='inputfn', action='store', default = None,
						help='Specify input data directory')
	parser.add_argument('-d', dest='datafn', action='store', default=None,
						help='Data file (pickle file for single bolo keyed by linear stage position)')
	parser.add_argument('-o', dest='outputdir', action='store', default='.',
						help='Specify output directory')
	parser.add_argument('-inttime', dest='inttime', action='store', default=2.0, type=float,
						help='Integration time for FTS step')

	args = parser.parse_args()

	return args



def get_interferogram_pts(data):
	"""Get the signal vs. frequency points from the raw data

	Params:
	data (dict) - the dictionary contained in the data file (keyed by linear stage position)
	
	Return:
	interferogram_pts (array) - an array where index [0] is an array of lin stage positions, [1] is
								an array of corresponding mean signal values, and [2] is an array of
								corresponding standard deviations on the mean signal value 
	"""

	#make lists to hold the position, signal, and std values
	posArray = []
	signalArray = []
	signalStdArray = []

	#for each position, calculated the mean and std of the signal and add these to the lists
	for pos in sorted(data.keys()):
		posArray.append(pos)

		mean_signal = np.mean(data[pos])
		std_signal = np.std(data[pos])

		signalArray.append(mean_signal)
		signalStdArray.append(std_signal)

	#convert to np arrays
	posArray = np.array(posArray)
	signalArray = np.array(signalArray)
	signalStdArray = np.array(signalStdArray)

	#returned the combined array
	interferogram_pts = np.array([posArray, signalArray, signalStdArray])
	return interferogram_pts


def analyze_spectra_NatesCryostat(interferogram_data, output_path, amode, steps_per_inch = 125000.):
	"""Main function to get an analyze the spectra from the measured data

	Parameters:
	interferogram_data (array) - array containing the lin stage positions and corresponding signal
								 values and standard deviations (errors) as produced by
								 get_interferogram_pts
	output_path (str) - the path to where the output files/plots should be saved
	amode (str) - 'd' if double-sided interferogram, 's' if single-sided
	steps_per_inch (float) - encoder steps per inch on the linear stage, default 125000.

	Return:
	dataout - dictionary containing the interferogram data and the spectral data
	"""

	# Some setup things ============================================================================
	#make the directory for the output plots and save a plot of zeros to test
	lFA.make_plots_dir(output_path)
	lFA.plot_blank(output_path)

	#make the dataout dictionary and add the interferogram data
	dataout = {}
	dataout['Interferograms'] = interferogram_data

	#get the position and signal arrays form interferogram_data
	positions = interferogram_data[0]
	interf_signals = interferogram_data[1]

	#determine the stepsize
	steps = np.diff(positions)
	mean_step = np.mean(steps)
	stepsize = mean_step/steps_per_inch
	std_step = np.std(steps)

	print 'Mean step size (in encoder counts): ' + str(mean_step)
	print 'Mean step size (in inches): ' + str(stepsize)
	if std_step != 0.0:
		print 'Stepsize not consistent.  Standard deviation is ' + str(std_step) + '.  May want to examine data.'
	else:
		print 'Stepsize is consistent.  Continuing.'

	# Some data verification things ================================================================
	#check that the 'energy' is sufficient after detrending (FROM FRED'S CODE - USE NOT TOTALLY CLEAR YET)
	energy = np.std(lFA.detrend(interf_signals))
	#if energy < 0.25:
	#	print 'Energy from detrending is less than 0.25.  Exiting.'
	#	sys.exit()
	#check that there is a maximum in the interferogram data
	max_data = np.max(interf_signals)
	if np.isnan(max_data):
		print 'NaN maximum in interferogram data.  Exiting.'
		sys.exit()



	# Get the spectral data using lFA.fts_spectrum() ===============================================
	spectral_data = lFA.fts_spectrum(xs_raw=interf_signals, name='Nate Cryostat Data',
									 window_name='triangle', plot_path=output_path,
									 stepsize = stepsize, mode=amode)


	# Plot the spectrum ============================================================================
	lFA.plot_spectrum(spectral_data,'o-')
	lFA.plot_conf()
	pl.savefig(os.path.join(output_path,'plots/fts-points.png'))
	pl.clf()


	# Add the spectral data to the output dictionary and return the dictionary =====================
	dataout['Spectra'] = spectral_data
	return dataout



def main(args):
	"""Run the FTS Analysis code using the parsed arguments"""

	#print args
	#print args.datafn
	
	#load the data file
	if args.datafn == None:
		print 'No data file specified.  Exiting analysis script.'
		sys.exit()
	full_input_filename = os.path.join(args.inputfn, args.datafn)
	data = lFA.load_dict(full_input_filename)
	print 'Data file loaded.'

	#get the interferogram points from the raw data
	interferogram_data = get_interferogram_pts(data)

	#run the analyze_spectra function
	outputdata = analyze_spectra_NatesCryostat(interferogram_data, output_path = args.outputdir,
											   amode = 'd')

	pl.close()

	#save the output data
	print 'Saving pkl file for raw data outputs...'
	
	datetime = args.datafn[:15]
	with open(os.path.join(args.outputdir,'FTS_data-%s.pkl'%(datetime)),'w') as f:
		pickle.dump(outputdata,f)

	print '\nDONE\n'



if __name__ == "__main__":
	#parse the arguments and run the main code
	print 'Starting FTS analysis'

	args = getargs()
	main(args)

