# - Library used with FTS_Analysis_PB2b_20190221.py
# - Adapted from the code sent to Lindsay by Fred that was used in Run38 testing at KEK



#!/usr/bin/env python


#from PyPolarbear.Mapping import *
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as pl
import scipy as sp
from scipy.optimize import leastsq, fmin_bfgs, minimize, curve_fit
import scipy.integrate as integ
import cPickle as pickle
import argparse
import os
import sys
import subprocess
import heapq
#import AnalysisBackend.misc.parse_flat_map as parsemap
#from AnalysisBackend.mapping import libpolmap
#from AnalysisBackend.unpack.datafolder import DataFolder
#from AnalysisBackend.deglitching import packetdrop
#from AnalysisBackend.misc import util
#from AnalysisBackend.misc.skyside_focalplane import gen_focalplane_offsets


def load_dict(fn):
	"""Loads a dictionary from a pickle file given a filename"""
	f = open(fn,'r')
	d = pickle.load(f)
	f.close()
	return d

    
def show_keys(d):
	"""Prints the keys in a dictionary and the datatype of each"""
	print "Keys:"
	keys = d.keys()
	keys.sort()
	for k in keys: print k,type(d[k])


""" CANT USE BECAUSE DEPENDS ON ANALYSIS BACKEND ----------------------------------------------------------------------------
def load_all(datafn,timepkl,chan=None):
	df = DataFolder(datafn)
	f = open(timepkl,'r')
	dt = pickle.load(f)
	f.close()
	if chan is not None:
		times = dt[chan+'_timestamp']
	else:
		#keys = dt.keys()
		#keys.remove('Position')
		#times = dt[keys[0]]
		times = dt
	return df,times
"""

""" PROBABLY DONT NEED THIS BECAUSE DIDNT TAKE TIME STAMPS ------------------------------------------------------------------
def get_step_times(p,deltat=0.72879752):
	#calculate the keys into Modified Julian Days
	#jarray is the times in one integration time step
	#starttimes stores the earliest time in each time step
	#endtime stores the latest time in each time step
	starttimes = []
	endtimes = []
	#offset = 56657. + 0.72879752/3600./24.
	print 'Adjusting step times directly from QuickRecord with time correction %s sec '%(deltat)
	offset = 56657. + deltat/3600./24.
    	for i in range(len(p)):
		#print ' Calculating step %s/%s'%(i,len(p)-1)
        	jarray = []
        	for j in range(len(p[i])):
            		p2 = p[i][j]
			days = p2['days']
			secdays = p2['seconds']*(1./60.)*(1./60.)*(1./24.)
			nsdays = p2['10nsTicks']*(1e-8)*(1./3600.)*(1./24.)
			hundays = p2['hundredsOfDays']*100.
			hrdays = p2['hours']*(1./24.)
			yrdays = p2['years']*365.
			mindays = p2['minutes']*(1./60.)*(1./24.)
			Jdays = days+secdays+nsdays+hundays+hrdays+yrdays+mindays+offset
			jarray.append(Jdays)
        	starttimes.append(min(jarray))
        	endtimes.append(max(jarray))
    	return jarray,starttimes,endtimes
"""

""" PROBABLY DONT NEED THIS BECAUSE DIDNT TAKE TIME STAMPS ------------------------------------------------------------------
def closest_value(df,starttimes,endtimes):
    	#c1 and c2 are used in order to determine the closest values to timestamp in bolotime
    	#c1 and c2 determine which arguments in bolotime are closest
	bolo_times = df.bolo_time
    	closestart = []
    	closeend = []
    	for i in range(len(starttimes)):
		c1 = abs(bolo_times-starttimes[i])
		c2 = abs(bolo_times-endtimes[i])
        	closestart.append(bolo_times[np.argmin(c1)])
        	closeend.append(bolo_times[np.argmin(c2)])
    	return closestart,closeend
"""

def step_idx(df,starttimes,endtimes,buff=0.,steplength=2.0):
	"""Create array of indices of where stepping and integrating
	
	Parameters:
	df () - ???
	starttimes () - ???
	endtimes () - ???
	buff (float) - ???, default 0.
	steplength (float) - ???, default 2.0

	Return:
	
	"""
	sample_rate = 25e6/2.**17.
	deltat = steplength
	extra = 1.0
	nt = np.round((deltat+extra)*sample_rate)
	steps = np.zeros((len(starttimes),int(nt)),dtype=int)
	steps_mask = np.zeros((len(starttimes),int(nt)),dtype=bool)
	bolo_times = df.bolo_time
	#print len(bolo_times), steps.shape
	for i in range(len(starttimes)):
		c1 = np.argmin(abs(bolo_times-starttimes[i]-buff/3600./24.))
		c2 = np.argmin(abs(bolo_times-endtimes[i]+buff/3600./24.))
		#print starttimes[i],c1,endtimes[i],c2
		steps[i][:c2+1-c1] = range(c1,c2+1)
		steps_mask[i][:c2+1-c1] = True
	return steps, steps_mask


""" DONT NEED THIS NOW =============================================================================
def load_drops(df):
	# Load arrays that contain where packet drops occurred so we can fill in lost data
	print 'Loading packet drop data...'
	mobotimes = df.load('bolo_time_all')
	mobodrops = packetdrop.find_all_drops(mobotimes)
	mobomasks = packetdrop.build_mobo_masks(mobotimes,mobodrops)
	return mobomasks
"""

""" DONT NEED THIS NOW =============================================================================
def find_packet_skips(times):
	# Determine packet skip times which is when data is completely lost in the timestream
	sample_rate = 25e6/2.**17.
	times_shifted = np.zeros(len(times))
	times_shifted[1:] = times[:-1]
	times_shifted[0] = times[0]
	timebefore = (times - times_shifted)*24.*3600.	# sec 
	skips = [timebefore>2.*1./sample_rate]
	skip_idx = np.where(skips[0])[0]
	skip_num = np.round(timebefore[skips[0]]*sample_rate)-1.
	print 'Found %s packet skips with total of %s samples lost.'%(len(skip_idx),np.sum(skip_num))
	return skips,skip_idx,skip_num
"""

""" DONT NEED THIS NOW =============================================================================
def calc_refIQ(refsig,skip_idx,skip_num,harm=1,deltat=0.72879752):
	# Calculate in-phase and quadrature components from chopper reference signal in order to demodulate signal
	print 'Adjusting quad pair for timing error of value %s sec '%(deltat)
	skip_tot = np.sum(skip_num)
	if skip_tot > 0.:
		print 'Adjusting quad pair for packet skips...'
		refIQ_ext = GenQuadPair2(refsig,harm=harm,extra=int(skip_tot),deltat=deltat)
		n = 0
		refIQ_cut = np.ones(len(refIQ_ext['X']),dtype=bool)
		for i in range(len(skip_idx)):
			refIQ_cut[skip_idx[i]+n:skip_idx[i]+int(skip_num[i])+n] = False 
			n += int(skip_num[i])
		refx = refIQ_ext['X'][refIQ_cut]
		refy = refIQ_ext['Y'][refIQ_cut]
		refIQ = {'X':refx,'Y':refy}
	else:
		refIQ = GenQuadPair2(refsig,harm=harm,deltat=deltat)
	return refIQ
"""

""" DONT NEED THIS NOW =============================================================================
def calc_refIQ_value(bolotime,skip_idx,skip_num,harm=1,deltat=0.72879752):
    # Calculate in-phase and quadrature components from chopper reference signal in order to demodulate signal
    print 'Adjusting quad pair for timing error of value %s sec '%(deltat)
	refIQ = GenQuadPairValue(bolotime,harm=harm,deltat=deltat)
    return refIQ
"""

""" DONT NEED THIS NOW =============================================================================
def calc_reftiming(bolotime,reftime,skip_idx,skip_num):
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION
	skip_tot = np.sum(skip_num)
    if skip_tot > 0.:
		n = 0
        cut = np.ones(len(reftime),dtype=bool)
        for i in range(len(skip_idx)):
            cut[skip_idx[i]+n:skip_idx[i]+int(skip_num[i])+n] = False
            n += int(skip_num[i])
		reftime_cut = reftime[cut]
		bolotime_cut = bolotime[:-int(np.sum(skip_num))]
	else:
		reftime_cut = np.copy(reftime)
		bolotime_cut = np.copy(bolotime)
	delta_t = np.average(bolotime_cut-reftime_cut)
	#delta_t *= -1.
	print 'Timing difference of bolo_time - reference_time = %.6f sec'%(delta_t*24.*3600.)
	return delta_t*24.*3600.
"""

""" DONT NEED THIS NOW =============================================================================
def demod_data(df,refsig,refIQ,closestart,closeend,mobomasks,buff=0.,harm=1,wafer=None,chan_itr=None,chans_ref=None,hwmapfn=None,save_dir=None,phase_corr=False,phase_opt=None,verbose=True,quick=False):
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION
	# Demodulate data and do phase corrections and etc 
	data = np.zeros((1512,len(closestart)))
	noise = np.zeros((1512,len(closestart)))
	drops = np.zeros((1512,len(closestart)),dtype=bool)
	phase_ang = np.zeros(1512)
	phase_ang[:] = np.NaN
	Xout = np.zeros((1512,len(closestart)))
	Yout = np.zeros((1512,len(closestart)))

    chans = range(1512)
	xmlmap = parsemap.build_index_maps(hwmapfn)
	out_chlist, xpos_fp, ypos_fp, polangles, crosspol, boloids = gen_focalplane_offsets(xmlmap,cal_source=None)
	if wafer is not None:
		wafers = ['1','2','3','4','5','6','7']
		wafernames = ['10.2','10.4','10.3','10.5','9.4','10.1','8.2.0']
		wafernums = np.array(xmlmap.wafer)
		wchs = np.where(wafernums==str(wafer))[0]
		chans = np.intersect1d(chans,wchs)
		#chans = [1503]#[1026,1030,1503,1507]#[1099,1103,1082,1086]##[1325,1329,1343,1347]#[2,3,6,7]
		if quick:
			print 'Demodulation in Quick mode with only ~90 channels'
			chans = chans[::2]
	elif (chan_itr is not None) and (phase_opt is not None):
		chans = [int(chan_itr)]
		phase_ang[int(chan_itr)] = phase_opt
	if verbose:
		print chans
	
	if (chan_itr is not None) and (chans_ref is not None):
		print 'Cannot use iteration and reference option simultaneously. Exiting...'
		sys.exit(1)
	elif (chans_ref is not None) and (phase_opt is not None):
		for chan_ref in chans_ref:
			delta = polangles - polangles[chan_ref]
			phase_ang[delta==0.] = phase_opt[chan_ref]
			#phase_ang = phase_opt - delta*np.pi/180.

	bolo_times = df.bolo_time
	#refIQ = GenQuadPair2(refsig,harm=harm) 

	for chan in chans:
	if verbose:
		print 'Demodulating channel %s'%(chan)

	mobo,moboch = util.chantomobo(int(chan))
	#print mobo
	#drops_idx = mobodrops[mobo]
	pktmask = mobomasks[mobo]

	#drop_count = 0
	drop_step = 0
	partial_step = 0
	ts = df.get_bolo(chan)
	#print len(ts), len(pktmask)
	if save_dir is not None:
		pl.figure()
		pl.plot(ts[::50])
		pl.savefig(os.path.join(save_dir,'plots/ts_%s.png'%(xmlmap.boloid[chan])))
		pl.close()

	#print len(bolo_times), len(ts), len(refsig), len(refIQ['X'])
	'''
	time_temp = bolo_times[(bolo_times>closestart[149])&(bolo_times<closeend[151])]
	pl.figure()
	pl.plot((time_temp-time_temp[0])*24.*3600.)
	pl.show()
	pl.clf()
	data_temp = ts[(bolo_times>closestart[149])&(bolo_times<closeend[151])]
	ref_temp = refsig[(bolo_times>closestart[149])&(bolo_times<closeend[151])]
	refX_temp = refIQ['X'][(bolo_times>closestart[149])&(bolo_times<closeend[151])]*1000.
	refY_temp = refIQ['Y'][(bolo_times>closestart[149])&(bolo_times<closeend[151])]*1000.
	mask_temp = pktmask[(bolo_times>closestart[149])&(bolo_times<closeend[151])]
	pl.plot(mask_temp)
	pl.show()
	pl.clf()
	pl.plot((data_temp-np.average(data_temp))*10.)
	pl.show()
	pl.clf()
	pl.plot(ref_temp)
	pl.show()
	pl.clf()
	pl.plot(refX_temp)
	pl.show()
	pl.clf()
	pl.plot(refY_temp)
	pl.show()
	pl.clf()
	pl.plot((data_temp-np.average(data_temp))*10.)
	pl.plot(ref_temp)
	pl.plot(refX_temp)
	pl.plot(refY_temp)
	pl.show()
	pl.close()
	'''

	for i in range(len(closestart)):
		#print closestart, closeend
		#print bolo_times
		start = closestart[i]-buff/3600./24.
		end = closeend[i]+buff/3600./24.
		data_step = ts[(bolo_times>start)&(bolo_times<end)]
		#ref_step = refsig[(bolo_times>start)&(bolo_times<end)]
		refIQ_step = {'X':refIQ['X'][(bolo_times>start)&(bolo_times<end)],'Y':refIQ['Y'][(bolo_times>start)&(bolo_times<end)]}
		mask_step = pktmask[(bolo_times>start)&(bolo_times<end)]
		#print i, start, end, len(data_step), len(ref_step)
		#nzeros = len(data_step[data_step==0.])
		#drop_data = [data_step!=0.][0]
		nzeros = len(data_step[-mask_step])

		#print nzeros
		sample_rate = 25e6/2.**17.
		step_time = (closeend[i]-closestart[i])*24.*3600.
		#print step_time
		if len(data_step) == 0: print 'bolo_time missing for this step.'
		if (len(data_step) == 0) or (nzeros > len(data_step)/2.) or (len(data_step) < int(step_time*sample_rate/2.)):
			drop_step += 1
			drops[chan,i] = True
			#data[chan,i] = np.average(data[chan,:i-1])
			#noise[chan,i] = np.average(noise[chan,:i-1])
			#print data[chan,i], noise[chan,i]
			data[chan,i] = 0.
			noise[chan,i] = 0.
			Xout[chan,i] = 0.
			Yout[chan,i] = 0.		
		elif not np.all(mask_step):
			partial_step += 1
			#print 'Found partial signal bin. Demod with partial data.'
			#data_demod = Demod(data_step,ref_step,harm=harm,drop_data=mask_step,phase_corr=phase_corr)#,phase_opt=phase_ang[chan])
			data_demod = Demod2(data_step,refIQ_step,drop_data=mask_step,phase_corr=phase_corr,phase_opt=phase_ang[chan])
			data[chan,i] = data_demod['Mag']
			noise[chan,i] = data_demod['Noise']
			Xout[chan,i] = data_demod['X']
			Yout[chan,i] = data_demod['Y']
		else:
			#data_demod = Demod(data_step,ref_step,harm=harm,phase_corr=phase_corr)#,phase_opt=phase_ang[chan])
			data_demod = Demod2(data_step,refIQ_step,phase_corr=phase_corr,phase_opt=phase_ang[chan])
			data[chan,i] = data_demod['Mag']
			noise[chan,i] = data_demod['Noise']
			Xout[chan,i] = data_demod['X']
			Yout[chan,i] = data_demod['Y']
		#print i, len(data_step),data_demod['X'],data_demod['Y'],data_demod['Mag']
		#pl.figure()
		#pl.plot(data_step)
		#pl.plot(refIQ_step['X'])
		#pl.plot(refIQ_step['Y'])
		#pl.show()
		#pl.close()

	if verbose:
		print 'Dropped signal bins %i / %i'%(drop_step,len(closestart))
		print 'Partial signal bins %i / %i'%(partial_step,len(closestart))
	if save_dir is not None:
		pl.figure()
		pl.plot(Xout[chan])
		pl.plot(Yout[chan])
		#pl.show()
		pl.savefig(os.path.join(save_dir,'plots/demodXY_%s.png'%(xmlmap.boloid[chan])))
		pl.close()
	return data,noise,chans,drops,Xout,Yout
"""

""" DONT NEED THIS NOW =============================================================================
def demod_data2(df,step_idx,step_mask,refIQ,mobomasks,buff=0.,harm=1,wafer=None,chan_itr=None,chans_ref=None,hwmapfn=None,save_dir=None,phase_corr=False,phase_opt=None,verbose=True,quick=False):
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION
	# Demodulate data and do phase corrections and etc 
	data = np.zeros((1512,step_idx.shape[0]))
	noise = np.zeros((1512,step_idx.shape[0]))
	drops = np.zeros((1512,step_idx.shape[0]),dtype=bool)
	Xout = np.zeros((1512,step_idx.shape[0]))
	Yout = np.zeros((1512,step_idx.shape[0]))	
	phase_ang = np.zeros(1512)
	phase_ang[:] = np.NaN

    chans = range(1512)
	xmlmap = parsemap.build_index_maps(hwmapfn)
	if wafer is not None:
		wafers = ['1','2','3','4','5','6','7']
		wafernames = ['10.2','10.4','10.3','10.5','9.4','10.1','8.2.0']
		wafernums = np.array(xmlmap.wafer)
		wchs = np.where(wafernums==str(wafer))[0]
		chans = np.intersect1d(chans,wchs)
		#chans = [1503]#[1026,1030,1503,1507]#[1099,1103,1082,1086]##[1325,1329,1343,1347]#[2,3,6,7]
		if quick:
			print 'Demodulation in Quick mode with only ~90 channels'
			chans = chans[::2]
	elif (chan_itr is not None) and (phase_opt is not None):
		chans = [int(chan_itr)]
		phase_ang[int(chan_itr)] = phase_opt
	if verbose:
		print chans
	
	if (chan_itr is not None) and (chans_ref is not None):
		print 'Cannot use iteration and reference option simultaneously. Exiting...'
		sys.exit(1)
	elif (chans_ref is not None) and (phase_opt is not None):
		out_chlist, xpos_fp, ypos_fp, polangles, crosspol, boloids = gen_focalplane_offsets(xmlmap,cal_source=None)
		for chan_ref in chans_ref:
			delta = polangles - polangles[chan_ref]
			phase_ang[delta==0.] = phase_opt[chan_ref]

	#chans = [957]

	bolo_times = df.bolo_time

	if save_dir is not None:
		pl.figure()
    	for chan in chans:
		if verbose:
			print 'Demodulating channel %s'%(chan)

		mobo,moboch = util.chantomobo(int(chan))
		pktmask = mobomasks[mobo]
		ts = df.get_bolo(chan)

		data_step = ts[step_idx]
		pktmask_step = pktmask[step_idx]
		refIQ_step = {'X':refIQ['X'][step_idx],'Y':refIQ['Y'][step_idx]} 
		#print data_step.shape,data_step
		#print pktmask_step.shape,pktmask_step
		#print refIQ_step['X'].shape,refIQ_step['X']
		#print refIQ_step['Y'].shape,refIQ_step['Y']

		if save_dir is not None:
			pl.plot(ts[::50])
			pl.savefig(os.path.join(save_dir,'plots/ts_%s.png'%(xmlmap.boloid[chan])))
			pl.clf()

		mask = step_mask & pktmask_step
		#print mask.shape, mask
		data_demod = Demod_chan(data_step,refIQ_step,mask,phase_opt=phase_ang[chan])
		
		data[chan][:] = data_demod['Mag']
		noise[chan][:] = data_demod['Noise']
		Xout[chan][:] = data_demod['X']
		Yout[chan][:] = data_demod['Y']

		ndata = np.sum(mask,axis=1)
		nzeros = np.sum(-mask,axis=1)
		ndata_pkt = np.sum(pktmask_step,axis=1)
		drops[chan][ndata<step_mask.shape[1]/2.] = True
		drop_step = np.sum(drops[chan])
		partial_step = np.sum((ndata_pkt>step_mask.shape[1]/2.)&(ndata_pkt<step_mask.shape[1]))

		if verbose:
			print 'Dropped signal bins %i / %i'%(drop_step,step_idx.shape[0])
			print 'Partial signal bins %i / %i'%(partial_step,step_idx.shape[0])
		if save_dir is not None:
			pl.plot(Xout[chan])
			pl.plot(Yout[chan])
			#pl.show()
			pl.savefig(os.path.join(save_dir,'plots/demodXY_%s.png'%(xmlmap.boloid[chan])))
			pl.clf()
	#if save_dir is not None:
	#	pl.close()

	return data,noise,chans,drops,Xout,Yout
"""

""" DONT NEED THIS NOW =============================================================================
def phase_iteration(df,refsig,refIQ,closestart,closeend,mobomasks,harm=1,chan_itr=None,hwmapfn=None):
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION
	# Calculate necessary phase correction using iteration procedure and maximization of interferogram STD
	def resid_phase(c):
		data, noise, chans, drops, X, Y = demod_data(df,refsig,refIQ,closestart,closeend,mobomasks,buff=0.,harm=harm,chan_itr=chan_itr,hwmapfn=hwmapfn,phase_opt=c,verbose=False)
		data_chan = data[chan_itr]
		#pl.figure()
		#pl.plot(data_chan)
		#pl.show()
		#pl.close()
		cut = -drops[chan_itr]
		#data_full, noise_full = interp_data(chans,data,noise,drops,verbose=False)
		data_d = detrend(data_chan[cut],poly=5)
		var = np.std(data_d)
		print var, 
		sys.stdout.flush()
		return -var
	print ' \n'
	p0 = 0.
    p = minimize(resid_phase,p0)
	print 'Optimized phase = %.3f \n'%(p.x[0]*180./np.pi)
	return p.x[0]
"""

""" DONT NEED THIS NOW =============================================================================
def phase_iteration2(df,step_idx,step_mask,refIQ,mobomasks,harm=1,chan_itr=None,hwmapfn=None):
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION
	# Calculate necessary phase correction using iteration procedure and maximization of interferogram STD
	def resid_phase(c):
		data, noise, chans, drops, X, Y = demod_data2(df,step_idx,step_mask,refIQ,mobomasks,buff=0.,harm=harm,chan_itr=chan_itr,hwmapfn=hwmapfn,phase_opt=c,verbose=False)
		data_chan = data[chan_itr]
		cut = -drops[chan_itr]
		data_d = detrend(data_chan[cut],poly=5)
		var = np.std(data_d)
		print var, 
		sys.stdout.flush()
		return -var
	print ' \n'
	p0 = 0.
    p = minimize(resid_phase,p0)
	print 'Optimized phase = %.3f \n'%(p.x[0]*180./np.pi)
	return p.x[0]
"""

""" DONT NEED THIS NOW =============================================================================
def phase_iteration2_minY(df,step_idx,step_mask,refIQ,mobomasks,harm=1,chan_itr=None,hwmapfn=None):
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION    
    # Calculate necessary phase correction using iteration procedure and maximization of interferogram STD
    def resid_phase(c):
        data, noise, chans, drops, X, Y = demod_data2(df,step_idx,step_mask,refIQ,mobomasks,buff=0.,harm=harm,chan_itr=chan_itr,hwmapfn=hwmapfn,phase_opt=c,verbose=False)
        Y_chan = Y[chan_itr]
        cut = -drops[chan_itr]
        Y_d = detrend(Y_chan[cut],poly=5)
		#Y_d = Y_chan[cut]
        return Y_d
    p0 = 0.
    p = leastsq(resid_phase,p0)
    print 'Optimized phase = %.3f \n'%(p[0][0]*180./np.pi)
    return p[0][0]
"""

""" DONT NEED THIS NOW =============================================================================
def interp_data(chans,data,noise,drops,verbose=True):
	# Interpolate and fill in data where lost
	for chan in chans:

		if np.all(drops[chan]):
			if verbose:
				print 'Channel %s is not functioning.'%(chan)
			continue

		elif (not np.all(drops[chan])) and np.any(drops[chan]):
			indrop = False
			drop_sec = []
			for i in range(len(drops[chan])):
				if (not indrop) and drops[chan,i]:
					indrop = True
					dropstart = i-1		# Last sample before the gap
				elif indrop and (not drops[chan,i]):
					indrop = False
					dropend = i			# First sample after the gap
					drop_sec.append((dropstart,dropend))
			if indrop:
				drop_sec.append((dropstart,len(drops[chan])))	

			n = len(data[chan])
			for start,end in drop_sec:
				if start < 0 and end < n:
					data[chan,start+1:end-1] = data[chan,end]
				elif start >= 0 and end >= n:
					data[chan,start+1:n] = data[chan,start]
				elif start >= 0 and end < n:
					p = np.polyfit([0,end-start],[data[chan,start],data[chan,end]],1)
					ii = np.arange(1,end-start)
					fit = np.polyval(p,ii)
					data[chan,start+1:end] = fit
			
		else:
			if verbose:
				print 'Channel %s has no packet drops'%(chan)
			continue

	return data,noise
"""
 
def est_zeropath(xs):
	"""Find the index of the interferogram data corresponding to zero path difference
	(and therefore max signal - white light fringe)

	Parameters:
	xs (array) - the interferogram signal data

	Return:
	the index where the signal is maximum
	"""
	return xs.argmax()

""" DONT NEED THIS NOW =============================================================================
def plot_phase(fs,xs,label,path):	
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION
	x = xs.real
	y = xs.imag
   	theta = np.arctan2(y,x)*180./np.pi
	#print max(theta), min(theta), np.average(theta)

	#pl.figure()
	pl.plot(fs,theta,'-',label=label)
	pl.title('Frequency vs Phase')
	pl.savefig(os.path.join(path,'plots/phase-%s.png'%(label)))
	pl.ylim([-45,45])
	pl.xlim([100,200])
	pl.savefig(os.path.join(path,'plots/phase-%s_zoom.png'%(label)))
	#pl.show()
	pl.clf()
	
	#print 'imag band AVG: ',np.average(y[(fs>100.)&(fs<200.)])
	#print 'imag band MAX: ',np.amax(abs(y[(fs>100.)&(fs<200.)]))
	#print 'imag band STD: ',np.std(y[(fs>100.)&(fs<200.)])
	#print 'imag out AVG,STD: ',np.average(y[fs>200.]),np.std(y[fs>200.])
	pl.plot(fs,y,'-')
	pl.savefig(os.path.join(path,'plots/phase-%s_imag.png'%(label)))
	pl.clf()
	#pl.close()
"""

def zeropath_search(xs):
	#?????????????????????? DONT FULLY UNDERSTAND YET ??????????????????????????????????????????????
	"""Find the optimal zero path point by minimizing the energy in the discrete sine transform
	Maybe maximizing zero path is a better heuristic? (Fred's comment)

	Steps: (Fred's comments)
	Take symmetric chunk of data preadjusted to zero estimate
	Fourier transform
	Phase shift until optimally lined up - minimized imag^2 
	
	Parameters:
	xs (array) - the interferogram points, centered around the zero path difference estimate

	Return:
	result from the leastsq minimization
    """
	fxs = np.fft.fft(xs)
    #print len(fxs)

	def phase_shift(xs,c):
		n = np.fft.fftfreq(len(xs))
		return np.exp(2.*np.pi*1j*c*n)*xs

	def resid(c):
		return np.imag(phase_shift(fxs,c))

	p0 = 0.0 #starting estimate for leastsq minimization
	p,ier = leastsq(resid,p0)

	return p

""" DONT NEED THIS NOW =============================================================================
def zeropath_search2(xs):
	'''Maximize zero path to find phase shift'''
	fxs = np.fft.fft(xs)
    #print len(fxs)
	
	def phase_shift2(x,c):
        n = np.fft.fftfreq(len(x))
        return np.exp(2.*np.pi*1j*c*n)

	def resid2(c):	
		filter = np.fft.ifft(phase_shift2(fxs,c))
		filter = np.fft.ifftshift(filter)
		window = np.hamming(len(fxs))
    	xs_shifted = np.convolve(xs,filter.real*window,mode='valid')
		return -1.*xs_shifted[0]

	p0 = 0.0
    p = fmin_bfgs(resid2,p0)

	return p[0]
"""

def zeropath_check(xs,stepsize,fl=50.,fr=150.):#fl=100.,fr=200.):
	#?????????????????????? DONT FULLY UNDERSTAND YET ??????????????????????????????????????????????
	"""Check that the zeroed interferogram produces a reasonable spectrum
	Parameters:
	xs (array) - interferogram data centered around the zero path point
	stepsize (float) - inches per step
	fl (float) - frequency band lower limit, default 50.
	fr (float) - frequency band upper limit, default 150.

	Return:
	True if the standard deviation of the inband data is less than 3*the standard deviation of the
	out-of-band data
	"""
	ys = np.fft.fft(xs)
	fs = np.fft.fftfreq(len(ys))
	y = ys.imag

	c = 11.8 #inches/ns
	inches_per_step = stepsize
	fs /= 2*inches_per_step/c

	#inband_avg = np.average(y[(fs>100.)&(fs<200.)])
	inband_max = np.amax(abs(y[(fs>fl)&(fs<fr)]))
	inband_std = np.std(y[(fs>fl)&(fs<fr)])
	outband_std = np.std(y[fs>fr])

	return (inband_std<3.*outband_std)

    
def shift_filter(delay):
	#?????????????????????? DONT FULLY UNDERSTAND YET ??????????????????????????????????????????????
	"""
	Parameters:
	delay (float) - 
	Return:

	"""
	taps = 16
	f = np.fft.fftfreq(taps)
	phase_shift = np.exp(2.*np.pi*1j*delay*f)

	filter = np.fft.ifft(phase_shift)
	filter = np.fft.ifftshift(filter)
	window = np.hamming(taps)

	def showfilt():
		pl.plot(filter.real*window)
		pl.plot(filter.imag*window)
		pl.show()
		#pl.close()

	return filter.real*window


def detrend(xs,poly=5,params=False):
	"""Remove an overall polynomial trend from the data, assuming equal spacing
	Uses the np.polyfit() function.

	Parameters:
	xs (array) - the data to be fit
	poly (int) - the degree of the polynomial with which to fit
	params (bool) - if True, return the residual and the fit, if False only return the resisual

	Return:
	xsd (array) - the original array minus the fit
	fit (array) - the polynomial fit (only if params = True)
	"""
	
	ns = np.arange(len(xs))
	
	#fit the data assuming equal spacing and get the residual
	fit = np.polyfit(ns,xs,poly)
	xsd = xs - np.polyval(fit,ns)

	#return the residual and/or the fit
	if params:
		return xsd, fit
	else:
		return xsd


""" DONT NEED THIS NOW =============================================================================
def get_array(lst,key):
	return np.array([x[key] for x in lst])
"""

""" DONT NEED THIS NOW ============================================================================= 
def single_phase_lockin(lst):
 	mags = get_array(lst,'Mag')
 	phis = get_array(lst,'Phi')
 	x = get_array(lst,'X')
 	y = get_array(lst,'Y')

 	maxi = mags.argmax()
 	phi_max = phis[maxi]*np.pi/180.0
 	locked_in = x*np.cos(phi_max) + y*np.sin(phi_max)

 	return locked_in
"""   

""" DONT NEED THIS NOW ============================================================================= 
def find(lst):
	for i in range(len(lst)):
    	if lst[i]: return i
"""

""" DONT NEED THIS NOW =============================================================================
def rfind(lst):
	for i in range(len(lst))[::-1]:
    	if lst[i]: return i
"""

""" DONT NEED THIS NOW =============================================================================
def est_band(fs,xs):
	threshold = 0.3
	select = xs>threshold
	l = find(select)-1
	r = rfind(select)+1
	integ_bandwidth = np.trapz(xs[l:r],fs[l:r])
	center_eff = np.trapz(fs[l:r]*xs[l:r],fs[l:r])/integ_bandwidth

	return integ_bandwidth,center_eff,l,r
"""

""" DONT NEED THIS NOW =============================================================================
def est_band2(fs,xs,atmdata,threshold=0.05,fmin_l=100.,fmin_r=150):
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION
	selecthigh = np.where(xs>threshold)[0]
	l = selecthigh[0]
	li = 0
	while fs[l] < fmin_l:
		l = selecthigh[li]
		if ((li+1) < (len(selecthigh)-1)):
			li += 1
		else:
			break 
	selectlow = np.where(xs<threshold)[0]
	#try:
	r = selectlow[0]
	#except:
	#	r = l+10
	ri = 0
	while fs[r] < fmin_r:
		r = selectlow[ri]
		if ((ri+1) < (len(selectlow)-1)):
            ri += 1
        else:
            break
	#print l, r
	integ_bandwidth = np.trapz(xs[l:r],fs[l:r])
	center_eff = np.trapz(fs[l:r]*xs[l:r],fs[l:r])/integ_bandwidth

	atm_fs = atmdata['Frequencies']
	atm_t = atmdata['Transmission']
	atm_temp = atmdata['Temperature']
	atm_temp_interp = np.interp(fs,atm_fs,atm_temp)
	temp_eff = np.trapz(atm_temp_interp[l:r]*xs[l:r],fs[l:r])/integ_bandwidth
	#print temp_eff

	return integ_bandwidth,center_eff,temp_eff,l,r
"""

""" DONT NEED THIS NOW =============================================================================
def calc_gain(data1,data2,atmdata):
	temp1 = data1['eff_temperature']
	bw1 = data1['integ_bandwidth']
	#print temp1*bw1

	xs2 = data2['throughput']
	fs2 = data2['frequencies']

	def resid_temp(c):
		try:
			bw2,cen2,temp2,l2,r2=est_band2(fs2,c*xs2,atmdata)
			#print temp2*bw2
			#print temp1*bw1-temp2*bw2
			return temp1*bw1-temp2*bw2
		except:
			print 'est_band2 Failed'
			return 1e3
			

	p0 = 1.0
	p,ier = leastsq(resid_temp,p0)

	return p[0]
"""

def calc_sn(xs,ezp=None):
	"""Calculate the S/N of interferogram data using the peak as the signal and the std of data
	points away from the peak as the noise

	Parameters:
	xs (array) - the data
	ezp (int) - the index of the peak, default None

	Return:
	sn (float) - the S/N
	"""
	if ezp is not None:
		peak = xs[ezp]
		tail = xs[2*ezp-100:2*ezp]
	else:
		peak = np.amax(xs)
		tail = xs[-110:-10]
	noise = np.std(tail)
	sn = peak/noise
	print 'S/N: ', sn
	return sn

""" DONT NEED THIS NOW =============================================================================
def diff_band(data1,data2,path,pix_mode=True,full_output=False,plotname=None):
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION
	window1 = data1['window']
    xs1 = data1['throughput']
    fs1 = data1['frequencies']
    name1 = data1['name']
    poly1 = data1['polynomial']
	l1 = data1['left_edge']
    r1 = data1['right_edge']
	bw1 = data1['integ_bandwidth']
	cen1 = data1['eff_frequency']
	sn1 = data1['sig2noise']
	slope1 = data1['slope']
	
	window2 = data2['window']
    xs2 = data2['throughput']
    fs2 = data2['frequencies']
    name2 = data2['name']
    poly2 = data2['polynomial']
	l2 = data2['left_edge']
    r2 = data2['right_edge']
    bw2 = data2['integ_bandwidth']
    cen2 = data2['eff_frequency']
    sn2 = data2['sig2noise']
    slope2 = data2['slope']

	xs2_interp = np.interp(fs1,fs2,xs2)
	diff_std = np.std((xs1-xs2_interp)[l1:r1])
	diff_avg = np.average((xs1-xs2_interp)[l1:r1])

	#pl.figure()
	pl.plot(fs1,xs1-xs2_interp,'-')
	pl.title(name1+'-'+name2)
	pl.ylim([-0.15,0.25])
	pl.figtext(0.6,0.75,'In-band STD = %.3f'%(diff_std))
	pl.figtext(0.6,0.7,'In-band AVG = %.3f'%(diff_avg))
	#pl.show()
	if pix_mode:
		pl.savefig(os.path.join(path,'plots/ftsdiff-%s.png'%(name1[:-1])))
	elif plotname is not None:
		pl.savefig(os.path.join(path,'plots/ftsdiff-%s.png'%(name1+'_'+plotname)))
	else:
		pl.savefig(os.path.join(path,'plots/ftspairdiff-'+name1+'-'+name2+'.png'))
	pl.clf()
	#pl.close()

	delta_l = l1-l2
	delta_r = r1-r2
	delta_bw = (bw1-bw2)/bw1
	delta_cen = (cen1-cen2)/cen1
	delta_sn = sn1-sn2
	delta_sl = slope1-slope2

	if full_output:
		return diff_avg,diff_std,delta_l,delta_r,delta_bw,delta_cen,delta_sn,delta_sl
"""

def dct(xs,mode='s',modulus=False,full=False):
	"""Get the fourier transform of an interferogram using np.fft.fft() and np.fft.fftfreq()

	Parameters:
	xs (array) - 
	mode (str) - the type of interferogram data, 's' for single-sided and 'd' for double-sided,
				 default 's'
	modulus (bool) - if True, take the absolute value of the fourier transform data, default False
	full (bool) - if True and modulus == False, take the real+imaginary part of the fourier
				  transform data, otherwise take the real part only, default False

	Return:
	freqs, ys - the frequency and signal data from the fourier transform
	"""
	print '\nxs data given to dct function:'
	print xs
	xs_dict = {}
	xs_dict['interferogram_data'] = xs
	with open('/Users/lindsay/Desktop/xs_data.pkl', 'w') as f:
		pickle.dump(xs_dict,f)

	#for double-sided interferogram
	if mode == 'd':
		#get the absolute value of the FT, if desired
		if modulus:
			ys = np.absolute(np.fft.fft(xs))[:len(xs)/2]
		#get the full FT data, if desired
		elif full:
			ys = np.fft.fft(xs)[:len(xs)/2]
		#get the real part of the FT otherwise
		else:
			ys = np.fft.fft(xs).real[:len(xs)/2]
		#get the frequencies of the FT
		freqs = np.fft.fftfreq(len(xs))[:len(xs)/2]
		print '\nFrequencies from dct function:'
		print freqs
	
	#for single-sided interferogram
	elif mode == 's':
		#copy the single-sided data on both sides to get a full interferogram
		ext = np.concatenate((xs,xs[1:][::-1]))
		#get the absolute value of the FT, if desired
		if modulus:
			ys = np.absolute(np.fft.fft(ext))[:len(ext)/2]
		#get the full FT data, if desired
		elif full:
			ys = np.fft.fft(ext)[:len(ext)/2]
		#get the real part of the FT otherwise
		else:
			ys = np.fft.fft(ext).real[:len(ext)/2]
		#get the frequencies of the FT
		freqs = np.fft.fftfreq(len(ext))[:len(ext)/2]
    
    #exit if mode is invalid
	else:
		print ("Invalid mode \'" + mode + "\'' specified.  Exiting.")
    
	#return the frequencies and the FT values
	return freqs,ys

    
def fts_spectrum(xs_raw,name,window_name,plot_path=None,stepsize=738./125000.,mode='s',
				 modulus=False,porder=5,ezptest=True,include_tail=True):#,band='150GHz'): #highres:844, highran:591
	"""
	
	Parameters:
	xs_raw (array) - interferogram signals for mirror positions spaced 'stepsize' apart
	name (str) - bolometer name, or other identifier
	window_name () - 
	plot_path (str) - path for where to save the plots, default None
	stepsize (float) - inches per step in interferogram data, default 738./125000
					   * CHANGED FROM FRED WHO HAD DEFAULT 591./100000
	mode (str) - the type of interferogram data, 's' for single-sided and 'd' for double-sided,
			  default 's'
	modulus (bool) - True if want the absolute value of the fft returned, False for the real part,
				 default False
	porder (int) - degree to which to do the polynomial fitting of the interferogram data, default 5
	ezptest (bool) - if True, do some verification plotting of the zero path difference point,
					 default True
	include_tail (bool) - True if including the asymmetric part of the interferogram, only matters
						  if mode == 'd', default True
						  ???????????? NEED TO CHECK THAT THIS IS THE CORRECT INTERPRETATION ???????
	band - IGNORING FOR NOW
	"""

	#plot and save the raw interferogram
	#pl.figure()
	if plot_path is not None:
   		#pl.plot(xs_raw,marker='o')
		#pl.title('FTS Interferogram: Raw')
    		#pl.savefig(os.path.join(plot_path,'plots/%s_points.png'%(name)))
		#pl.clf()
	
		pl.plot(xs_raw)
		pl.title('FTS Interferogram: Raw')
		pl.xlabel('Linear Stage Step')
		pl.ylabel('Interferogram Raw Signal')
		pl.savefig(os.path.join(plot_path,'FTSAnalysis_Plots/%s.png'%(name)))
		pl.clf()

		#pl.plot(detrend(xs_raw,poly=porder))
		#pl.title('FTS Interferogram: %s th-order Polynomial Subtracted'%(porder))
                #pl.savefig(os.path.join(plot_path,'plots/%s_detrend.png'%(name)))
		#pl.clf()

	#find the index of the maximum signal value of the interferogram
	i_max = np.argmax(xs_raw)
	if i_max < 200:
		print 'Maximum index too low.  Recalculating using only higher indices.'
		i_max = 200 + np.argmax(xs_raw[200:])
		print ('Name: '+name + '\tMax Value: '+str(xs_raw[i_max])) #i_max+1#, xs_raw

	#make a dictionary to hold the data
	data = {}
	data['name'] = name

	#use detrend to remove a polynomial fit from the interferogram data and get the residual and fit
	xs_full, pfit = detrend(xs_raw,poly=porder,params=True)

	#estimate the index of zero path difference using the fitted interferogram data and split data
	ezp = est_zeropath(xs_full) #index corresponding to zero path differnce
	left = xs_full[:ezp]
	right = xs_full[ezp:2*ezp]
	tail = xs_full[2*ezp:] #any remaining data

	#if desired, plot the data to verify the zero path difference point
	if (plot_path is not None) and (ezptest):
		pl.plot(xs_full)
		pl.xlim([np.round(ezp*3./4.),np.round(ezp*5./4.)])
		pl.xlabel('Linear Stage Step')
		pl.ylabel('Interferogram Signal (polyfit subtracted)')
		pl.axvline(ezp,color='r',linestyle='--')
		pl.savefig(os.path.join(plot_path,'FTSAnalysis_Plots/ezp_%s_0.png'%(name)))
		pl.clf()    

	#assume the zero path estimate from above is correct and get the signal to noise
	if modulus:
		left_after = left
		right_after = right
		tail_after = tail
		#get the signal to noise using the located zero path difference point
		sn = calc_sn(xs_full)

	#get a better estimate for the zero path difference point
	#split the data and get the signal to noise using the new estimate
	else:
		#combine the right and left sides of the interferogram so the white light fringe is at the
		#beginning and then it goes through the left and right sides individually
		main_lobe = np.concatenate((right,left))
		#find the necessary shift from the zero point estimate
		sample_shift = zeropath_search(main_lobe)
		print "shift: ",sample_shift
		#get the filter using the calculated shift as the delay and convolve it with the
		#interferogram data
		filter = shift_filter(sample_shift)
		shifted_full = np.convolve(xs_full,filter,mode='valid')
		#find the new (better estimated) index of the zero path difference using the shifted data
		ezp = est_zeropath(shifted_full)

		#split the data using the new zero path point
		left_after = shifted_full[:ezp]
		if 2*ezp > len(shifted_full):
			right_after = np.zeros(ezp)
			end = len(shifted_full[ezp:])
			right_after[:end] = shifted_full[ezp:]
		else:
			right_after = shifted_full[ezp:2*ezp]
		tail_after = shifted_full[2*ezp:]


		#check the data after finding the zero path again
		xs_check = np.concatenate((right_after,left_after))
		#if abs(min(xs_check)) > max(xs_check):
		#	print 'Min larger than max'
		#	ok = False
		#else:

		try:
			ok = zeropath_check(xs_check,stepsize)
		except:
			ok = False


		#if desired, plot the data to verify the zero path difference point, again
		if (plot_path is not None) and (ezptest):
			pl.plot(np.concatenate((left_after,right_after)))
			pl.xlim([np.round(ezp*3./4.),np.round(ezp*5./4.)])
			pl.xlabel('Linear Stage Step')
			pl.ylabel('Interferogram Signal (polyfit subtracted, shifted)')
			pl.axvline(ezp,color='r',linestyle='--')
			pl.savefig(os.path.join(plot_path,'FTSAnalysis_Plots/ezp_%s_1.png'%(name)))
			pl.show()
			pl.clf()

			#pl.plot(np.concatenate((right_after[1:][::-1],right_after)))
			#pl.xlim([np.round(ezp*3./4.),np.round(ezp*5./4.)])
			#pl.axvline(ezp-1,color='r',linestyle='--')
			#pl.savefig(os.path.join(plot_path,'plots/ezp_%s_2.png'%(name)))
			#pl.clf()

		#if the zeropath_check fails, flip the interferogram and repeat the process to find the
		#zero path point
		if not ok:
			print 'May have missed WLF. Flipping interferogram...'
			#peaks = heapq.nlargest(2,xs_full)
			#i1 = np.where(xs_full==peaks[0])[0]
			#i2 = np.where(xs_full==peaks[1])[0]
			#ezp = int(np.round((i1+i2)/2.))
			#left_after = xs_full[:ezp]
                	#right_after = xs_full[ezp:2*ezp]
                	#tail_after = xs_full[2*ezp:]

			#filter = shift_filter(sample_shift/2.)
			#shifted_full = np.convolve(xs_full,filter,mode='valid')
			#peaks = heapq.nlargest(2,shifted_full)
                        #i1 = np.where(shifted_full==peaks[0])[0]
                        #i2 = np.where(shifted_full==peaks[1])[0]
                        #ezp = int(np.round((i1+i2)/2.))
                	#left_after = shifted_full[:ezp]
                	#right_after = shifted_full[ezp:2*ezp]
                	#tail_after = shifted_full[2*ezp:]

			xs_full *= -1.
			ezp = est_zeropath(xs_full)
			left = xs_full[:ezp]
			right = xs_full[ezp:2*ezp]
			tail = xs_full[2*ezp:]

			main_lobe = np.concatenate((right,left))
			sample_shift = zeropath_search(main_lobe)
			print "Re-shift: ",sample_shift
			filter = shift_filter(sample_shift)
			shifted_full = np.convolve(xs_full,filter,mode='valid')
			ezp = est_zeropath(shifted_full)
			#print ezp

			left_after = shifted_full[:ezp]
			right_after = shifted_full[ezp:2*ezp]
			tail_after = shifted_full[2*ezp:]

			if (plot_path is not None) and (ezptest):
				pl.plot(np.concatenate((left_after,right_after)))
				pl.xlim([np.round(ezp*3./4.),np.round(ezp*5./4.)])
				pl.axvline(ezp-1,color='r',linestyle='--')
				pl.savefig(os.path.join(plot_path,'FTSAnalysis_Plots/ezp_%s_3.png'%(name)))
				pl.clf()


		#get the signal to noise using the new zero path difference point
		sn = calc_sn(shifted_full,ezp)


		#???????????????? WHY IS THIS HERE ?????????????????????????????????????????????????????????
		def show_sample_delay():
			pl.subplot(211)
			pl.plot(right)
			pl.plot(arange(1,len(right)+1),left[::-1])
			pl.show()
			exit()
    
	#for a double-sided interferogram
	if mode == 'd':
		#combine the two sides of the interferogram (right then left side), with tail if desired
		if include_tail:
			right_after = np.concatenate((right_after,tail_after))
			left_after = np.concatenate((tail_after[::-1],left_after))
		xs = np.concatenate((right_after,left_after))
		#conbine the two sides of the interferogram in step order for plotting
		xs_save = np.concatenate((left_after,right_after))

		#plot the final interferogram if desired
		if plot_path is not None: 
			pl.plot(xs_save)
			pl.title('FTS Interferogram: %s th-order Polynomial Subtracted preFT'%(porder))
			pl.xlabel('Linear Stage Step')
			pl.ylabel('Interferogram Signal (polyfit subtracted, shifted)')
			pl.savefig(os.path.join(plot_path,'FTSAnalysis_Plots/%s_preFT.png'%(name)))
			pl.clf()

		#if sn is NaN and there is an odd number of data points, delete the last one
		#????????????????????????????? CONFUSED ABOUT THE ISNAN CONDITION ??????????????????????????
		if np.isnan(sn) and (ezp>len(xs)/2):
			print 'Found NaN interferogram with odd number of elements.'
			xs = xs[:len(xs)/2*2]
			print ("Interferogram data array length = " + str(len(xs)))

		#multiply the interferogram data by the hamming window function unless window == 'rectangle'
		#?????????????????????????? CONFUSED WHY ONLY RECTANGLE IS SPECIFIED ???????????????????????
		if window_name != 'rectangle':
			window_l = np.hamming(len(xs))[:len(xs)/2]
			window_r = window_l[::-1]
			window = np.concatenate((window_r,window_l))
			print 'Window Length: ' + str(len(window))
			old_xs = xs
			xs*= window
	
	#for a single-sided interferogram
	elif mode == 's':
		#combine the single-sided interferogram data
		#right_avg = (right_after+left_after[::-1])/2.
		xs = np.concatenate((right_after,tail_after))
		xs_save = np.copy(xs)

		#multiply the interferogram data by the hamming window function unless window == 'rectangle'
		#?????????????????????????? CONFUSED WHY ONLY RECTANGLE IS SPECIFIED ???????????????????????
		if window_name != 'rectangle':
			window = np.hamming(2*len(xs))[:len(xs)][::-1]
			xs*= window
    
	#exit if mode is invalid
	else:
		print ("Invalid mode \'" + mode + "\'' specified.  Exiting.")



	#take the fourier transform and get the data
	#fs = freqs, ys = FT signal
	fs,ys = dct(xs,mode=mode,modulus=modulus)

	#normalize the FT signal values against the maximum
	#ys /= max(ys[10:])
	
	#convert the frequencies from inverse steps to GHz (1/ns)
	c = 11.8 #speed of light in inches/ns
	inches_per_step = stepsize
	fs /= 2*inches_per_step/c


	"""COMMENTING THIS OUT NOW BUT MAY WANT TO PUT BACK LATER ?????????????????

	if (plot_path is not None) and (not modulus):
		fs_temp,ys_temp = dct(xs,mode=mode,modulus=False,full=True)
		ys_temp /= max(ys_temp.real)
		fs_temp /= 2*inches_per_step/c
		plot_phase(fs_temp,ys_temp,name,plot_path)
	elif (plot_path is not None) and (modulus):
		fs_temp,ys_temp = dct(xs,mode=mode,modulus=False,full=True)
		ys_temp /= max(ys)
		fs_temp /= 2*inches_per_step/c
		plot_phase(fs_temp,ys_temp,name+'_nocorr',plot_path)
	"""

	
	"""COMMENTING OUT THE BAND STUFF FOR NOW ???????????????????????????????????????????????????????

	#atmData = getAtmData('/Users/fmatsuda/FTSFinalTest_Berkeley/FTS/ATM_PWV1p0mm_90Deg.out')
	#atmData = getAtmData('/Users/fmatsuda/FTS_Chile/FTS_berkrepo/FTS/AtmosphericModel/Chaj_zadegrees30_pwvmicrons1000.out')
	atmData = getAtmData('/home/cmb/fmatsuda/FTS_run38/AtmosphericModel/Chaj_zadegrees30_pwvmicrons1000.out')

	#try:   
	if band == '150GHz':
 		integ_bandwidth,center,temp,l,r = est_band2(fs,ys,atmData,fmin_l=100.,fmin_r=150)
	elif band == '90GHz':
 		integ_bandwidth,center,temp,l,r = est_band2(fs,ys,atmData,fmin_l=30.,fmin_r=90)
	else:
 		integ_bandwidth,center,temp,l,r = est_band2(fs,ys,atmData)
	#except:
	#	integ_bandwidth = 0.
	#	center = 0.
	#	temp = 0.
	#	l = 0
	#	r = 0
	"""

	#pl.close()   
 
 	#fill in the dictionary to return
 	#COMMENTING OUT THE BAND STUFF AGAIN ???????????????????????????????????????????????????????????
	data['interferogram'] = xs_raw
	data['interferogram_detrend'] = xs_save
	data['throughput'] = ys
	data['frequencies'] = fs
	#data['integ_bandwidth'] = integ_bandwidth
	#data['eff_frequency'] = center
	#data['left_edge'] = l
	#data['right_edge'] = r
	data['window'] = window_name
	data['polynomial'] = porder
	data['sig2noise'] = sn    
	data['slope'] = pfit[1]
	data['pfit'] = pfit
	#data['eff_temperature'] = temp
	data['gain'] = np.NaN

	return data


def make_plots_dir(output):
	"""Make a directory to hold the generated output plots if one does not yet exist.
	The new directory will be called 'output_path/plots'.

	Parameters:
	output_path (str) - the desired path to the output plots
	"""
	if not os.path.exists(os.path.join(output,'FTSAnalysis_Plots')):
		print 'Creating plots output directory...'
		cmd = 'mkdir '+os.path.join(output,'FTSAnalysis_Plots')
		ret = subprocess.call(cmd,shell=True)
		if ret != 0:
			print 'Failed creating plots directory in specified output directory. Exiting...'
			sys.exit(1)
		else:
			print ('Plots output directory ' + os.path.join(output,'FTSAnalysis_Plots') +
				   ' successfully created.')
	else:
		print 'Plots output directory already exists.  No new directory created.'


def plot_spectrum(data, style='-'):
	"""Plot the spectrum using the info from the fts_spectrum() function

	Parameters:
	data (dict) - dictionary with spectral data returned by fts_spectrum()
	style (str) - marker style for plotting, default '-'
	"""

	#load all the data from the dictionary
	window = data['window']
	ys = data['throughput']
	fs = data['frequencies']
	#bw = data['integ_bandwidth']
	#l = data['left_edge']
	#r = data['right_edge']
	name = data['name']
	poly = data['polynomial']

	#plot the data, ignoring the band info for now
	#pl.plot(fs,ys,style,label=name+' '+window+' %.1f GHz '%bw+'%sth-order'%(poly))
	pl.plot(fs,ys,style,linewidth=1,label=name+' '+window+' %sth-order'%(poly))
    

def plot_conf():
	"""Configure the plot"""
	pl.title('FTS Spectrum')
	pl.xlabel('Frequency (GHz)')
	leg = pl.legend(loc=3)
	for t in leg.get_texts():
		t.set_fontsize('small')
	pl.grid()

""" DONT NEED THIS NOW =============================================================================
def plot_stats(data):
	xs = data['throughput']
    fs = data['frequencies']
	bw = data['integ_bandwidth']
	c = data['eff_frequency']
    l = data['left_edge']
    r = data['right_edge']
	t = data['eff_temperature']

	pl.axvline(fs[l],color='k',linestyle='--')
   	pl.axvline(fs[r],color='k',linestyle='--')
	#pl.axvline(c,color='k',linestyle='--')

	pl.figtext(0.7,0.75,'BW = %.1f GHz'%(bw))
	pl.figtext(0.7,0.7,'Center = %.1f GHz'%(c))
	pl.figtext(0.7,0.65,'Temp = %.3f K'%(t))
	pl.figtext(0.7,0.6,'Tot Temp = %.3f'%(t*bw))
"""

""" DONT NEED THIS NOW =============================================================================
def analyze_all_chans_berk(hwmap_path,pklfn,output_path,amode):
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION
	make_plots_dir(output_path)

	os.environ['HWMAP_PATH'] = hwmap_path
    hwmap = os.environ['HWMAP_PATH']
    devman = DeviceManager(hwmap)
	
	print 'Loading data pickle file...'
	data = load_dict(pklfn)
	chans = []
    for k in data.keys():
        try:
            chans.append(str(sqch_of_wtl(devman,k)))
        except:
            pass
	chans.sort()

	print 'Start calculating and plotting spectra for all channels...'
	for chan in chans:
    	print chan
    	wtl_name = wtl_of_sqch(devman,chan)
    	xs_raw = get_array(data[wtl_name],'Mag')
    	energy = np.std(detrend(xs_raw))
    	if energy < 2: continue

    	spectral_data = fts_spectrum(xs_raw,chan,'triangle',plot_path=output_path,mode=amode,ezptest=True)

		pl.figure()
    	plot_spectrum(spectral_data,'o-')
    	plot_conf()
    	pl.savefig(os.path.join(output_path,'plots/fts-'+chan+'_points.png'))
		pl.clf()		

		plot_spectrum(spectral_data)
		#atmData = getAtmData('/home/fmatsuda/FTSFinalTest_Berkeley/FTS/ATM_PWV1p0mm_90Deg.out')
		#atmData = getAtmData('/Users/fmatsuda/FTS_Chile/FTS_berkrepo/FTS/AtmosphericModel/Chaj_zadegrees30_pwvmicrons1000.out')
		atmData = getAtmData('/global/homes/f/fmatsuda/FTS/AtmosphericModel/Chaj_zadegrees30_pwvmicrons1000.out')
    	pl.plot(atmData['Frequencies'],atmData['Transmission'],linestyle = '--',color = 'r')
		pl.xlim([0,250])
		pl.ylim([-0.25,1.0])
		plot_stats(spectral_data)
		plot_conf()
		pl.savefig(os.path.join(output_path,'plots/fts-'+chan+'.png'))
    	pl.close()

		spectral_data_m = fts_spectrum(xs_raw,chan,'triangle',plot_path=output_path,mode=amode,modulus=True)

        pl.figure()
        plot_spectrum(spectral_data_m,'o-')
        plot_conf()
        pl.savefig(os.path.join(output_path,'plots/fts-'+chan+'_points_modulus.png'))
        pl.clf()

        plot_spectrum(spectral_data_m)
        pl.xlim([0,250])
        pl.ylim([-0.25,1.0])
        plot_conf()
        pl.savefig(os.path.join(output_path,'plots/fts-'+chan+'_modulus.png'))
        pl.close()

		pl.figure()
		plot_spectrum(spectral_data)
		plot_spectrum(spectral_data_m)
		plot_conf()
		pl.savefig(os.path.join(output_path,'plots/ftscomp-'+chan+'.png'))
        pl.close()

	pairs = [('K01Sq1Ch2','K01Sq1Ch6'),('K01Sq5Ch2','K01Sq5Ch6')]
	plot_pairs_berk(devman,data,pairs,output_path,amode)
"""

""" DONT NEED THIS NOW =============================================================================
def analyze_all_chans_apex(pklfn,output_path,amode):
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION
    make_plots_dir(output_path)

    print 'Loading data pickle file...'
    data = load_dict(pklfn)
    chans = ['Mb97Mux1Ch1','Mb97Mux1Ch4']#,'Mb99Mux1Ch3']

    print 'Start calculating and plotting spectra for all channels...'
    for chan in chans:
        print chan
        wtl_name = chan
        xs_raw = get_array(data[wtl_name],'Mag')
        energy = np.std(detrend(xs_raw))
        if energy < 2: continue

        spectral_data = fts_spectrum(xs_raw,chan,'triangle',plot_path=output_path,mode=amode,ezptest=True)

        pl.figure()
        plot_spectrum(spectral_data,'o-')
        plot_conf()
        pl.savefig(os.path.join(output_path,'plots/fts-'+chan+'_points.png'))
        pl.clf()

        plot_spectrum(spectral_data)
        pl.xlim([0,250])
        pl.ylim([-0.25,1.0])
        plot_stats(spectral_data)
        plot_conf()
        pl.savefig(os.path.join(output_path,'plots/fts-'+chan+'.png'))
        pl.close()

        spectral_data_m = fts_spectrum(xs_raw,chan,'triangle',plot_path=output_path,mode=amode,modulus=True)

        pl.figure()
        plot_spectrum(spectral_data_m,'o-')
        plot_conf()
        pl.savefig(os.path.join(output_path,'plots/fts-'+chan+'_points_modulus.png'))
        pl.clf()

        plot_spectrum(spectral_data_m)
        pl.xlim([0,250])
        pl.ylim([-0.25,1.0])
        plot_conf()
        pl.savefig(os.path.join(output_path,'plots/fts-'+chan+'_modulus.png'))
        pl.close()

		pl.figure()
        plot_spectrum(spectral_data)
        plot_spectrum(spectral_data_m)
        plot_conf()
        pl.savefig(os.path.join(output_path,'plots/ftscomp-'+chan+'.png'))
        pl.close()

		if chan == 'Mb97Mux1Ch1':
			fs = spectral_data['frequencies']
			xs = spectral_data['throughput']
			xs_low = xs[fs<50.]
			xs_low = detrend(xs_low,poly=1)
			xs = np.concatenate((xs_low,xs[fs>=50.]))

			g = gauss_kern(4)
			xs = np.convolve(xs,g,mode='same')

			xs /= max(xs)
			spectral_data['frequencies'] = fs
			spectral_data['throughput'] = xs
			spectral_data['name'] = '150 GHz'
			spectral_data1 = spectral_data
		else:
			fs = spectral_data['frequencies']
            xs = spectral_data['throughput']

			g = gauss_kern(4)
            xs = np.convolve(xs,g,mode='same')

			xs /= max(xs)
			spectral_data['frequencies'] = fs
			spectral_data['throughput'] = xs
			spectral_data['name'] = '90 GHz'
			spectral_data2 = spectral_data

	pl.figure()
	plot_spectrum(spectral_data1,'-')
	plot_spectrum(spectral_data2,'-')
	#atmData = getAtmData('/home/fmatsuda/FTSFinalTest_Berkeley/FTS/ATM_PWV1p0mm_90Deg.out')
	#atmData = getAtmData('/Users/fmatsuda/FTS_Chile/FTS_berkrepo/FTS/AtmosphericModel/Chaj_zadegrees30_pwvmicrons1000.out')
	atmData = getAtmData('/global/homes/f/fmatsuda/FTS/AtmosphericModel/Chaj_zadegrees30_pwvmicrons1000.out')
    pl.plot(atmData['Frequencies'],atmData['Transmission'],linestyle = '--',color = 'r')
	pl.xlim([0,250])
	pl.ylim([-0.2,1.1])
	plot_conf()
	pl.savefig(os.path.join(output_path,'plots/fts-90_150.png'))
	pl.close()
"""

""" DONT NEED THIS NOW =============================================================================
def analyze_all_chans_chile(data,hwmapfn,output_path,amode,select_wafer=None,gain_cal=False):
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION
	make_plots_dir(output_path) 
	plot_blank(output_path)

	wafers = ['1','2','3','4','5','6','7']
    wafernames = ['10.2','10.4','10.3','10.5','9.4','10.1','8.2.0']
	xmlmap = parsemap.build_index_maps(hwmapfn)
	wafernums = np.array(xmlmap.wafer)

	QLout = {}
	dataout = {}
	dataout['Interferograms'] = data

	chans = range(1512)
	if select_wafer is not None:
		wchs = np.where(wafernums==str(select_wafer))[0]
   		chans = np.intersect1d(chans,wchs)
		#nochans = np.setdiff1d(chans,wchs)

		#pair_list0 = libpolmap.pairlist(nochans,xmlmap.boloid)
		#for pair in pair_list0:
		#	boloid0 = xmlmap.boloid[pair[0]]
		#	dlist0 = {'Spectrum':None,'Diff_Band':None}
		#	list0 = {'Interferogram':None,'Spectrum':None,'Phase':None,'SNR':None,'Eff_Freq':None,'BandWidth':None,'Integ_BandWidth':None}
		#	QLout[boloid0] = {'Dual': dlist0,'Top': list0,'Bottom': list0}

	#atmData = getAtmData('/Users/fmatsuda/FTSFinalTest_Berkeley/FTS/ATM_PWV1p0mm_90Deg.out')
	#atmData = getAtmData('/Users/fmatsuda/FTS_Chile/FTS_berkrepo/FTS/AtmosphericModel/Chaj_zadegrees30_pwvmicrons1000.out')
	atmData = getAtmData('/global/homes/f/fmatsuda/FTS/AtmosphericModel/Chaj_zadegrees30_pwvmicrons1000.out')
	pair_list = libpolmap.pairlist(chans,xmlmap.boloid)
	
	for pair in pair_list:
		t_ok = True
		b_ok = True

		tchan = pair[0]
		bchan = pair[1]
		boloid_t = xmlmap.boloid[tchan]
		boloid_b = xmlmap.boloid[bchan]
	
		energy_t = np.std(detrend(data[tchan]))
		energy_b = np.std(detrend(data[bchan]))
    	if energy_t < 0.25: 
			t_ok = False
			print '%s energy less than 0.25'%(boloid_t)
		if energy_b < 0.25: 
			b_ok = False 
			print '%s energy less than 0.25'%(boloid_b)

		#if len(data[tchan]) < 100:
		#	t_ok = False
		#	print '%s length of data less than 100'%(boloid_t)
		#if len(data[bchan]) < 100:
		#	b_ok = False
		#	print '%s length of data less than 100'%(boloid_b)

		max_t = np.max(data[tchan])
		max_b = np.max(data[bchan])
		if np.isnan(max_t): 
			t_ok = False
			print '%s NaN max with %s'%(boloid_t,len(np.isnan(data[tchan])))
		if np.isnan(max_b): 
			b_ok = False
			print '%s NaN max with %s'%(boloid_b,len(np.isnan(data[tchan])))
	
		#pl.figure()

		if t_ok:
			spectral_data_t = fts_spectrum(data[tchan],boloid_t,'triangle',plot_path=output_path,mode=amode)

			plot_spectrum(spectral_data_t,'o-')
			plot_conf()
			pl.savefig(os.path.join(output_path,'plots/fts-%s_points.png'%(boloid_t)))
			pl.clf()

			plot_spectrum(spectral_data_t)
        		pl.plot(atmData['Frequencies'],atmData['Transmission'],linestyle = '--',color = 'r')
			#pl.xlim([0,250])
			pl.ylim([-0.25,1.0])
			plot_stats(spectral_data_t)
			plot_conf()
			pl.savefig(os.path.join(output_path,'plots/fts-%s.png'%(boloid_t)))
			pl.clf()

			bw_t = spectral_data_t['frequencies'][spectral_data_t['left_edge']] - spectral_data_t['frequencies'][spectral_data_t['right_edge']]

		if b_ok:
			spectral_data_b = fts_spectrum(data[bchan],boloid_b,'triangle',plot_path=output_path,mode=amode)
		       
			plot_spectrum(spectral_data_b,'o-')
			plot_conf()
			pl.savefig(os.path.join(output_path,'plots/fts-%s_points.png'%(boloid_b)))
			pl.clf()

			plot_spectrum(spectral_data_b)
        		pl.plot(atmData['Frequencies'],atmData['Transmission'],linestyle = '--',color = 'r')
			#pl.xlim([0,250])
			pl.ylim([-0.25,1.0])
			plot_stats(spectral_data_b)
			plot_conf()
			pl.savefig(os.path.join(output_path,'plots/fts-%s.png'%(boloid_b)))
			pl.clf()

			bw_b = spectral_data_b['frequencies'][spectral_data_b['left_edge']] - spectral_data_b['frequencies'][spectral_data_b['right_edge']]

		if t_ok and b_ok:
			plot_spectrum(spectral_data_t)
			plot_spectrum(spectral_data_b)
			#pl.xlim([0,300])
			pl.ylim([-0.25,1.0])
			plot_conf()
			pl.savefig(os.path.join(output_path,'plots/fts-%s.png'%(boloid_t[:-1])))
			pl.clf()

			diff_band(spectral_data_t,spectral_data_b,output_path)
			if gain_cal:
				print 'Calculating relative gain parameter between pixel pair'
				g = calc_gain(spectral_data_t,spectral_data_b,atmData)
				print 'Gain = ',g
				spectral_data_t['gain'] = 1.0

				integ_bandwidth,center,temp,l,r = est_band2(spectral_data_b['frequencies'],g*spectral_data_b['throughput'],atmData)
				spectral_data_b['throughput'] *= g
				spectral_data_b['integ_bandwidth'] = integ_bandwidth
				spectral_data_b['eff_frequency'] = center
				spectral_data_b['left_edge'] = l
				spectral_data_b['right_edge'] = r
				spectral_data_b['eff_temperature'] = temp
				spectral_data_b['gain'] = g

		#pl.close()

		if t_ok and b_ok:
			dlist = {'Spectrum':'fts-%s.png'%(boloid_t[:-1]),'Diff_Band':'ftsdiff-%s.png'%(boloid_t[:-1])}
		else:
			dlist = {'Spectrum':'blank.png','Diff_Band':'blank.png'}

		if t_ok:
        		tlist = {'Interferogram':'%s.png'%(boloid_t),'Spectrum':'fts-%s.png'%(boloid_t),'Phase':'phase-%s.png'%(boloid_t),'SNR':spectral_data_t['sig2noise'],'Eff_Freq':spectral_data_t['eff_frequency'],'BandWidth':bw_t,'Integ_BandWidth':spectral_data_t['integ_bandwidth'],'Interfero_Slope':spectral_data_t['slope'],'Eff_Temp':spectral_data_t['eff_temperature']}
		else:
			tlist = {'Interferogram':'blank.png','Spectrum':'blank.png','Phase':'blank.png','SNR':0.,'Eff_Freq':0.,'BandWidth':0.,'Integ_BandWidth':0.,'Interfero_Slope':-1.,'Eff_temp':0.}

		if b_ok:
        		blist = {'Interferogram':'%s.png'%(boloid_b),'Spectrum':'fts-%s.png'%(boloid_b),'Phase':'phase-%s.png'%(boloid_b),'SNR':spectral_data_b['sig2noise'],'Eff_Freq':spectral_data_b['eff_frequency'],'BandWidth':bw_b,'Integ_BandWidth':spectral_data_b['integ_bandwidth'],'Interfero_Slope':spectral_data_b['slope'],'Eff_Temp':spectral_data_b['eff_temperature']}
		else:
			blist = {'Interferogram':'blank.png','Spectrum':'blank.png','Phase':'blank.png','SNR':0.,'Eff_Freq':0.,'BandWidth':0.,'Integ_BandWidth':0.,'Interfero_Slope':-1.,'Eff_temp':0.}

        	QLout[boloid_t[:-1]] = {'Dual': dlist,'Top': tlist,'Bottom': blist}
		dataout[boloid_t[:-1]] = {}
		if t_ok:
			dataout[boloid_t[:-1]]['Top'] = spectral_data_t
		if b_ok:
			dataout[boloid_t[:-1]]['Bottom'] = spectral_data_b

	return QLout, dataout
"""

#def show_all():
#    	all_xs /= max(all_xs)
#    	l,r,bw = est_band(all_fs,all_xs)
#    	pl.plot(all_fs,all_xs, label='Average of all pixels bw=%.1f GHz'%bw)
#    	pl.plot(all_fs,all_xs>0.35)
#    	pl.axvline(all_fs[l])
#    	pl.axvline(all_fs[r])


def plot_blank(output_path):
	"""Make a plot of an array of zeros and save it to the output plots directory"""
	#pl.figure()
	pl.plot(np.zeros(100))
	pl.savefig(os.path.join(output_path,'FTSAnalysis_Plots/blank.png'))
	pl.clf()
	#pl.close()

""" DONT NEED THIS NOW =============================================================================
def plot_pairs_berk(devman,data,pairs,output_path,amode):
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION
	for pair in pairs:
		wtl_name0 = wtl_of_sqch(devman,pair[0])
        xs_raw0 = get_array(data[wtl_name0],'Mag')
		wtl_name1 = wtl_of_sqch(devman,pair[1])
        xs_raw1 = get_array(data[wtl_name1],'Mag')

		spectral_data0 = fts_spectrum(xs_raw0,pair[0],'triangle',mode=amode)
      	spectral_data1 = fts_spectrum(xs_raw1,pair[1],'triangle',mode=amode)

        pl.figure()
        plot_spectrum(spectral_data0)
        plot_spectrum(spectral_data1)
		pl.xlim([0,300])
        pl.ylim([-0.25,1.0])
        plot_conf()
        pl.savefig(os.path.join(output_path,'plots/ftspair-'+pair[0]+'-'+pair[1]+'.png'))
        pl.close()
	
		diff_band(spectral_data0,spectral_data1,output_path,pix_mode=False)
"""

""" DONT NEED THIS NOW =============================================================================
def gauss_kern(sigma):
	isize = int(sigma)
	x = np.mgrid[-2*isize:2*isize+1]
	g = np.exp(-(np.float_(x**2)/(2.*sigma**2)))
	return g/max(g)
"""

""" DONT NEED THIS NOW =============================================================================
def getAtmData(atmFilename):
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION
	atmFile = open(atmFilename,'r')
	atmData = dict()
	atmData['Frequencies'] = list()
	atmData['Transmission'] = list()
	atmData['Temperature'] = list()
	for line in atmFile:
    	lineSplit = line.split()
    	atmData['Frequencies'].append(float(lineSplit[0]))
    	#atmData['Transmission'].append(float(lineSplit[3]))
		atmData['Transmission'].append(float(lineSplit[1]))
		atmData['Temperature'].append(float(lineSplit[2]))
	atmFile.close()
	return atmData
"""

""" DONT NEED THIS NOW =============================================================================
def GenQuadPair(ref,harm=1):#,phase_corr=np.NaN):
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION
    '''
    First count rising edges in ref, calculate freq from spacing
    Expects either a square wave or pulse input, triggers on rising edge
    harm is harmonic of input reference signal to fit quad pair to
    '''
    maxref=max(ref)   # Max level in signal
    minref=min(ref)   # Min level in signal
    highlevel=minref + 0.5*(maxref-minref)
    lowlevel=minref + 0.25*(maxref-minref)
    edgelist=[]
    curstate='high'
    for i in range(len(ref)):
        if (curstate=='high'):
            if ref[i] < lowlevel:
                curstate='low'
            else:
                if ref[i] > highlevel:
                    curstate='high'
                    edgelist.append(i)

        if len(edgelist)>2:
            period=(np.float(edgelist[-1]-edgelist[0]))/(len(edgelist)-1.)
            freq=1./period*harm
            #print freq, period
            phase=np.average(np.array(edgelist)%period)
		#if not np.isnan(phase_corr):
		#	X=np.sin(freq*2.*np.pi*(pl.frange(len(ref)-1)-phase-phase_corr))
                #	Y=-np.cos(freq*2.*np.pi*(pl.frange(len(ref)-1)-phase-phase_corr))
		#else:
            X=np.sin(freq*2.*np.pi*(pl.frange(len(ref)-1)-phase))
            Y=-np.cos(freq*2.*np.pi*(pl.frange(len(ref)-1)-phase))
        else:
            X = 0*np.array(range(len(ref)) )
            Y = X
    return {'X':X,'Y':Y}
"""

""" DONT NEED THIS NOW =============================================================================
def GenQuadPair2(ref,harm=1,extra=0,deltat=0.72879752):
	sample_rate = 25e6/2.**17.
	nt = len(ref)
	ts = np.arange(nt) / sample_rate

	#def model(ts,ampx,ampy,f):
	#	phi = 2.0*np.pi*f*harm*ts
	#	return ampx*np.cos(phi) + ampy*np.sin(phi)

	def model(ts,ampx,ampy,f,theta):
		phi = 2.0*np.pi*f*harm*ts + theta
		return ampx*np.cos(phi) + ampy*np.sin(phi)

	#popt,pcov = curve_fit(model,ts,ref,(0.,0.,4.0005))
	popt,pcov = curve_fit(model,ts,ref,(0.,0.,4.0005,0.))
	print popt
	
	if extra > 0:
		ts = np.arange(nt+extra) / sample_rate
	#print deltat
	ts += deltat
	ref = model(ts,*popt)
	#refx = model(ts,1.0,0.0,popt[2])
	#refy = model(ts,0.0,1.0,popt[2])
	refx = model(ts,1.0,0.0,popt[2],popt[3])
	refy = model(ts,0.0,1.0,popt[2],popt[3])

	return {'X':refx,'Y':refy}
"""

""" DONT NEED THIS NOW =============================================================================
def GenQuadPairValue(time,harm=1,extra=0,deltat=0.):
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION
    #sample_rate = 25e6/2.**17.
    #nt = len(time)
    #ts = np.arange(nt) / sample_rate
	ts = time

    def model(ts,ampx,ampy,f,theta):
        phi = 2.0*np.pi*f*harm*ts + theta
        return ampx*np.cos(phi) + ampy*np.sin(phi)

    refx = model(ts,1.0,0.0,4.0,0.0)
    refy = model(ts,0.0,1.0,4.0,0.0)

    return {'X':refx,'Y':refy}
"""

""" DONT NEED THIS NOW =============================================================================
def channel_noise(timestream,fmin,fmax,fs):
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION
    timestream = np.array(timestream)
    n = len(timestream)
	#print n
    fourier = np.fft.fft(timestream)/n**0.5        # normalize
    nmin = int(n*(fmin/fs))
    nmax = int(n*(fmax/fs))
    chanbw = fs/n

    # integrate noise energy in the power spectrum
    sum = 0.0
    total_bw = 0
    for i in np.arange(nmin,nmax):
        chan_energy = chanbw*abs(fourier[i])**2
        sum += 2*chan_energy if i!= 0 else chan_energy  # include negative frequency contribution for non-DC

        total_bw += 2 if i!= 0 else 1

        noise_energy = sum**0.5

	#print chanbw, total_bw
    if total_bw == 0: 
		total_bw = 1
		print '  Zero'

    noise_power = noise_energy/((chanbw**0.5)*total_bw)

    return noise_power
"""

""" DONT NEED THIS NOW =============================================================================
def Demod(data,refsig,harm=1,drop_data=None,phase_corr=False,phase_opt=None):#phase_opt=False):
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION 
    '''
    Data is a dictionary of test data, refname is the dict key for the ref channel
    harm is harmonic of reference data to fit quadrature pair to
    '''
    ref=GenQuadPair(refsig,harm=harm)

	if drop_data is not None:
		data = data[drop_data]
		ref['X'] = ref['X'][drop_data]
		ref['Y'] = ref['Y'][drop_data]

        x = 2.*np.array(data-np.mean(data))*np.array(ref['X'])
        y = 2.*np.array(data-np.mean(data))*np.array(ref['Y'])

        Xout = np.average(x)
        Yout = np.average(y)
        Phi = np.arctan2(Yout,Xout)

        #noise = channel_noise(data[item],1.0,20.0,381.47)
        noise = channel_noise(data,35,55,190.735)

	if phase_corr and (phase_opt is not None):
		print 'Both phase correction and optimization are selected. Defaulting to phase correction.'
		phase_opt=None

	if phase_corr:
		#print 'Apply phase correction'
		#z = Xout+Yout*1j
		#Mag = (z*np.exp(-Phi*1j)).real
		#Mag = (z*np.exp(-Phi*1j)).imag
		Mag = Xout*np.cos(Phi) + Yout*np.sin(Phi)	
		#print Mag, Phi
	#elif phase_opt:
	#	#print 'Apply phase optimization'
	#	def resid3(c):
	#		Xopt = x*np.cos(c) + y*np.sin(c)
	#		return -np.std(Xopt)
	#	p0 = 0.
        #	p = minimize(resid3,p0)
	#	#print p.x
	#	Mag = Xout*np.cos(p.x[0]) + Yout*np.sin(p.x[0])
	elif phase_opt is not None:
		Mag = Xout*np.cos(phase_opt) + Yout*np.sin(phase_opt)
		Phi = phase_opt 
	else:
    	Mag = (Xout**2 + Yout**2)**0.5

    out={'X':Xout,'Y':Yout,'Mag':Mag,'Phi':Phi*180/np.pi, 'Noise':noise}

    return out
"""

""" DONT NEED THIS NOW =============================================================================
def Demod2(data,refIQ,drop_data=None,phase_corr=False,phase_opt=np.NaN):
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION
	if drop_data is not None:
		data = data[drop_data]
		refIQ['X'] = refIQ['X'][drop_data]
		refIQ['Y'] = refIQ['Y'][drop_data]

        x = 2.*np.array(data-np.mean(data))*np.array(refIQ['X'])
        y = 2.*np.array(data-np.mean(data))*np.array(refIQ['Y'])

        Xout = np.average(x)
        Yout = np.average(y)
        Phi = np.arctan2(Yout,Xout)

        #noise = channel_noise(data[item],1.0,20.0,381.47)
        noise = channel_noise(data,35,55,190.735)

	if phase_corr and (not np.isnan(phase_opt)):
		print 'Both phase correction and optimization are selected. Defaulting to phase correction.'
		phase_opt=np.NaN

	if phase_corr:
		#print 'Apply phase correction'
		#z = Xout+Yout*1j
		#Mag = (z*np.exp(-Phi*1j)).real
		#Mag = (z*np.exp(-Phi*1j)).imag
		Mag = Xout*np.cos(Phi) + Yout*np.sin(Phi)	
		#print Mag, Phi
	#elif phase_opt:
	#	#print 'Apply phase optimization'
	#	def resid3(c):
	#		Xopt = x*np.cos(c) + y*np.sin(c)
	#		return -np.std(Xopt)
	#	p0 = 0.
        #	p = minimize(resid3,p0)
	#	#print p.x
	#	Mag = Xout*np.cos(p.x[0]) + Yout*np.sin(p.x[0])
	elif not np.isnan(phase_opt):
		Mag = Xout*np.cos(phase_opt) + Yout*np.sin(phase_opt)
		Phi = phase_opt 
	else:
    	Mag = (Xout**2 + Yout**2)**0.5

    out={'X':Xout,'Y':Yout,'Mag':Mag,'Phi':Phi*180/np.pi, 'Noise':noise}

    return out
"""

""" DONT NEED THIS NOW =============================================================================
def Demod_chan(data,refIQ,drop_data,phase_opt=np.NaN):
	#NEED TO CHEK THE INDENTATION IF ACTUALLY USING THIS FUNCTION
	data_mean = np.zeros((data.shape[0],1))
	for i in range(data.shape[0]):
		data_mean[i] = np.mean(data[i][drop_data[i]])
	#print data_mean.shape, data_mean
	x = 2.*np.array(data-data_mean)*np.array(refIQ['X'])
    y = 2.*np.array(data-data_mean)*np.array(refIQ['Y'])
	#print x.shape, x
	#print y.shape, y

	#pl.plot((data-data_mean)[653])
	#pl.plot(refIQ['X'][653])
	#pl.plot(refIQ['Y'][653])
	#pl.show()

	Xout = np.zeros(x.shape[0])
	Yout = np.zeros(y.shape[0])
	noise = np.zeros(x.shape[0])
	for j in range(x.shape[0]):
		if len(data[j][drop_data[j]]) < 100.:
			Xout[j] = 0.
			Yout[j] = 0.
			noise[j] = 0.
		else:
			Xout[j] = np.average(x[j][drop_data[j]])
			Yout[j] = np.average(y[j][drop_data[j]])
			noise[j] = channel_noise(data[j][drop_data[j]],35,55,190.735)

	if not np.isnan(phase_opt):
		#Mag = Xout*np.cos(phase_opt) + Yout*np.sin(phase_opt)
		Zout = Xout+1.j*Yout
		Z_rot = Zout*np.exp(-1.j*phase_opt)
    	Mag = Z_rot.real
    	iMag = Z_rot.imag
		Xout = Mag
		Yout = iMag
		Phi = np.average(np.arctan2(Yout,Xout))
	else:
    	Mag = (Xout**2 + Yout**2)**0.5
		Phi = np.average(np.arctan2(Yout,Xout))
		iMag = np.NaN

	out={'X':Xout,'Y':Yout,'Mag':Mag,'iMag':iMag,'Phi':Phi*180./np.pi,'Noise':noise}

    return out
"""

if __name__=="__main__":

	parser = argparse.ArgumentParser(description='Analyze FTS data and calculate spectra')
	parser.add_argument('-hm', dest='hwmap', action='store', help='Specify input hwmap directory')
	parser.add_argument('-i', dest='inputfn', action='store', help='Specify input data pickle file')
	parser.add_argument('-o', dest='outputdir', action='store', default='.', help='Specify output directory')
	parser.add_argument('-m', dest='mode', action='store', default='d', help='Specify if we want to do a single-sided analysis or a double-sided analysis')
	args = parser.parse_args()
    
	analyze_all_chans_berk(args.hwmap,args.inputfn,args.outputdir,args.mode)
	#analyze_all_chans_apex(args.inputfn,args.outputdir,args.mode)
    
     
    
