#!/usr/bin/env python

import os
import sys
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as pl
import cPickle as pickle
import argparse
#from scipy.optimize import leastsq, fmin_bfgs, curve_fit
from scipy.optimize import leastsq, fmin_bfgs, minimize, curve_fit
#import pygetdata
import time
import calendar
#import netCDF4

#sys.path.append('/home/polarbear/fmatsuda/FTS')
#sys.path.append('/global/homes/f/fmatsuda/FTS_KEK')
sys.path.append('/home/cmb/fmatsuda/FTS_run38')
sys.path.append('/home/cmb/fmatsuda/g3tests')
print 'Load libg3'
import libg3
print 'Load FTS Analysis Library'
import lib_FTSAnalysis as lFA
print 'Load SPT3G software'
from spt3g import core

#sys.path.append('/group/cmb/polarbear/usr/nishino/opt2/h5py2.7/lib/python2.7/site-packages')
import h5py


def read_refdata(refpath,cut=True):
	# Load reference signal and timestamps of rotating output polarizer
	# Recorded time is in seconds
	# Reference signal data should be contained in HDF5 file format

	#refkey = {}
	f = h5py.File(refpath,'r')
	refdata_all = f['EncoderData']['Raw'].value
	reftime = f['EncoderData']['Tsec'].value
	refdata = refdata_all & 0x2
	#reftime = refdata_all & 0x2	
	#reftime = refdata_all & 0x1
	refdata_all = np.NaN
	f = np.NaN
	#print reftime

	mask = int(80000*5)

	if not cut:
		print 'Not cutting initial data'
		return np.array(reftime), np.array(refdata,dtype='float64')
	else:
		print 'Cut initial 5 seconds of data: %i datapoints'%(mask)
		print len(np.array(reftime)), len(np.array(reftime)[mask:])
		return np.array(reftime)[mask:], np.array(refdata,dtype='float64')[mask:]


time_keys = {'y': '_irig_y',
           'd': '_irig_d',
           'h': '_irig_h',
           'm': '_irig_m',
           's': '_irig_s',
           'ss': '_irig_ss',
           't_s': '_irig_test_s',
           't_ss': '_irig_test_ss'
       	  }


def get_timekey(board,tstr):
	return '192_168_1_%i'%(board)+time_keys[tstr]


def get_timekey_gen(board,tstr):
        return board+time_keys[tstr]


def get_datakey(board,squid,chan):
	keybase = '192_168_1_%i_bolo_m%i_c%02i'%(board,squid,chan)
	return keybase,keybase+'_i',keybase+'_q'


def get_datakey_gen(board,squid,chan,ice=False):
	if ice:
        	keybase = '%s_m%i_bolo_c%02i'%(board,squid,chan)
	else:
        	keybase = '%s_bolo_m%i_c%02i'%(board,squid,chan)
        return keybase,keybase+'_i',keybase+'_q'

'''
def read_timedata(datapath,boloBoard,ice=False):
	# Load timestamps from specified DfMux board
	# Returns time in seconds
	# Convert timestamps into UTC time in units of seconds

	if ice:
        	dirfile = datapath
	else:
		dirfile = os.path.join(datapath,boloBoard)

	df = pygetdata.dirfile(dirfile)
	#ts = df.getdata(get_timekey(boloBoard,'t_s'),pygetdata.INT32,num_frames=df.nframes)
	#tss = df.getdata(get_timekey(boloBoard,'t_ss'),pygetdata.INT32,num_frames=df.nframes)
	y = df.getdata(get_timekey_gen(boloBoard,'y'),pygetdata.INT32,num_frames=df.nframes)
        d = df.getdata(get_timekey_gen(boloBoard,'d'),pygetdata.INT32,num_frames=df.nframes)
        h = df.getdata(get_timekey_gen(boloBoard,'h'),pygetdata.INT32,num_frames=df.nframes)
        m = df.getdata(get_timekey_gen(boloBoard,'m'),pygetdata.INT32,num_frames=df.nframes)
        s = df.getdata(get_timekey_gen(boloBoard,'s'),pygetdata.INT32,num_frames=df.nframes)
        ss = df.getdata(get_timekey_gen(boloBoard,'ss'),pygetdata.INT32,num_frames=df.nframes)

	times = []
	for i in range(len(y)):
		st = time.strptime('%d:%d:%d:%d:%d'%(y[i],d[i],h[i],m[i],s[i]),'%y:%j:%H:%M:%S')
		timestamp = calendar.timegm(st) + ss[i]*1e-8
		times.append(timestamp)
	
	return np.array(times)
'''
'''
def read_bolodata(datapath,boloBoard,ice=False,IQ=False,mmax=1,cnum=4):
	# Load bolometer data from DfMux board
	# Data is taken using DAN firmware that contains in-phase and quadrature components
	# Currently assuming the signal is just the magnitude of I and Q components
	# Need to understand better the mux demodulation to determine what is best

	if ice:
        	dirfile = datapath
		crange = range(cnum)
	else:
		dirfile = os.path.join(datapath,boloBoard)
		crange = range(1,cnum+1)
	mrange = range(1,mmax+1)

	df = pygetdata.dirfile(dirfile)
	names = []
	data = []
	data_I = []
	data_Q = []
	for i in mrange:
		for j in crange:
			datakey = get_datakey_gen(boloBoard,i,j,ice=ice)
			print datakey
			idata = df.getdata(datakey[1],pygetdata.INT32,num_frames=df.nframes)
			qdata = df.getdata(datakey[2],pygetdata.INT32,num_frames=df.nframes)
			if IQ:
				data_I.append(idata)
				data_Q.append(qdata)
			names.append(datakey[0])
			data.append(np.sqrt(idata**2.+qdata**2.))

	if IQ:
		return np.array(names),np.array(data_I),np.array(data_Q)
	else:
		return np.array(names),np.array(data)
'''
'''
def read_bolodata_2(datapath,IQ=False,sort=False):
        # Load bolometer data from ICEboard with phase2 data structure
        # Data is taken using DAN firmware that contains in-phase and quadrature components
        # Currently assuming the signal is just the magnitude of I and Q components
        # Need to understand better the mux demodulation to determine what is best

        d = netCDF4.Dataset(datapath)
	keys = d.variables.keys()
	keys.remove('Time')
	bolonames = [x[:-2] for x in keys if x[-1]!='Q']
	times = d.variables['Time'][:]

	if sort:
		csvfile = '/scratch2/scratchdirs/fmatsuda/FTS_KEK/phase2/FTSData/PB20_08_01.csv'
		lc,id,comb = np.genfromtxt(csvfile,dtype=str,skip_header=1,unpack=True)
		#combname = ['006/16/1/1/']
		#combname = ['006/16/2/2/']#,'006/16/1/1/']
		id_comb = []
		#for j in range(len(combname)):
		#	for i in range(len(comb)):
		#		if combname[j] in comb[i]:
		#			id_comb.append(id[i].replace('/','_'))
		id_comb.append('others_others.others-02-B.2')
		print 'Only using bolos: ', id_comb

        names = []
        data = []
        data_I = []
        data_Q = []
        for name in bolonames:
		if sort:
			if name not in id_comb:
				continue
		idata = d.variables[name+'_I'][:]
		qdata = d.variables[name+'_Q'][:]
		if IQ:
			data_I.append(idata)
			data_Q.append(qdata)
		names.append(name)
		data.append(np.sqrt(idata**2.+qdata**2.))

        if IQ:
                return np.array(names),np.array(times),np.array(data_I),np.array(data_Q)
        else:
                return np.array(names),np.array(times),np.array(data)
'''

def read_bolodata_3(datapath,wafer,combs=None):
	# Load bolometer data for run 38
	# Data is now in g3 format

	hwmpath = '/home/cmb/fmatsuda/PB2a_run38/'
	g3c = libg3.G3Compressed(datapath,hwm=hwmpath,loadtime=False)
	bolonames_all = g3c.bolonames_all
	wafer_header = 'PB20.'
	prop = g3c.boloprop

	find_wafer = wafer_header+wafer
	wafernames = [x for x in bolonames_all if os.path.split(x)[0]==find_wafer]

	bolonames = []
	if not combs:
		combs = range(1,29)
	for i in combs:
		combname = 'Comb%02d'%(i)
		bolonames.extend([x for x in wafernames if combname in os.path.split(x)[1]])
	print bolonames
	
	frequency = []
	for iname in bolonames:
		maps = re.split(r'\.|\s|_', "%s"%prop[iname])
		#print maps
		if len(maps) < 9: 
			frequency.append(np.NaN)
		else:
			frequency.append(maps[8])

	names = [os.path.split(name)[1] for name in bolonames]
	g3c.loadbolo(name=bolonames)
        data = g3c.bolo
	times = g3c.bolotime / 1e8

	#time_mjd = []
	#for j in range(len(times)):
	#	times_mjd.append(core.G3Time(times[j]).mjd)
	
	#return np.array(names),np.array(times_mjd),np.array(data)
	return np.array(names),np.array(times),np.array(data),np.array(frequency)


def get_steptimes(timepath,steplength=2.0):
	# Load and extract start and end times of each integration step (i.e. interferogram datapoint)
	# This timestamp should also be in UTC seconds

	with open(timepath,'r') as f:
		tpkl = pickle.load(f)
	starttimes = []
	endtimes = []
	#boardkey = 'Mb63_timestamp'
	boardkey = 'Mb0_timestamp'
	#print boardkey
	steptimes = tpkl[boardkey]	
	drop = np.zeros(len(steptimes),dtype=bool)

	for i in range(len(steptimes)):
		starttime_array = np.array(steptimes[i][0])
		endtime_array = np.array(steptimes[i][1])
		if np.std(starttime_array) > 1. :
			ok = ((starttime_array-starttimes[i-1]) < 2.*steplength)&((starttime_array-starttimes[i-1]) > 0.)
			starttime_array = starttime_array[ok]
		if np.std(endtime_array) > 1. :
			ok = ((endtime_array-endtimes[i-1]) < 2.*steplength)&((endtime_array-endtimes[i-1]) > 0.)
                        endtime_array = endtime_array[ok]
		start_avg = np.average(starttime_array)
		end_avg = np.average(endtime_array)
		if np.isnan(start_avg):
			start_avg = end_avg - steplength
			print 'Step %i start time is NaN. Defaulting to standard 2 seconds interval based off end time...'%(i)
		if np.isnan(end_avg):
			end_avg = start_avg + steplength
			print 'Step %i end time is NaN. Defaulting to standard 2 seconds interval based off start time...'%(i)
		#if (end_avg-start_avg) > 3.:
		#	end_avg = start_avg + 2.
		#	print 'Step %i is too long: %s. Defaulting to standard 2 seconds...'%(i,end_avg-start_avg)
		delta = end_avg-start_avg
		if (delta > steplength+1.) or (delta < steplength-1.):
			print 'Step %i is irregular length: %s. Flagging'%(i,delta)
			drop[i] = True
		starttimes.append(start_avg)
		endtimes.append(end_avg)
			
	return np.array(starttimes),np.array(endtimes),drop


def calc_refIQ(refsig,reftime,bolotime,harm=1,deltat=0.):
        # Calculate in-phase and quadrature components from chopper reference signal in order to demodulate signal
        print 'Adjusting quad pair for timing error of value %s sec '%(deltat)
	refIQ = GenQuadPair(refsig,reftime,bolotime,harm=harm,deltat=deltat)
        return refIQ


def GenQuadPair(ref,reftime,bolotime,harm=1,deltat=0.):
        sample_rate = 25e6/2.**17.
	#ref -= (np.max(ref)-np.min(ref))/2.
	#print np.max(ref), np.min(ref)
	print np.average(ref)
	ref -= np.average(ref)

        def model(ts,ampx,ampy,f,theta):
                phi = 2.0*np.pi*f*harm*ts + theta
                return ampx*np.cos(phi) + ampy*np.sin(phi)

	def model_cos(ts,ampx,ampy,f,theta):
                phi = 2.0*np.pi*f*harm*ts + theta
                return ampx*np.cos(phi) 

	ref_fit = ref[::100]
	reftime_fit = reftime[::100] - reftime[0]
	nt = len(ref_fit)

	nsplit = 5
        amp = []
        freq = []
        phase = []
        print 'Total length of downsampled reference = %i'%(nt)
        for i in range(nsplit):
                ibegin = i*nt/nsplit
                iend = (i+1)*nt/nsplit
                print 'Fitting reference signal range: %s : %s'%(ibegin,iend)
                popt_sec,pcov_sec = curve_fit(model_cos,reftime_fit[ibegin:iend],ref_fit[ibegin:iend],(0.,0.,4.0005,0.))
                if popt_sec[0] < 0.:
                        print '   Found negative amplitude, making positive and shifting phase by pi'
                        popt_sec[0] = np.absolute(popt_sec[0])
                        popt_sec[3] = popt_sec[3] - np.pi
                print '   Fitted parameters for range: ', popt_sec
                amp.append(popt_sec[0])
                freq.append(popt_sec[2])
                phase.append(popt_sec[3])
        popt = [np.average(amp),0.,np.average(freq),np.average(np.unwrap(phase))]
        print 'Averaged fitted parameters for reference signal: ',popt

        #print deltat
        bolotime += deltat
        #ref = model(bolotime-bolotime[0],*popt)
        #refx = model(ts,1.0,0.0,popt[2])
        #refy = model(ts,0.0,1.0,popt[2])
        refx = model(bolotime-bolotime[0],1.0,0.0,popt[2],popt[3])
        refy = model(bolotime-bolotime[0],0.0,1.0,popt[2],popt[3])

        return {'X':refx,'Y':refy}


def demod_data(df,step_idx,step_mask,refIQ,buff=0.,harm=1,chan_itr=None,save_dir=None,phase_opt=None,verbose=True):
	# Demodulate interferogram data

	nchans = len(df.bolo)
	data = np.zeros((nchans,step_idx.shape[0]))
	noise = np.zeros((nchans,step_idx.shape[0]))
	drops = np.zeros((nchans,step_idx.shape[0]),dtype=bool)
	Xout = np.zeros((nchans,step_idx.shape[0]))
	Yout = np.zeros((nchans,step_idx.shape[0]))	
	phase_ang = np.zeros(nchans)
	phase_ang[:] = np.NaN

	if (chan_itr is not None) and (phase_opt is not None):
		chans = [int(chan_itr)]
		phase_ang[int(chan_itr)] = phase_opt
	else:
        	chans = range(nchans)
		phase_ang[:] = phase_opt 

	bolo_times = df.bolo_time

    	for chan in chans:
		if verbose:
			print 'Demodulating channel %s'%(chan)

		ts = df.bolo[chan]

		data_step = ts[step_idx]
		refIQ_step = {'X':refIQ['X'][step_idx],'Y':refIQ['Y'][step_idx]} 

		if save_dir is not None:
			pl.plot(ts[::50])
			pl.savefig(os.path.join(save_dir,'plots/ts_%s.png'%(df.boloid[chan])))
			pl.clf()

		mask = step_mask 
		data_demod = lFA.Demod_chan(data_step,refIQ_step,mask,phase_opt=phase_ang[chan])
		
		data[chan][:] = data_demod['Mag']
		noise[chan][:] = data_demod['Noise']
		Xout[chan][:] = data_demod['X']
		Yout[chan][:] = data_demod['Y']

		#ndata = np.sum(mask,axis=1)
		#nzeros = np.sum(-mask,axis=1)
		#drops[chan][ndata<step_mask.shape[1]/2.] = True
		#drop_step = np.sum(drops[chan])

		#if verbose:
		#	print 'Dropped signal bins %i / %i'%(drop_step,step_idx.shape[0])
		if save_dir is not None:
			pl.plot(Xout[chan])
			pl.plot(Yout[chan])
			pl.savefig(os.path.join(save_dir,'plots/demodXY_%s.png'%(df.boloid[chan])))
			pl.clf()

    	return data,noise,chans,Xout,Yout


def demod_single(bolo_time,ts,step_idx,step_mask,refIQ,buff=0.,harm=1,phase_opt=None):
        # Demodulate interferogram data

        phase_ang = phase_opt

	data_step = ts[step_idx]
	refIQ_step = {'X':refIQ['X'][step_idx],'Y':refIQ['Y'][step_idx]}

	mask = step_mask
	data_demod = lFA.Demod_chan(data_step,refIQ_step,mask,phase_opt=phase_ang)

	data = data_demod['Mag']
	noise = data_demod['Noise']
	Xout = data_demod['X']
	Yout = data_demod['Y']

        return data,noise,Xout,Yout


def phase_iteration_minY(df,step_idx,step_mask,refIQ,harm=1,chan_itr=None):
        # Calculate necessary phase correction using iteration procedure and minimize quadrature component
        def resid_phase(c):
                data, noise, chans, X, Y = demod_data(df,step_idx,step_mask,refIQ,buff=0.,harm=harm,chan_itr=chan_itr,phase_opt=c,verbose=False)
                Y_chan = Y[chan_itr]
                Y_d = lFA.detrend(Y_chan,poly=12)
                return Y_d
        p0 = 0.
        p = leastsq(resid_phase,p0)
        print 'Optimized phase = %.3f \n'%(p[0][0]*180./np.pi)
        return p[0][0]


def phase_iteration_single_minY(df,step_idx,step_mask,refIQ,chan_itr,harm=1):
        # Calculate necessary phase correction using iteration procedure and minimize quadrature component
	bolo_time = df.bolo_time
	ts = df.bolo[chan_itr]
        def resid_phase(c):
                data, noise, X, Y = demod_single(bolo_time,ts,step_idx,step_mask,refIQ,buff=0.,harm=harm,phase_opt=c)
                Y_d = lFA.detrend(Y,poly=12)
                return Y_d
        p0 = 0.
        p = leastsq(resid_phase,p0)
        print 'Optimized phase = %.3f \n'%(p[0][0]*180./np.pi)
        return p[0][0]


def interp_data(data,drop):
	# Interpolate spectra where there was irregular data

	drop_idx = np.where(drop)[0]
	print drop_idx
	for idx in drop_idx:
		for chan in range(data.shape[0]):
			fit = np.average(np.array([data[chan,idx-1],data[chan,idx+1]]))
			data[chan,idx] = fit


def analyze_spectra(data,boloid,output_path,amode,freq):

	lFA.make_plots_dir(output_path) 
	lFA.plot_blank(output_path)	

	dataout = {}
	dataout['Interferograms'] = data
	spectra_list = []
	chans = range(data.shape[0])
	bands = np.empty(len(chans),dtype='S10')
	for f in range(len(freq)):
		if freq[f] == '90':
			bands[f] = '90GHz'
		elif freq[f] == '150':
			bands[f] = '150GHz'
		else:
			bands[f] = '150GHz'

	atmData = lFA.getAtmData('/home/cmb/fmatsuda/FTS_run38/AtmosphericModel/Chaj_zadegrees30_pwvmicrons1000.out')

	for i in chans:
		energy = np.std(lFA.detrend(data[i]))
		if energy < 0.25: 
			print '%s energy less than 0.25'%(boloid[i])
			spectra_list.append(None)
			continue
		max_data = np.max(data[i])
		if np.isnan(max_data): 
			print '%s NaN max with %s'%(boloid[i],len(np.isnan(data[i])))
			spectra_list.append(None)
			continue

		spectral_data = lFA.fts_spectrum(data[i],boloid[i],'triangle',plot_path=output_path,mode=amode,band=bands[i])

		lFA.plot_spectrum(spectral_data,'o-')
		lFA.plot_conf()
		pl.savefig(os.path.join(output_path,'plots/fts-%s_points.png'%(boloid[i])))
		pl.clf()

		lFA.plot_spectrum(spectral_data)
		pl.plot(atmData['Frequencies'],atmData['Transmission'],linestyle = '--',color = 'r')
		#pl.xlim([0,250])
		pl.ylim([-0.25,1.0])
		lFA.plot_stats(spectral_data)
		lFA.plot_conf()
		pl.savefig(os.path.join(output_path,'plots/fts-%s.png'%(boloid[i])))
		pl.clf()

		bw = spectral_data['frequencies'][spectral_data['left_edge']] - spectral_data['frequencies'][spectral_data['right_edge']]
		spectra_list.append(spectral_data)

	dataout['Spectra'] = spectra_list

	return dataout


def getargs():

	parser = argparse.ArgumentParser(description='Demodulate and analyze FTS data and calculate spectra')
	parser.add_argument('-i', dest='inputfn', action='store', help='Specify input data directory')
	parser.add_argument('-t', dest='timefn', action='store', help='Specify input time pkl file')
	parser.add_argument('-r', dest='reffn', action='store', help='Reference signal directory')
	parser.add_argument('-x', dest='hwmapfn', action='store', help='Hwmap file')
	parser.add_argument('-d', dest='demodfn', action='store', default=None, help='Demodulated data file')
	parser.add_argument('-o', dest='outputdir', action='store', default='.', help='Specify output directory')
	parser.add_argument('-hm', dest='harmonic', action='store', default=1, type=int, help='Specify harmonic to demodulate data with')
	parser.add_argument('-popt', dest='popt', action='store_true', default=False, help='Specify to do phase optimization in demodulation')
	parser.add_argument('-pfn', dest='pfn', action='store', default=None, help='Previously calculated optimized phase file')
	parser.add_argument('-g', dest='gain', action='store_true', default=False, help='Calculate relative gain between pixel pairs')
	#parser.add_argument('-ice', dest='ice', action='store_true', default=False, help='Use ICEBoards')
	#parser.add_argument('-boardIP', dest='boardIP', action='store', help='IP address of board to analyze')
	parser.add_argument('-inttime', dest='inttime', action='store', default=2.0, type=float, help='Integration time for FTS step')
	args = parser.parse_args()

	return args


def main(args):

	pl.figure()

	tstr = os.path.split(args.timefn)[1][14:-4]

	if args.demodfn is None:
		lFA.make_plots_dir(args.outputdir)
		print 'Loading all necessary files...'
		reftimes, refsig = read_refdata(args.reffn)
		print 'Loaded reference data.'
		waferstr = '11.08'
		comblist = range(13,29)
		boloid,bolo_time,ts,freq = read_bolodata_3(args.inputfn,waferstr,combs=comblist)
		print boloid, bolo_time.shape,ts.shape
		print 'Loaded bolometer data.'
		starttimes,endtimes,drops = get_steptimes(args.timefn,steplength=args.inttime)
		print 'Loaded time pickle data.'

		pl.plot(endtimes-starttimes,'o')
		pl.savefig(os.path.join(args.outputdir,'plots/steptimes.png'))
		pl.clf()


		print len(refsig), len(reftimes)
		print ts.shape,len(bolo_time)

		print starttimes[0], endtimes[-1]
		print bolo_time[0], bolo_time[-1]
		print reftimes[0], reftimes[-1]

		# Cut data where step data began and ended to make data size smaller
		cutbuff = 1.0	# sec
		refcut = (reftimes > starttimes[0]-cutbuff)&(reftimes < endtimes[-1]+cutbuff)
		bolocut = (bolo_time > starttimes[0]-cutbuff)&(bolo_time < endtimes[-1]+cutbuff)

		reftimes = reftimes[refcut]
		refsig = refsig[refcut]
		bolo_time = bolo_time[bolocut]
		ts = ts[:,bolocut]

		print len(refsig), len(reftimes)
                print ts.shape,len(bolo_time)

		print 'Step: ',starttimes[0]
                print 'Ref: ',reftimes[0]
                print 'Bolo: ',bolo_time[0]

		class df: pass
		df.bolo_time = bolo_time
		df.bolo = ts
		#df.boloid = boloid
		boloid2 = []
		for i in range(len(boloid)):
			boloid2.append(boloid[i]+'_'+freq[i])
		df.boloid = np.array(boloid2)
		print df.boloid

		skip_idx = []
		skip_num = []
		delta_t = bolo_time[0] - reftimes[0]

		step_idx, step_mask = lFA.step_idx(df,starttimes,endtimes,steplength=args.inttime)
		print 'Fitting for reference signal IQ...'
		refIQ = calc_refIQ(refsig,reftimes,df.bolo_time,harm=args.harmonic,deltat=delta_t)
		#refIQ = lFA.calc_refIQ_value(df.bolo_time,skip_idx,skip_num,harm=args.harmonic,deltat=delta_t)

		downsample = False
		if downsample:
			print 'Downsampling only for phase optimization'
			downsample = 2
			class df_opt: pass
			df_opt.bolo_time = bolo_time[::downsample]
			df_opt.bolo = ts[:][::downsample]
			df_opt.boloid = df.boloid
			print 'Downsampling: ',len(df.bolo_time),len(df_opt.bolo_time)
			step_idx_opt, step_mask_opt = lFA.step_idx(df_opt,starttimes,endtimes,steplength=args.inttime)
			refIQ_opt = calc_refIQ(refsig,reftimes,df_opt.bolo_time,harm=args.harmonic,deltat=delta_t)

		if args.popt:
			if args.pfn is None:
				phases = np.zeros(len(df.bolo))
				for i in range(len(df.bolo)):
					print 'Determining phase through phase optimization of channel %s'%(df.boloid[i])
					if downsample:
						phase_opt = phase_iteration_single_minY(df_opt,step_idx_opt,step_mask_opt,refIQ_opt,i,harm=args.harmonic)
					else:
						phase_opt = phase_iteration_minY(df,step_idx,step_mask,refIQ,harm=args.harmonic,chan_itr=i)
						#phase_opt = phase_iteration_single_minY(df,step_idx,step_mask,refIQ,i,harm=args.harmonic)
					phases[i] = phase_opt
				print 'Saving phase optimization data...'
				phase_dict = {'BoloID':boloid,'Phases':phases}
				with open(os.path.join(args.outputdir,'Phase_opt-%s.pkl'%(tstr)),'w') as f:
					pickle.dump(phase_dict,f)
			else:
				print 'Loading phase optimization data...'
				with open(args.pfn,'r') as f:
					phase_dict = pickle.load(f)
				phases = phase_dict['Phases']
			print 'Demodulating each step with optimized phase using harmonic %s...'%(args.harmonic)
			data, noise, chans, X, Y = demod_data(df,step_idx,step_mask,refIQ,buff=0.,harm=args.harmonic,save_dir=args.outputdir,phase_opt=phases)
		else:
			print 'Demodulating each step using harmonic %s...'%(args.harmonic)
			data, noise, chans, X, Y = demod_data(df,step_idx,step_mask,refIQ,buff=0.,harm=args.harmonic,save_dir=args.outputdir)

		print 'Saving demodulated data...'	
		with open(os.path.join(args.outputdir,'Demod-%s.pkl'%(tstr)),'w') as f:
			pickle.dump(data,f)
	else:
		with open(args.demodfn,'r') as f:
			data = pickle.load(f)

	#print data
	print 'Filling in data with interpolation for steps with irregular lengths...'
	interp_data(data,drops)

	print 'Start full spectral analysis...'
	print df.boloid
	outputdata = analyze_spectra(data,df.boloid,args.outputdir,'d',freq)

	pl.close()

	print 'Saving pkl file for raw data outputs...'
	with open(os.path.join(args.outputdir,'FTS_data-%s.pkl'%(tstr)),'w') as f:
		pickle.dump(outputdata,f)


if __name__=="__main__":

	print 'Start FTS analysis'
	args = getargs()
	main(args)
