####################################################################################################
#
# FTS Scan for UCSD Testing with Nate's Cryostat
#
# Lindsay Ng Lowry
# llowry@ucsd.edu
# 20180827
#
# This program was modified from the FTSScan_UCSD.py code to be used with the UCSD FTS during a test
# scan using Nate's cryostat at UCSD.  It runs an FTS scan using a specified integration time and 
# linear stage steps indicated by the FTS linear stage map file, defaulted to be named
# 'FTSScanPoints.txt', but does NOT spin the output polarizer.  The test was performed by Lindsay
# and Praween in 201808.
#
# This code uses Praween's LabJack for getting the timestamp and the bolometer data, and saves this
# along with the linear stage position to a specified file.  It runs an FTS scan using a specified
# integration time and linear stage steps indicated by the FTS linear stage map file, defaulted to
# be named 'FTSScanPoints.txt', and saves data to a specified file.
#
####################################################################################################


#!/usr/bin/env python

import os
import calendar
import sys, time, numpy, pylab
import argparse
import cPickle as pickle
#import pygetdata
#import pydfmux

sys.path.append('/home/polarbear/FTS_ControlComputer/FTS')
import FTSControl_UCSD
from T7python import T7


def LoadBoard(Board=int(1)):
	"""Dummy function for loading an iceboard.  Returns string 'DUMMY ICEBOARD'."""
	return 'DUMMY ICEBOARD'

def ReadTime(b):
	"""Read the time from an iceboard converted using timegm.

	Parameters:
	b (dfmux object)-- the dfmux object corresponding to the desired iceboard

	Returns:
	tepoch (float) -- the timestamp of the iceboard in seconds since 1970???
	"""
	t = b.get_timestamp()
	d = t.d
	st = time.strptime('%d:%d:%d:%d:%d'%(d['y'],d['d'],d['h'],d['m'],d['s']),'%y:%j:%H:%M:%S')
	tepoch = calendar.timegm(st)
	tepoch += d['ss']*1.e-8

	return tepoch


def FTS_Scan(dataname,mapfile,stageIP,stagePort,vel=1.0,accel=2.0,intTime=3,pauseTime=0.5,save=None):
	"""Perform an FTS Scan using specified linear stage positions.

	Using mapfile, retrieve a list of specified linear stage positions.  Move the stage to each
	position, waiting between each motion to collect bolometer data (saved/recorded elsewhere) and
	storing timestamp and linear stage position data.  Return either a dictionary containing the
	name of the pkl file to which the data was saved, or a dictionary containing the data.

	Parameters:
	dataname (str) -- 
	mapfile (str) -- the filename containing the linear stage scan points
	stageIP (str) -- the IP address associated with the linear stage motor
	stagePort (int) -- the port address associated with the linear stage motor
	vel (float) -- the velocity for the linear stage motor (default 1.0)
	accel (int) -- the acceleration for the linear stage motor (default 2.0)
	intTime (float) -- the integration time in seconds for each linear stage position (default 3)
	pauseTime (float) -- time in seconds to wait between moving the stage and recording the timestamps (default 0.5)
	save (str) -- the path to the directory in which to save the data file, with None indicating no data is saved (default None)

	Returns:
	TimeMap (if save is not None) -- dictionary indicating the pkl file to which the data was saved
	BeamMap (if save is None) -- dictionary indicating the timestamp and position data
	"""


	"""
	BeamMap={}
	chList = []

	# Load an iceboard from the hardware map (second board in hwm?)
	b = LoadBoard()
	"""
	

	# **********************************************************************************************
	# Set up a dictionary to hold the data and timestamps for each linear stage position
	dataDict = {}

	# Set the parameters for the labjack DAQ
	scanRate = 1000 # (roughly) number of samples labjack will take per channel in 1 second
	aScanListNames = ['AIN0'] # the port on the labjack you want data from
	aNames = ['AIN0_RANGE', 'AIN0_SETTLING_US','AIN0_RESOLUTION_INDEX']
	aValues = [10.0, 0, 0]
	
	maxRequest = intTime*2

	# Initialize the labjack
	labjack = T7()
	# **********************************************************************************************
	
	"""
	# ????
	chList.append('Mb'+str(0))
	print 'Recording timestamp from boards ' + str(0)
	for ch in chList: BeamMap[ch+'_timestamp']=[]
	BeamMap['Position'] = []
	"""

	# Get the map points from mapfile IN REVERSE ORDER and print the integration time
	mappoints = numpy.loadtxt(mapfile) # generated from WriteFTSScanFile to go from positive to negative
	mappoints = numpy.array(sorted(mappoints)) # reverse the order so that it scans from negative to positive
	print 'Integration time set to %s sec'%(intTime)

	# Initialize the linear stage motor FTSControl object
	StageMotor=FTSControl_UCSD.FTSControl(stageIp=stageIP,stagePort=stagePort) 
	StageMotor.setVelocity(motor=1,velocity=vel)
	StageMotor.setAcceleration(motor=1,accel=accel)
	StageMotor.setMotorEnable(motor=1,enable=True)

	# For each map point in the mapfile...
	for i in range(len(mappoints)):

		# Move the linear stage to the specified position
   		#StageMotor.setMotorEnable(1,True)
   		print ('  Position : ' + str(mappoints[i]) + ' steps' )
   		StageMotor.moveAxisToPosition(motor=1,pos=mappoints[i],linStage=True)

   		StageMotor.blockWhileMoving(motor=1)
   		#StageMotor.setMotorEnable(1,False)


   		time.sleep(pauseTime)		# Wait for transients to die down


   		"""
   		# Read the timestamp from the iceboard before taking data
			# Apparently errors can occur when fetching timestamps using this setup
			# Just to be sure get at least 5 timestamps each time
		timestamp_start = []
		timestamp_end = []
		loopnum = 10

		for m in range(loopnum):
			timestamp_start.append(ReadTime(b))

		# Wait the specified integration time so data can be taken
		time.sleep(intTime)		# Data taken here

		# Read the timestamp from the iceboard after taking data
		for n in range(loopnum):
			timestamp_end.append(ReadTime(b))


		print timestamp_start, timestamp_end

		# Store the timestamp and position data
		BeamMap[chList[0]+'_timestamp'].append([timestamp_start,timestamp_end])
		BeamMap['Position'].append(mappoints[i])
		"""

		# ******************************************************************************************
		# Get the bolometer and timestamp data from the labjack
		data = labjack.stream_array(maxRequest, scanRate, aScanListNames, aNames, aValues)
		labjack.stopstream()

		# Save the data to a dictionary keyed by linear stage position
		dataDict[mappoints[i]] = data
		# ******************************************************************************************




	# Save the timestamp and position data from the full scan to a file with filename indicated by current gmtime
	if save is not None:
		tstr = time.strftime('%Y%m%d_%H%M%S',time.gmtime(time.time()))
		print tstr
		"""
		print '  Saving Interferogram_%s.pkl'%(tstr)
		with open(os.path.join(save,'Interferogram_'+tstr+'.pkl'),'w') as f:
			pickle.dump(BeamMap,f)
		"""
		# ******************************************************************************************
		print ' Saving %s_Interferogram_NateCryostat.pkl'%(tstr)
		with open(os.path.join(save, tstr+'_Interferogram_NateCryostat.pkl'), 'w') as f:
			pickle.dump(dataDict,f)
		# ******************************************************************************************

	# Return the linear stage to its zero position
	#StageMotor.setMotorEnable(1,True)
	StageMotor.moveAxisToPosition(1,0,linStage=True)
	StageMotor.blockWhileMoving(1)
	#StageMotor.setMotorEnable(1,False)

	# close the connection to the linear stage motor
	StageMotor.closeConnection(1)


	"""
	# Return either the full timestamp and position dataset (if the data was not saved) or TimeMap which indicates the filename
	if save is not None:
		TimeMap = {}
		TimeMap['Time_pkl'] = 'Interferogram_'+tstr+'.pkl'
		return TimeMap
	else:
		return BeamMap
	"""

	# **********************************************************************************************
	# Return either the full dataset as a dictionary (if the data was not saved) or a string indicating the filename to which the data was saved
	if save is not None:
		saved_filename = tstr+'_Interferogram_NateCryostat.pkl'
		return saved_filename
	else:
		return dataDict
	# **********************************************************************************************



if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-d',dest='dataname',action='store',default=None,help='Name of dataset')
	parser.add_argument('-f',dest='ftsmapfn',action='store',default='FTSScanPoints.txt',help='FTS linear stage map file name')
	parser.add_argument('-int',dest='inttime',action='store',type=int,default=3,help='Integration time per step in seconds')
	parser.add_argument('-center',dest='center',action='store_true',help='Center linear stage before scan')
	parser.add_argument('-test',dest='test',action='store_true',help='Test run')
	args = parser.parse_args()

	# Set the integration time per step assuming 0.5s if specified as test run, using -int argument otherwise
	if args.test:
		intTime = 0.5
	else:
		intTime = args.inttime


	# Set the file paths for where to save the data and from where to access the linear stage map file
	save_path = '/home/polarbear/FTS_Data'

	map_dir = '/home/polarbear/FTS_ControlComputer/FTS/MapFiles'
	map_file = os.path.join(map_dir,args.ftsmapfn)


	# Get the IP addresses and port numbers for the two motors
	Mirror_Pol_IP = '192.168.1.103'
	Mirror_Stage_IP = Mirror_Pol_IP
	Mirror_Stage_Port = 4001
	Pol_Rot_IP = Mirror_Pol_IP
	Pol_Rot_Port = 4002

	# Center the linear stage if specified
	if args.center:
		# Initialize the linear stage motor
		StageMotor=FTSControl_UCSD.FTSControl(stageIp=Mirror_Stage_IP,stagePort=Mirror_Stage_Port) 
		StageMotor.setVelocity(motor=1,velocity=1.0)
		StageMotor.setAcceleration(motor=1,accel=2.0)
		StageMotor.setMotorEnable(motor=1,enable=True)

		# Command the stage motor to go to its home position
		StageMotor.seekHomeLinearStage(motor=1)

		# Move the stage to approximately the middle
		StageMotor.moveAxisByLength(motor=1,pos=9.0,posIsInches=True,linStage=True)    # NEED TO CHANGE THIS TO BE EXACT ONCE WE GET AN ALIGNMENT
		StageMotor.blockWhileMoving(1)

		# Set this position to be the temporary zero point
		StageMotor.setZero(motor=1)
		newpos = StageMotor.getPosition(motor=1)[0]
		if newpos !=0:
			print ('Not zeroed after initial centering.  Exiting.')
			StageMotor.closeConnection()
			exit()

		# Move the stage to the middle as determined by interferogram with Nate's Cryostat
		StageMotor.moveAxisToPosition(motor=1,pos=-6000,linStage=True)
		StageMotor.blockWhileMoving(1)

		# Set this position to be the zero point
		StageMotor.setZero(motor=1)

		# Check that the position is zero
		newpos = StageMotor.getPosition(motor=1)[0]
		if newpos != 0:
			print ('Not Zeroed')
			print newpos
			print ('Something went wrong with centering the stage - stage not zeroed.  Current position is %i (in counts).  Exiting.' % newpos)
			StageMotor.closeConnection()
			exit()
		else:
			"""
			print ('Stage zeroed.  Closing connection while output polarizer starts.')
			"""
			# **************************************************************************************
			print ('Stage zeroed.  Closing connection.')
			StageMotor.closeConnection()
			# **************************************************************************************



	"""
	# Initialize the output polarizer FTSControl object
	PolRotator = FTSControl_UCSD.FTSControl(polarizerIp=Pol_Rot_IP,polarizerPort=Pol_Rot_Port)

	# Start up output polarizer rotation
	vel_pol = 12.0      # 3.0 = 1Hz, 6.0 = 2Hz, ...
	accel_pol = 1.0

	print ' Starting up output rotating polarizer'
	PolRotator.setMotorEnable(motor=2,enable=True)
	PolRotator.startRotation(motor=2,velocity=vel_pol,accel=accel_pol)
	time.sleep(20)
	"""

	# Run the FTS Scan and save file to save_path using default values for velocity, acceleration, and pauseTime
	data = FTS_Scan(args.dataname,map_file,Mirror_Stage_IP,Mirror_Stage_Port,intTime=intTime,save=save_path)

	"""
	# Stop the output polarizer rotation once the FTS Scan is complete
	print ' Stopping output rotating polarizer'
	PolRotator.stopRotation(motor=1)
	PolRotator.closeConnection(1)
	"""

	print 'Finished interferogram scan'

