####################################################################################################
#
# FTS Scan for UCSD
#
# Lindsay Ng Lowry
# llowry@ucsd.edu
# 20190220
#
# This program was modified from the FTSScan.py code used in KEK during Run38 testing of PB2a
# (written by Fred Matsuda) to work with the FTSs at UCSD.  It runs an FTS scan using a specified
# integration time and linear stage steps indicated by the FTS linear stage map file, defaulted to
# be named 'FTSScanPoints.txt', and saves data to a specified file.
#
#
# Run this code as:
# python FTSScan_UCSD.py -hwm '[path to hwm yaml file]' -f '[name of map file with file extension]'
# -int [integration time in seconds as float] -center
#
# mapfile should be saved in: '/home/polarbear/FTS_ControlComputer/FTS/MapFiles'
# data file will be saved in: '/home/polarbear/FTS_Data'
# data file will be named: '[timestamp at end of scan']_InterferogramPosData.pkl'
#
# Data file is a dictionary with keys 'Position' and 'Mb0_timestamp', where each is an array with
# the same number of elements.  'Position' contains the linear stage positions in encoder counts (as
# in the mapfile).  Each element in 'Mb0_timestamp' contains 2 arrays, the first with timestamps
# recorded at the start of data taking for that position, and the second with timestamps recorded at
# the end of data taking for that position.
#
####################################################################################################


#!/usr/bin/env python

import os
import calendar
import sys, time, numpy, pylab
import argparse
import cPickle as pickle
import pygetdata
import pydfmux

sys.path.append('/home/polarbear/FTS')
import FTSControl_UCSD


def LoadBoard(hwm_loc,Board=int(1)):
	"""Get a dfmux (iceboard) object for a particular board from a pydfmux hwm query.

	Parameters:
	hwm_loc (str) -- the path to the hardware map yaml file
	Board (int) -- the integer indicating the desired iceboard according to the hardware map (default 1)

	Returns:
	b (dfmux object) -- the desired dfmux object
	"""
	hwm_location = hwm_loc
	try:
		session = pydfmux.load_session(file(hwm_location))
		hwm = session['hardware_map']
		ds = hwm.query(pydfmux.Dfmux)
		b = ds[Board]
	except:
		print ('\n\tSomething went wrong with loading the hwm.  Exiting.')
		exit()

	return b

def ReadTime(b):
	"""Read the time from an iceboard converted using timegm.

	Parameters:
	b (dfmux object) -- the dfmux object corresponding to the desired iceboard

	Returns:
	tepoch (float) -- the timestamp of the iceboard in seconds since 1970???
	"""
	t = b.get_timestamp()
	d = t.d
	st = time.strptime('%d:%d:%d:%d:%d'%(d['y'],d['d'],d['h'],d['m'],d['s']),'%y:%j:%H:%M:%S')
	tepoch = calendar.timegm(st)
	tepoch += d['ss']*1.e-8

	return tepoch


def FTS_Scan(dataname,hwm_loc,mapfile,stageIP,stagePort,vel=1.0,accel=2.0,intTime=3,pauseTime=0.5,save=None):
	"""Perform an FTS Scan using specified linear stage positions.

	Using mapfile, retrieve a list of specified linear stage positions.  Move the stage to each
	position, waiting between each motion to collect bolometer data (saved/recorded elsewhere) and
	storing timestamp and linear stage position data.  Return either a dictionary containing the
	name of the pkl file to which the data was saved, or a dictionary containing the data.

	Parameters:
	dataname (str) -- THIS DOES NOT SEEM TO BE USED
	hwm_loc (str) -- the full path to the hwm yaml file
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

	# Create a dictionary to hold the stage position and timestamp data
	BeamMap={}
	# Create a list to hold the channels (boards) from which you are getting the timestamps
	# (in this case using only one board so list probably isn't necessary)
	chList = []

	# Load an iceboard from the hardware map (second board in hwm?)
	b = LoadBoard(hwm_loc)

	# Set up the dictionary to hold timestamps (from each board being used) and stage positions
	chList.append('Mb'+str(0))
	print '\n\tRecording timestamp from boards ' + str(0)
	for ch in chList:
		BeamMap[ch+'_timestamp']=[]
	BeamMap['Position'] = []



	# Get the map points from mapfile and print the integration time
	mappoints = numpy.loadtxt(mapfile)
	print '\n\tIntegration time set to %s sec'%(intTime)



	# Initialize the linear stage motor FTSControl object
	StageMotor=FTSControl_UCSD.FTSControl(stageIp=stageIP,stagePort=stagePort) 
	StageMotor.setVelocity(motor=1,velocity=vel)
	StageMotor.setAcceleration(motor=1,accel=accel)
	StageMotor.setMotorEnable(motor=1,True)



	try: #allow for keyboard interrupt
		# For each map point in the mapfile...
		print ('\n\tBeginning scan')
		for i in range(len(mappoints)):

			# Move the linear stage to the specified position
	   		#StageMotor.setMotorEnable(1,True)
	   		print ('\t\tPosition: ' + str(mappoints[i]) + ' steps' )
	   		StageMotor.moveAxisToPosition(motor=1,mappoints[i],linStage=True)

	   		StageMotor.blockWhileMoving(motor=1)
	   		#StageMotor.setMotorEnable(1,False)


	   		time.sleep(pauseTime)		# Wait for transients to die down


	   		# Set up lists to hold the start and end timestamps for this position
			timestamp_start = []
			timestamp_end = []
			loopnum = 10

			# Read the timestamp from the iceboard before taking data
				# Apparently errors can occur when fetching timestamps using this setup
				# Just to be sure get at least 5 timestamps each time
			for m in range(loopnum):
				timestamp_start.append(ReadTime(b))

			# Wait the specified integration time so data can be taken
			time.sleep(intTime)		# DATA TAKEN HERE

			# Read the timestamp from the iceboard after taking data
			for n in range(loopnum):
				timestamp_end.append(ReadTime(b))

			print ('\t\tStart and End timestamps: ')
			print timestamp_start, timestamp_end

			# Store the timestamp and position data to the dictionary
			BeamMap[chList[0]+'_timestamp'].append([timestamp_start,timestamp_end])
			BeamMap['Position'].append(mappoints[i])

	except KeyboardInterrupt:
		print 'Keyboard Interrupt found.  Exiting mapfile loop.  Last position was '+ str(mappoints[i]) + '.  Saving data acquired so far.'
		#continue

	finally:
		# Save the timestamp and position data from the full scan to a file with filename indicated by current gmtime
		if save is not None:
			tstr = time.strftime('%Y%m%d_%H%M%S',time.gmtime(time.time()))
			print ('\n\tTimestamp for data file: ' + tstr)
			print ('\n\tSaving %s_InterferogramPosData.pkl'%(tstr))

			with open(os.path.join(save,'%s_InterferogramPosData.pkl'),'w') as f:
				pickle.dump(BeamMap,f)

		# Return the linear stage to its zero position
		#StageMotor.setMotorEnable(1,True)
		StageMotor.moveAxisToPosition(1,0,linStage=True)
		StageMotor.blockWhileMoving(1)
		#StageMotor.setMotorEnable(1,False)

		# close the connection to the linear stage motor
		StageMotor.closeConnection(1)

		# Return either the full timestamp and position dataset (if the data was not saved) or TimeMap which indicates the filename
		if save is not None:
			TimeMap = {}
			TimeMap['Time_pkl'] = tstr+'_InterferogramPosData.pkl'
			return TimeMap
		else:
			return BeamMap


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-d',dest='dataname',action='store',default=None,help='Name of dataset')
	parser.add_argument('-hwm',dest='hwm_loc',action='store',default=None,help='Path to hardware map yaml file')
	parser.add_argument('-f',dest='ftsmapfn',action='store',default='FTSScanPoints.txt',help='FTS linear stage map file name')
	parser.add_argument('-int',dest='inttime',action='store',type=int,default=3,help='Integration time per step in seconds')
	parser.add_argument('-center',dest='center',action='store_true',help='Center linear stage before scan')
	parser.add_argument('-test',dest='test',action='store_true',help='Test run (0.5s integration time)')
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

	# Exit if no hwm is given
	if hwm_loc == None:
		print ('\n\tNo hardware map given.  Exiting.')
		exit()


	# Get the IP addresses and port numbers for the two motors
	Mirror_Pol_IP = '192.168.1.113'
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
		StageMotor.setMotorEnable(motor=1,True)

		# Command the stage motor to go to its home position
		StageMotor.seekHomeLinearStage(motor=1)

		# Move the stage to approximately the middle
		StageMotor.moveAxisByLength(motor=1,pos=9.0,posIsInches=True,linStage=True)
		StageMotor.blockWhileMoving(1)

		# Set this position to be the temporary zero point
		StageMotor.setZero(motor=1)
		newpos = StageMotor.getPosition(motor=1)[0]
		if newpos !=0:
			print ('\n\tStage not zeroed after initial centering.  Exiting.')
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
			print ('\n\tSomething went wrong with centering the stage - stage not zeroed.  Current position is %i (in counts).  Exiting.' % newpos)
			StageMotor.closeConnection()
			exit()
		else:
			print ('\n\tStage centered and zeroed.  Closing connection.')
			StageMotor.closeConnection()


	# Initialize the output polarizer FTSControl object
	PolRotator = FTSControl_UCSD.FTSControl(polarizerIp=Pol_Rot_IP,polarizerPort=Pol_Rot_Port)

	# Start up output polarizer rotation
	vel_pol = 12.0      # 3.0 = 1Hz, 6.0 = 2Hz, ...
	accel_pol = 1.0

	print '\n\tStarting up output rotating polarizer'
	PolRotator.setMotorEnable(motor=2,enable=True)
	PolRotator.startRotation(motor=2,velocity=vel_pol,accel=accel_pol)
	time.sleep(20)

	# Run the FTS Scan and save file to save_path using default values for velocity, acceleration, and pauseTime
	data = FTS_Scan(args.dataname,map_file,Mirror_Stage_IP,Mirror_Stage_Port,intTime=intTime,save=save_path)

	# Stop the output polarizer rotation once the FTS Scan is complete
	print '\n\tStopping output rotating polarizer'
	PolRotator.stopRotation(axis=1)
	PolRotator.closeConnection(1)


	print '\n\tFinished interferogram scan'

