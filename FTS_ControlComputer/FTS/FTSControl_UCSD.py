####################################################################################################
#
# FTS Control for UCSD
#
# Kevin Crowley and Lindsay Ng Lowry
# llowry@ucsd.edu
# 20180320
#
# This program was modified from the FTSControl.py code used in KEK during Run38 testing of PB2a
# (written by Fred Matsuda) to work with the FTSs at UCSD.  It sets up the FTSControl objects that
# are used to control the two motors (output polarizer and linear stage) that are used by the FTS.
# This assumes the motor controllers have been configured properly using the ST Configurator Windows
# application.
#
# NEED TO DOCUMENT HOW THE MOTOR CONTROLLERS SHOULD BE CONFIGURED.
#
# Commands for communicating with the motor controllers can be found here:
# https://www.applied-motion.com/sites/default/files/hardware-manuals/Host-Command-Reference_920-0002P.PDF
#
####################################################################################################

import sys
sys.path.append('/Users/lindsay/Documents/UCSD/POLARBEAR_Research/FourierTransformSpectrometer/Code/LabTesting/FromPraween_2017')
from FTSMoxaSerial import Serial_TCPServer
from time import sleep
from pylab import load

# Time to wait for the user to power on the controllers and to see the power-on signature from the serial port
DEFAULT_WAIT_START_TIME = 15.0 #seconds

# Motor Names
STAGE = 1
POLARIZER = 2
ALL = 3

# Conversions for stages
AXIS_THREADS_PER_INCH_STAGE = 10.0 #Conversion for the FTS Linear Stage - Check for factor of two later
AXIS_THREADS_PER_INCH_XYZ = 10.0 #Measured on stepper, Check for factor of two later
#AXIS_THREADS_PER_INCH = 10.0 #Measured on stepper
MR_CODE_TO_STEPS_PER_REV = {
	'3': 2000.0,
	'4': 5000.0,
	'5': 10000.0,
	'6': 12800.0,
	'7': 18000.0,
	'8': 20000.0,
	'9': 21600.0,
	'10': 25000.0,
	'11': 25400.0,
	'12': 25600.0,
	'13': 36000.0,
	'14': 50000.0,
	'15': 50800.0
}




class FTSControl(object):
	"""FTSControl object for controlling up to 2 motors - linear stage and output polarizer"""
	def __init__(self, stageIp=None, stagePort=None, polarizerIp=None, polarizerPort=None, mRes=False):
		"""Initialize an FTSControl object.
		
		Parameters:
		
		stageIp (str) -- the IP address associated with the linear stage motor
		stagePort (int) -- the port address associated with the linear stage motor
		polarizerIp (str) -- the IP address associated with the output polarizer motor
		polarizerPort (int) -- the port address associated with the output polarizer motor
		mRes (bool) -- True if manual resolution, False if default (res=8) ???
		"""	

		# Set up the connection to the stage motor
		# Initialized so that the startup position is set to zero
		if not (stageIp and stagePort):
			print("Invalid stage connection information.  No stage control.")
			self.stage = None
		else:
			self.stage = Serial_TCPServer((stageIp, stagePort))
			if mRes:
				self.stage.propDict = {
					'name': 'stage',
					'motor': STAGE,
					'res': 'manual', 
					'pos': 0, #Position in counts (should always be an integer)
					'realPos': 0.0, #Position in inches
					'sPRev': 8000.0 #Steps per revolution (thread)
				}
			else:
				self.stage.propDict = {
					'name': 'stage',
					'motor': STAGE,
					'res': '8', #Corresponds to mapping above
					'pos': 0, #Position in counts (should always be an integer)
					'realPos': 0.0, #Position in inches
					'sPRev': MR_CODE_TO_STEPS_PER_REV['8'] #Steps per revolution (thread)
				}

		# Set up the connection to the polarizer motor
		# Initialized so that the startup position is set to zero
		if not (polarizerIp and polarizerPort):
			print("Invalid polarizer connection information.  No polarizer control.")
			self.polarizer = None
		else:
			self.polarizer = Serial_TCPServer((polarizerIp, polarizerPort))
			if mRes:
				self.polarizer.propDict = {
					'name': 'polarizer',
					'motor': POLARIZER,
					'res': 'manual', 
					'pos': 0, #Position in counts
					'realPos': 0, #Position in inches
					'sPRev': 8000.0 #Steps per revolution (thread)
				}
			else:
				self.polarizer.propDict = {
					'name': 'polarizer',
					'motor': POLARIZER,
					'res': '8', #Corresponds to mapping above
					'pos': 0, #Position in counts
					'realPos': 0, #Position in inches
					'sPRev': MR_CODE_TO_STEPS_PER_REV['8'] #Steps per revolution (thread)
				}

		for motor in [self.stage, self.polarizer]:
			if motor:
				# Check to make sure the device is in receive mode and reset if necessary
				msg = motor.writeread('RS\r') #RS = Request Status
				motor.flushInput()
				print msg
				if (msg == 'RS=R'):
					print("%s in receive mode." % (motor.propDict['name']))
				elif (msg != 'RS=R'):
					print("%s not in receive mode.  Resetting." % (motor.propDict['name']))
					print "Message was: ",msg
					self.killAllCommands(motor.propDict['motor'])
					if (msg == 'RS=AR'):
						amsg = motor.writeread('AL\r') #AL = Alarm Code
						print 'Alarm message is: ',amsg
						print "Alarm was found. Resetting."
						motor.write('AR\r') #AR = Alarm Reset
						motor.flushInput()
					else:
						print('Irregular message received.')
						sys.exit(1)
				#axis.write('DL2\r')    # Define the limit switch as normally closed
				#axis.write('DL1\r')    # Define the limit switch as normally open (Unplug cables, not in use)

				# Check the microstep resolution
				if mRes:
					motor.write('EG8000\r') #EG = Electronic Gearing
					motor.write('SA\r') #SA = Save Parameters
					motor.flushInput()
					sleep(0.1)
					msg = motor.writeread('EG\r')
					motor.flushInput()
					if(len(msg) <= 4):    # Need at least MR=X + \r, which is 5 characters
						print("Couldn't get microstep resolution for %s.  Assuming 8." % (motor.propDict['name']))    # keeps params from initialization?
					else:
						print msg
						msInfo = msg.rstrip('\r')[3:]
						#print msInfo
						motor.propDict['sPRev'] = float(msInfo)
				else:
					#axis.write('MR\r')    # Request microstep resolution
					msg = motor.writeread('MR\r') #MR = Microstep Resolution
					motor.flushInput()
					#msInfo = axis.read(6)    # Wait to get back microstep info
					if(len(msg) <= 3):    # Need at least MR=X + \r, which is 5 characters
						print("Couldn't get microstep resolution for %s.  Assuming 8." % (motor.propDict['name']))    # keeps params from initialization?
					else:
						msInfo = msg.rstrip('\r')[3:]
						print msInfo
						motor.propDict['res'] = msInfo
						motor.propDict['sPRev'] = MR_CODE_TO_STEPS_PER_REV[msInfo]

			# Set up the limit switches (as normally closed) and check the operating current for the stage motor
			if (motor == self.stage) and (motor != None):
				msg = motor.writeread('DL\r') #DL = Define Limits
				if msg != 'DL=2':
					print "Limits not defined as normally open. Resetting..."
					motor.write('DL2\r')
					sleep(0.1)
					motor.flushInput()
				msg = motor.writeread('CC\r') #CC = Change Current
				print msg
				current = float(msg[3:])
				if current < 1.5:
					print "Operating current insufficient. Resetting..."
					motor.write('CC1.5\r')
			else:
				if motor != None:
					motor.write('JE\r') #JE = Jog Enable


	
	def genMotorList(self, motor):
		"""Get a list of the motors in an FTSControl object.

		Parameters:
		motor (int/motor name) -- STAGE, POLARIZER, or ALL

		Returns:
		mList (list) -- list of desired motors
		"""
		mList = []
		if motor==STAGE or motor==ALL: mList.append(self.stage)
		if motor==POLARIZER or motor==ALL: mList.append(self.polarizer)
		return mList
	

	def startJogging(self, motor=ALL):
		"""Starts jogging control for specified motors."""
		mList = self.genMotorList(motor)
		for motor in mList:
			if not motor:
				print ("Specified motor is invalid - not starting jogging.")
				continue
			motor.write('JE\r') #JE = Jog Enable
			motor.write('WI4L\r') #WI = Wait for Input - Set into wait mode on empty input pin
			motor.flushInput()
			

	def stopJogging(self, motor=ALL):
		"""Stop jogging control."""
		self.killAllCommands(motor)


	def seekHomeLinearStage(self, motor=STAGE):
		"""Move the linear stage to its home position (using the home limit switch)."""
		mList = self.genMotorList(motor)
		for motor in mList:
			if not motor:
				print ("Specified motor is invalid - no motion.")
				continue
			motor.write('VE2.0\r') #VE = Velocity
			motor.write('AC2.0\r') #AC = Acceleration Rate
			motor.write('DE2.0\r') #DE = Deceleration
			motor.write('DI-1\r') #DI = Distance/Position (sets direction)
			motor.write('SHX3L\r') #SH = Seek Home
			motor.flushInput()
			print("FTS linear stage homing...")
			self.blockWhileMoving(motor.propDict['motor'],verbose=True)
		print("FTS linear stage home found.")


#	def seekHomeRotator(self, axis=AXIS_X):
#		aList = self.genAxisList(axis)
#		for axis in aList:
#			axis.write('SHX3H\r')
#			axis.flushInput()
#			print("FTS rotator homing...")
#			self.blockWhileMoving(axis.propDict['axis'])
#		print("FTS rotator home found.")

	
	def setZero(self, motor=ALL):
		"""Tell the motor to set the current position as the zero point."""
		mList = self.genMotorList(motor)
		# for motor in mList:
		# 	motor.propDict['pos'] = 0
		# 	motor.propDict['realPos'] = 0.0
		for motor in mList:
			if not motor:
				print ("Specified motor is invalid.")
				continue
			motor.propDict['pos'] = 0
			motor.propDict['realPos'] = 0.0
			motor.write('SP0\r') #SP = Set Position
			motor.flushInput()


	def getPosition(self, motor=ALL):
		"""Get the position of the motor in counts, relative to the set zero point (or starting point).

		Parameters:
		motor (int/motor name) -- STAGE, POLARIZER, or ALL (default ALL)

		Returns:
		positions (list) -- the positions in counts of the specified motors
		"""
		mList = self.genMotorList(motor)
		positions = []
		for motor in mList:
			if not motor:
				print ("Specified motor is invalid - no position info.")
				positions.append(None)
			else:
				positions.append(motor.propDict['pos'])
		return positions


	def getPositionInInches():
		"""Get the position of the motor in inches, relative to the set zero point (or starting point).

		Parameters:
		motor (int/motor name) -- STAGE, POLARIZER, or ALL (default ALL)

		Returns:
		realPositions (list) -- the positions in inches of the specified motors
		"""
		mList = self.genMotorList(motor)
		realPositions = []
		for motor in mList:
			if not motor:
				print ("Specified motor is invalid - no position info.")
				realPositions.append(None)
			else:
				realPositions.append(motor.propDict['pos'])
		return realPositions

	
	def moveAxisToPosition(self, motor=STAGE, pos=0, posIsInches=False, linStage=True):
		"""Move the axis to the given absolute position in counts or inches.
		
		Parameters:
		motor (int/motor name) - STAGE, POLARIZER, or ALL (default STAGE)
		pos (float) - the desired position in counts or in inches, positive indicates away from the motor (default 0)
		posIsInches (bool) - True if pos was specified in inches, False if in counts (default False)
		linStage (bool) - True if the specified motor is for the linear stage, False if not (default True)
		"""
		mList = self.genMotorList(motor)
		for motor in mList:
			if not motor:
				print ("Specified motor is invalid - no motion.")
				continue
			# Set the threads per inch based on if the motor controls the FTS linear stage
			if linStage:
				AXIS_THREADS_PER_INCH = AXIS_THREADS_PER_INCH_STAGE
			else:
				AXIS_THREADS_PER_INCH = AXIS_THREADS_PER_INCH_XYZ

			# Convert from inches if necessary
			if(posIsInches): unitPos = int(pos*AXIS_THREADS_PER_INCH*motor.propDict['sPRev']/2.0) #2.0 is because threads go twice the distance for one revolution
			else: unitPos = int(pos)

			# Set the new pos and realPos parameters of the motor object
			motor.propDict['pos'] = unitPos
			motor.propDict['realPos'] = 2.0*unitPos/(AXIS_THREADS_PER_INCH*motor.propDict['sPRev']) #See 2.0 note above

			# Move the motor
			motor.write('DI%i\r' % (unitPos)) #DI = Distance/Position
			motor.write('FP\r') #FL = Feed to Position
			motor.flushInput()

			
	def moveAxisByLength(self, motor=STAGE, pos=0, posIsInches=False, linStage=True):
		"""Move the axis relative to the current position by the specified number of counts or inches.

		Parameters:
		motor (int/motor name) - STAGE, POLARIZER, or ALL (default STAGE)
		pos (float) - the desired number of counts or inches to move from current position, positive indicates away from the motor (default 0)
		posIsInches (bool) - True if pos was specified in inches, False if in counts (default False)
		linStage (bool) - True if the specified motor is for the linear stage, False if not (default True)
		"""
		mList = self.genMotorList(motor)
		for motor in mList:
			if not motor:
				print ("Specified motor is invalid - no motion.")
				continue
			# Set the threads per inch based on if the motor controls the FTS linear stage
			if linStage:
				AXIS_THREADS_PER_INCH = AXIS_THREADS_PER_INCH_STAGE
			else:
				AXIS_THREADS_PER_INCH = AXIS_THREADS_PER_INCH_XYZ

			# Convert from inches if necessary
			if(posIsInches): unitPos = int(pos*AXIS_THREADS_PER_INCH*motor.propDict['sPRev']/2.0) #See 2.0 note above
			else: unitPos = int(pos)

			# Set the new pos and realPos parameters of the motor object
			motor.propDict['pos'] += unitPos
			motor.propDict['realPos'] += 2.0*unitPos/(AXIS_THREADS_PER_INCH*motor.propDict['sPRev']) #See 2.0 note above

			# Move the motor
			motor.write('DI%i\r' % (unitPos)) #DI = Distance/Position
			motor.write('FL\r') #FL = Feed to Length
			motor.flushInput()
			print("Final position: ",pos)

			
	def setVelocity(self, motor=ALL, velocity=1.0):
		"""Set velocity in revolutions/second.  Range is .025 - 50.  Accepts floating point values."""
		mList = self.genMotorList(motor)
		for motor in mList:
			if not motor:
				print ("Specified motor is invalid - no velocity set.")
				continue
			motor.write('VE%1.3f\r' % (velocity)) #VE = Velocity
			motor.flushInput()
	

	def setAcceleration(self, motor=ALL, accel=5):
		"""Set acceleration in revolutions/second/second.  Range is 1-3000.  Accepts only integer values."""
		mList = self.genMotorList(motor)
		for motor in mList:
			if not motor:
				print ("Specified motor is invalid - no acceleration set.")
				continue
			motor.write('AC%i\r' % (accel)) #AC = Acceleration Rate
			motor.write('DE%i\r' % (accel)) #DE = Deceleration
			motor.flushInput()

	
	def killAllCommands(self, motor=ALL):
		"""Stop all active commands on the device."""
		mList = self.genMotorList(motor)
		for motor in mList:
			if not motor:
				print ("Specified motor is invalid - no motion.")
				continue
			motor.write('SK\r') #SK = Stop & Kill - Stop/kill all commands, turn off waiting for input
			motor.flushInput()

	
	def blockWhileMoving(self, motor=ALL, updatePeriod=.1, verbose=False):
		"""Block until the specified axes have stop moving.  Checks each axis every updatePeriod seconds.
		
		Parameters:
		motor (int/motor name) -- STAGE, POLARIZER, or ALL (default ALL)
		updatePeriod (float) -- time after which to check each motor in seconds (default .1)
		verbose (bool) -- prints output from motor requests if True (default False)
		"""
		mList = self.genMotorList(motor)
		count = 0
		while(mList):
			count += 1
			for motor in mList:
				motor.flushInput()

				# Get the status of the motor and print if verbose = True
				msg = motor.writeread('RS\r') #RS = Request Status
				motor.flushInput()
				if verbose:
					print msg
					sys.stdout.flush()
				# Remove the motor from mList (so that the while loop continues) only if the status is not "Ready"
				if(msg == 'RS=R'):
					mList.remove(motor) #Should only be in the list once
			# Break if too many while loop iterations - indicates potential problem
			if count > 2000:
				print 'Motion taking too long, there may be a different failure or alarm...'
				break

			# Wait the specified amount of time before rechecking the status
			sleep(updatePeriod)
		print ''

	
	def setMotorEnable(self, motor=ALL, enable=True):
		"""Set motor enable to true or false for given axis.  Should disable motor when stopped for lower noise data acquisition."""
		mList = self.genMotorList(motor)
		for motor in mList:
			if not motor:
				print ("Specified motor is invalid - cannot enable.")
				continue
			if enable==True: motor.write('ME\r') #ME = Motor Enable
			else: motor.write('MD\r') #MD = Motor Disable
			motor.flushInput()


	def retrieveEncoderInfo(self, motor=ALL):
		"""Retrieve encoder step count to verify movement."""
		mList = self.genMotorList(motor)
		for motor in mList:
			if not motor:
				print ("Specified motor is invalid - no encoder info.")
				continue
			ePos = motor.writeread('EP\r') #EP = Encoder Position
			sleep(0.1)
			motor.flushInput()
			print ePos
		return int(ePos.rstrip('\r')[3:])


	def setEncoderValue(self, motor=ALL, value=0):
		"""Set the encoder value in order to keep track of absolute position"""
		mList = self.genMotorList(motor)
		for motor in mList:
			if not motor:
				print ("Specified motor is invalid - encoder value not set.")
				continue
			# Set the encoder position
			ePosSet = motor.writeread('EP%i\r'%(value)) #EP = Encoder Position
			sleep(0.1)
			motor.flushInput()
			# Read and return the new encoder position
			ePos = motor.writeread('EP\r') #EP = Encoder Position
			sleep(0.1)
			motor.flushInput()
			print ePos
		return int(ePos.rstrip('\r')[3:])


	def startRotation(self, motor=POLARIZER, velocity=12.0, accel=1.0):
		"""Starts jogging specifically for the rotation of the output polarizer in the FTS.
		
		Parameters:
		motor (int/motor name) -- desired motor (default POLARIZER)
		velocity (float) -- the rotation velocity in revolutions/second (default 12.0)
		accel (float) -- the acceleration in revolutions/second/second (default 1.0)
		"""
		mList = self.genMotorList(motor)
		for motor in mList:
			if not motor:
				print ("Specified motor is invalid - not starting jogging.")
				continue

			# Set the jog parameters
			motor.write('JS%1.3f\r' % (velocity)) #JS = Jog Speed
			motor.write('JA%i\r' % (accel)) #JA = Jog Acceleration
			motor.write('JL%i\r' % (accel)) #JL = Jog Decel

			# Start rotation
			motor.write('CJ\r') #CJ = Commence Jogging
			motor.flushInput()


	def stopRotation(self, motor=POLARIZER):
		"""Stops jogging specifically for the rotation of the output polarizer in the FTS."""
		mList = self.genMotorList(motor)
		for motor in mList:
			if not motor:
				print ("Specified motor is invalid.")
				continue
			motor.write('SJ\r') #SJ = Stop Jogging
			motor.flushInput()


	def closeConnection(self, motor=ALL):
		"""Close the connection to the serial controller for the specified motor."""
		mList = self.genMotorList(motor)
		for motor in mList:
			if not motor:
				print ("Specified motor is invalid - no connection to close.")
				continue
			motor.sock.close()
		print("Connection to serial controller disconnected.")



class PowerControl(object):
	def __init__(self, PowerIP=None, PowerPort=None):
		self.Power = None
		if not (PowerIP and PowerPort):
			print("Invalid power connection information.  No power control.")
			self.PowerIP = None
			self.PowerPort = None
			self.coninfo = False
		else:
			self.PowerIP = PowerIP
			self.PowerPort = PowerPort
			self.coninfo = True

	def setPower(self, enable=False):
		if not self.coninfo:
			print("Invalid power connection information.  No power control.") 
		elif self.coninfo and not self.Power:
			if enable:
				self.Power = Serial_TCPServer((self.PowerIP, self.PowerPort))
			else:
				print("Power already disabled.")
		elif self.coninfo and self.Power:
			if not enable:
				self.Power.sock.close()
				self.Power = None
			else:
				print("Power already enabled.")
