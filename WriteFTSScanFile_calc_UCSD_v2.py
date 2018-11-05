####################################################################################################
#
# Write FTS Scan File Calculation for UCSD - Version 2
#
# Lindsay Ng Lowry
# llowry@ucsd.edu
# 20181017
#
# This program, like the previous versions, creates a .txt file containing a list of linear stage
# positions at which the FTS should stop and collect data.  The positions are calculated based on
# a frequency resolution and range as well as a microstep resolution which are given as inputs.
# The file includes a header which indicates the resolution, range, microstep resolution, # steps
# per inch, # of positions, and default scan time.
#
# This program differs from previous versions in that it is written such that the motor-defined zero
# point (generally the starting point of the scan) is ALWAYS included as one of the positions at
# which to collect data.
#
####################################################################################################



import numpy, pylab
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('-fn',dest='fn',action='store',type=str,default='',
					help='Suffix for to FTSScanPoints for saved file - include underscore')
parser.add_argument('-int',dest='inttime',action='store',type=float,default=3.0,
					help='Integration time per step in seconds')
parser.add_argument('-pause',dest='pausetime',action='store',type=float,default=0.5,
					help='Time in seconds to wait between moving the stage and recording the data')
parser.add_argument('-res',dest='res',action='store',type=float,default=1.,help='Resolution in GHz')
parser.add_argument('-range',dest='maxrange',action='store',type=float,default=500.,
					help='Maximum range in GHz')
parser.add_argument('-mres',dest='mres',action='store',type=int,default=10,
					help='Microstep resolution of the motor, default set for UCSD FTS')
args = parser.parse_args()



# Get the filename
filename = 'MapFiles/FTSScanPoints' + args.fn + '.txt'

# Get the integration time and the pause time
intTime = args.inttime    #seconds
pauseTime = args.pausetime    #seconds

# Get the resolution and the range
res = args.res    #GHz
maxRange = args.maxrange    #GHz

# Get the microstep resolution and convert to steps per inch
mres = args.mres
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
steps_per_rev = MR_CODE_TO_STEPS_PER_REV[str(mres)]
AXIS_THREADS_PER_INCH_STAGE = 10.0

steps_per_inch = AXIS_THREADS_PER_INCH_STAGE*steps_per_rev/2.0



#print intTime
#print pauseTime
#print res
#print maxRange


# Calculate the extreme points and step size of the roof mirror from the res and range
# See FTS Primer by Hans for this equation

# \delta\nu = c/(2*\pi*d)
#	\delta\nu = resolution [Hz]
#	c = speed of light [m/s]
#	d = max distance roof mirror is away from zero point [m]
d = 3e8/2./(res*1e9)

print 'Max distance [m]: ' + str(d)

# 2*\nu_max = c/(2*s)
#	\nu_max = max frequency from which we get info [Hz]
#	c = speed of light [m/s]
#	s = step size of roof mirror [m]
s = 3e8/4./(maxRange*1e9)

print 'Step Size [m]: ' + str(s)

# Convert the extreme points and step size from meters to encoder counts
min = -1*numpy.int(numpy.round(d/0.0254*steps_per_inch))
max = numpy.int(numpy.round(d/0.0254*steps_per_inch))
step = numpy.int(numpy.round(s/0.0254*steps_per_inch))

print 'Min, Max, Step [encoder counts]: %s, %s, %s'%(min,max,step)


# Create a list of the positions by creating a list of the positive positions starting at zero and
# then negating it in a new list and combining the two
x_pos = numpy.round(numpy.arange(0,max,step))
x_neg = -1*x_pos[1:][::-1]
x = numpy.concatenate((x_neg,x_pos))
print x


# Print the total number of steps and the total time
print('Number of Positions: ' + str(len(x)) + ' points')
scanTime = numpy.float(len(x))*(intTime+pauseTime)/60.0
#print scanTime
print('Time for Scan: ' + str(scanTime) + ' minutes')
scanTime_default = numpy.float(len(x))*(3.0+0.5)/60.0
print('Default Time for Scan: ' + str(scanTime_default) + ' minutes')

# Write the header for the generated txt file
header = 'FTS Scan Points File'
header += '\n\tResolution: ' + str(res) + ' GHz'
header += '\n\tRange: ' + str(maxRange) + ' GHz'
header += '\n\tMicrostep Resolution: ' + str(mres) + ', Steps per inch: ' + str(steps_per_inch)
header += '\n\tNumber of Positions: ' + str(len(x))
header += '\n\tDefault Scan Time: ' + str(scanTime) + ' minutes'



# Save the file
numpy.savetxt(filename,(x), header = header)
