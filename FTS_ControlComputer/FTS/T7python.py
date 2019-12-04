#update 08/24/2018
from labjack import ljm
import time
import sys
import numpy as np
from datetime import datetime

class T7():
	def __init__(self):
		try:
			self.handle = ljm.openS("T7", "ANY", "ANY")
			info = ljm.getHandleInfo(self.handle)
			print ("Opened a LabJack with Device type: %i, Connection type: %i,\n" \
                    "Serial number: %i, IP address: %s, Port: %i,\nMax bytes per MB: %i" % \
                    (info[0], info[1], info[2], ljm.numberToIP(info[3]), info[4], info[5]))
            
		except:
			print("Error opening Labjack")
			self.handle = None
        	
	def __del__(self):
        #Destructor. Ensures handle is closed properly.

		from labjack import ljm
		if self.handle:
			ljm.close(self.handle)
		
	def stream_txtrecord(self,maxRequest,scanRate,aScanListNames,aNames,aValues):

		tstr = time.strftime('%Y%m%d_%H%M%S',time.gmtime(time.time()))
		savefile = tstr+'_timestream.txt'
		
		f=open(savefile,'wb')

		NUM_IN_CHANNELS = len(aScanListNames)
		scansPerRead = int(scanRate/2)
		handle = self.handle
		aScanList = ljm.namesToAddresses(NUM_IN_CHANNELS, aScanListNames)[0]
    	# Ensure triggered stream is disabled.
		ljm.eWriteName(handle, "STREAM_TRIGGER_INDEX", 0)
    	 # Enabling internally-clocked stream.
		ljm.eWriteName(handle, "STREAM_CLOCK_SOURCE", 0)

		# stream settling time and stream resolution configuration.
		numFrames = len(aNames)
		ljm.eWriteNames(handle, numFrames, aNames, aValues)
		scanRate = ljm.eStreamStart(handle, scansPerRead, NUM_IN_CHANNELS, aScanList, scanRate)
		print("\nStream started with a scan rate of %0.0f Hz." % scanRate)
		print("\nPerforming %i stream reads." % maxRequest)
    
		start = datetime.now()
		totScans = 0
		totSkip = 0  # Total skipped samples
		i = 1
    	
		streamdata=[]
		try:
			while i <= maxRequest or maxRequest == -1:
    
				ret = ljm.eStreamRead(handle)
				data = ret[0][0:(scansPerRead * NUM_IN_CHANNELS)]
				scans = len(data) / NUM_IN_CHANNELS
				totScans += scans
        		# Count the skipped samples which are indicated by -9999 values. Missed
        		# samples occur after a device's stream buffer overflows and are
        		# reported after auto-recover mode ends.
				curSkip = data.count(-9999.0)
				totSkip += curSkip
				for j in range(0, scansPerRead):
					recordStr = " "
					for k in range(0, NUM_IN_CHANNELS):
						recordStr += str(data[j * NUM_IN_CHANNELS + k])+','
					f.write(recordStr+'\n')
				i += 1
        		
					
			end = datetime.now()
			f.close()
			print("\nTotal scans = %i" % (totScans))
			tt = (end - start).seconds + float((end - start).microseconds) / 1000000
			print("Time taken = %f seconds" % (tt))
			print("LJM Scan Rate = %f scans/second" % (scanRate))
			print("Timed Scan Rate = %f scans/second" % (totScans / tt))
			print("Timed Sample Rate = %f samples/second" % (totScans * NUM_IN_CHANNELS/tt))
			print("Skipped scans = %0.0f" % (totSkip / NUM_IN_CHANNELS))
		except ljm.LJMError:
			ljme = sys.exc_info()[1]
			print(ljme)
		except Exception:
			e = sys.exc_info()[1]
			print(e)
		try:
			print("\nStop Stream")
			

			ljm.eStreamStop(handle)
		except ljm.LJMError:
			ljme = sys.exc_info()[1]
			print(ljme)
		except Exception:
			e = sys.exc_info()[1]
			print(e)
    		
		ljm.close(handle)
		
		return np.array(streamdata)
            
	def stream_array(self,maxRequest,scanRate,aScanListNames,aNames,aValues):

		
		NUM_IN_CHANNELS = len(aScanListNames)
		scansPerRead = int(scanRate/2)
		handle = self.handle
		aScanList = ljm.namesToAddresses(NUM_IN_CHANNELS, aScanListNames)[0]
    	# Ensure triggered stream is disabled.
		ljm.eWriteName(handle, "STREAM_TRIGGER_INDEX", 0)
    	 # Enabling internally-clocked stream.
		ljm.eWriteName(handle, "STREAM_CLOCK_SOURCE", 0)

		# stream settling time and stream resolution configuration.
		numFrames = len(aNames)
		ljm.eWriteNames(handle, numFrames, aNames, aValues)
		scanRate = ljm.eStreamStart(handle, scansPerRead, NUM_IN_CHANNELS, aScanList, scanRate)
#		print("\nStream started with a scan rate of %0.0f Hz." % scanRate)
#		print("\nPerforming %i stream reads." % maxRequest)
    
		start = datetime.now()
		totScans = 0
		totSkip = 0  # Total skipped samples
		i = 1
    	
		streamdata=[]
		try:
			while i <= maxRequest or maxRequest == -1:
    
				ret = ljm.eStreamRead(handle)
				data = ret[0][0:(scansPerRead * NUM_IN_CHANNELS)]
				scans = len(data) / NUM_IN_CHANNELS
				totScans += scans
        		# Count the skipped samples which are indicated by -9999 values. Missed
        		# samples occur after a device's stream buffer overflows and are
        		# reported after auto-recover mode ends.
				curSkip = data.count(-9999.0)
				totSkip += curSkip
				for j in range(0, scansPerRead):
					recordStr = []
					for k in range(0, NUM_IN_CHANNELS):
						recordStr.append(data[j * NUM_IN_CHANNELS + k])
						streamdata.append(recordStr)
				i += 1
		except ljm.LJMError:
			ljme = sys.exc_info()[1]
			print(ljme)
		except Exception:
			e = sys.exc_info()[1]
			print(e)
			
		return np.array(streamdata)
	def stopstream(self):
		handle = self.handle
		ljm.eStreamStop(handle)
#		ljm.close(handle)
		
	def forceclose(self):
		ljm.close(self.handle)
        	
        	
        	
        	
        	
    	