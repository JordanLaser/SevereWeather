import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import os
import pylab
import sys
from netCDF4 import Dataset
from optparse import OptionParser

##################################################
#Input directory containing observation statistics
##################################################
parser = OptionParser()
parser.add_option("-o", dest="indir", type="string", default= None, help="Input directory")

(options, args) = parser.parse_args()

if (options.indir == None):
   parser.print_help()
   print()
   sys.exit(1)
else:
   indir = options.indir

################################################
#For loop generating plots#####################
###############################################
files = os.listdir(indir)
files.sort()
print('Here are all the files in the folder you designated: ' + str(files))

for f, infile in enumerate(files):
	infile = os.path.join(indir, infile)
	Stats = Dataset(infile)
	date = infile[-16:-8]
	time = infile[-7:-3]
	print("Date of Observation: " + str(date))
	print("Time of Observation: " + str(time))
	print("Current file: " + str(infile))

	##########Defining variables###############

	diffU0 = Stats.variables['DU0']
	diffU1 = Stats.variables['DU1']
	diffU2 = Stats.variables['DU2']
	diffU3 = Stats.variables['DU3']
	diffV0 = Stats.variables['DV0']
	diffV1 = Stats.variables['DV1']
	diffV2 = Stats.variables['DV2']
	diffV3 = Stats.variables['DV3']
	hWoFS = Stats.variables['h']
#	deciU = Stats.variables['deciU'][0,:]
#	deciV = Stats.variables['deciV'][0,:]



	###### Average by ensemble #######
	diffUM0 = np.average(diffU0, axis = 0)
	diffUM1 = np.average(diffU1, axis = 0)
	diffUM2 = np.average(diffU2, axis = 0)
	diffUM3 = np.average(diffU3, axis = 0)
	
	diffVM0 = np.average(diffV0, axis = 0)
	diffVM1 = np.average(diffV1, axis = 0)
	diffVM2 = np.average(diffV2, axis = 0)
	diffVM3 = np.average(diffV3, axis = 0)
	
	hWoFSM = np.average(hWoFS, axis = 0)

	#####PBL indicies#######
	member_list = [1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 2, 3, 4, 5, 6, 7, 8, 9]
	#ysu = [1,2,7,8,13,14]
	#myj = [3,4,9,10,15,16]            #pbl schemes based off of member
	#mynn = [5,6,11,12,17,18]
	ysu = [0,4,5,10,15,16]
	myj = [1,6,7,11,12,17]             #pbl indicies
	mynn = [2,3,8,9,13,14]

	###### Average by PBL scheme ######

# 	diffUysu = (diffU[0,:] + diffU[4,:] + diffU[5,:] + diffU[10,:] + diffU[15,:] + diffU[16,:])/6
# 	diffUmyj = (diffU[1,:] + diffU[6,:] + diffU[7,:] + diffU[11,:] + diffU[12,:] + diffU[17,:])/6
# 	diffUmynn = (diffU[2,:] + diffU[3,:] + diffU[8,:] + diffU[9,:] + diffU[13,:] + diffU[14,:])/6
# 
# 	diffVysu = (diffV[0,:] + diffV[4,:] + diffV[5,:] + diffV[10,:] + diffV[15,:] + diffV[16,:])/6
# 	diffVmyj = (diffV[1,:] + diffV[6,:] + diffV[7,:] + diffV[11,:] + diffV[12,:] + diffV[17,:])/6
# 	diffVmynn = (diffV[2,:] + diffV[3,:] + diffV[8,:] + diffV[9,:] + diffV[13,:] + diffV[14,:])/6


	####### Average all levels ######
	# diffUMlvl = round(np.average(diffUM, axis = 0),2)
	# diffVMlvl = round(np.average(diffVM, axis = 0),2)
	# diffTMlvl = round(np.average(diffTM, axis = 0),2)
	# diffqMlvl = round(np.average(diffqM, axis = 0),2)


# 	#Wind Differences plot
 	vertical_levels = hWoFSM[:]
	heightticks = [0,400,800,1200,1600,2000]
	LS = 20
	ls = 18
	a = 0.8
	lw = 0.7
	fig = pylab.figure(figsize = (16,13))  
	ax1 = fig.add_subplot(1, 2, 1)
	ax1.set_title('Zonal Wind', fontsize = LS)
	plt.xlabel('WoFS U - Obs U (m/s)', fontsize = ls)
	ax1.set_yticks(heightticks)
	plt.ylabel('Alt (m)', fontsize = ls, rotation = 'horizontal')
	ax1.axvline(x = 0, color = 'k')
	ax1.set_xlim(-10,10)
	ax1.set_ylim(0,1200)
	ax2 = fig.add_subplot(1 ,2 ,2)
	ax2.set_title('Meridional Wind', fontsize = LS)
	ax2.axvline(x = 0, color = 'k')
	plt.xlabel('WoFS V - Obs V (m/s)', fontsize = ls)
	ax2.set_yticks(heightticks)
	ax2.set_xlim(-10,10)
	ax2.set_ylim(0,1200)

	p1, = ax1.plot(diffUM0[:],vertical_levels, color = 'b', label = 'Lidar 2230', linewidth = 5.0, alpha = 0.4)
	p1, = ax2.plot(diffVM0[:],vertical_levels, color = 'b', linewidth = 5.0, alpha = 0.4)
	
	p1, = ax1.plot(diffUM1[:],vertical_levels, color = 'b', label = 'Lidar 2245',linewidth = 5.0, alpha = 0.6)
	p1, = ax2.plot(diffVM1[:],vertical_levels, color = 'b', linewidth = 5.0, alpha = 0.6)
	
	p1, = ax1.plot(diffUM2[:],vertical_levels, color = 'b', label = 'Lidar 2300',linewidth = 5.0, alpha = 0.8)
	p1, = ax2.plot(diffVM2[:],vertical_levels, color = 'b', linewidth = 5.0, alpha = 0.8)
	
	p1, = ax1.plot(diffUM3[:],vertical_levels, color = 'b', label = 'Lidar 2315',linewidth = 5.0, alpha = 1.0)
	p1, = ax2.plot(diffVM3[:],vertical_levels, color = 'b', linewidth = 5.0, alpha = 1.0)
	
	slegend = ax2.legend(loc='lower left', shadow=True)
	plt.savefig('plots/Wind_Plots/Lidar_winddiffplt_' + str(date) + '_' + str(time) + '.png', dpi = 300)
	plt.show()
	
	plt.rcParams.update({'figure.max_open_warning': 0})
	
	
