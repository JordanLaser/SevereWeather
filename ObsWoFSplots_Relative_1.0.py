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

	diffU = Stats.variables['DU']
	diffV = Stats.variables['DV']
	diffT = Stats.variables['DT'] 
	diffq = Stats.variables['Dq']
	hWoFS = Stats.variables['h']
	tWoFS = Stats.variables['tWoFS']
	deciU = Stats.variables['deciU'][0,:]
	deciV = Stats.variables['deciV'][0,:]
	decit = Stats.variables['decit'][0,:]
	deciq = Stats.variables['deciq'][0,:]


	###### Average by ensemble #######
	diffUM = np.average(diffU, axis = 0)
	diffVM = np.average(diffV, axis = 0)
	diffTM = np.average(diffT, axis = 0)
	diffqM = np.average(diffq, axis = 0)
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

	diffUysu = (diffU[0,:] + diffU[4,:] + diffU[5,:] + diffU[10,:] + diffU[15,:] + diffU[16,:])/6
	diffUmyj = (diffU[1,:] + diffU[6,:] + diffU[7,:] + diffU[11,:] + diffU[12,:] + diffU[17,:])/6
	diffUmynn = (diffU[2,:] + diffU[3,:] + diffU[8,:] + diffU[9,:] + diffU[13,:] + diffU[14,:])/6

	diffVysu = (diffV[0,:] + diffV[4,:] + diffV[5,:] + diffV[10,:] + diffV[15,:] + diffV[16,:])/6
	diffVmyj = (diffV[1,:] + diffV[6,:] + diffV[7,:] + diffV[11,:] + diffV[12,:] + diffV[17,:])/6
	diffVmynn = (diffV[2,:] + diffV[3,:] + diffV[8,:] + diffV[9,:] + diffV[13,:] + diffV[14,:])/6

	diffTysu = (diffT[0,:] + diffT[4,:] + diffT[5,:] + diffT[10,:] + diffT[15,:] + diffT[16,:])/6
	diffTmyj = (diffT[1,:] + diffT[6,:] + diffT[7,:] + diffT[11,:] + diffT[12,:] + diffT[17,:])/6
	diffTmynn = (diffT[2,:] + diffT[3,:] + diffT[8,:] + diffT[9,:] + diffT[13,:] + diffT[14,:])/6

	diffqysu = (diffq[0,:] + diffq[4,:] + diffq[5,:] + diffq[10,:] + diffq[15,:] + diffq[16,:])/6
	diffqmyj = (diffq[1,:] + diffq[6,:] + diffq[7,:] + diffq[11,:] + diffq[12,:] + diffq[17,:])/6
	diffqmynn = (diffq[2,:] + diffq[3,:] + diffq[8,:] + diffq[9,:] + diffq[13,:] + diffq[14,:])/6

	tWoFSysu = (tWoFS[0,:] + tWoFS[4,:] + tWoFS[5,:] + tWoFS[10,:] + tWoFS[15,:] + tWoFS[16,:])/6
	tWoFSmyj = (tWoFS[1,:] + tWoFS[6,:] + tWoFS[7,:] + tWoFS[11,:] + tWoFS[12,:] + tWoFS[17,:])/6
	tWoFSmynn = (tWoFS[2,:] + tWoFS[3,:] + tWoFS[8,:] + tWoFS[9,:] + tWoFS[13,:] + tWoFS[14,:])/6


	####### Average all levels ######
	# diffUMlvl = round(np.average(diffUM, axis = 0),2)
	# diffVMlvl = round(np.average(diffVM, axis = 0),2)
	# diffTMlvl = round(np.average(diffTM, axis = 0),2)
	# diffqMlvl = round(np.average(diffqM, axis = 0),2)

	##### Find PBL height with 0.6+ K rule
	PBLind = 0
	theta0 = np.average(decit[:3])
	for i in range(len(decit)-1):
		if (theta0 + 0.6) < decit[i+1]:
			PBLind = i + 1
			break
	ysuind = 0
	theta0ysu = np.average(tWoFSysu[:3])
	for i in range(len(decit)-1):
		if (theta0ysu + 0.6) < tWoFSysu[i+1]:
			ysuind = i + 1
			break
	myjind = 0
	theta0myj = np.average(tWoFSmyj[:3])
	for i in range(len(decit)-1):
		if (theta0myj + 0.6) < tWoFSmyj[i+1]:
			myjind = i + 1
			break
	mynnind = 0
	theta0mynn = np.average(tWoFSmynn[:3])
	for i in range(len(decit)-1):
		if (theta0mynn + 0.6) < tWoFSmynn[i+1]:
			mynnind = i + 1
			break


# 	#Wind Differences plot
	vertical_levels = hWoFSM[:]
	heightticks = [0,400,800,1200,1600,2000]
	LS = 20
	ls = 18
	a = 0.8
	lw = 0.7
	fig = pylab.figure(figsize = (16,13))  
	fig.suptitle('Radiosonde Observation May 17th ' + time + ' UTC')
	ax1 = fig.add_subplot(1, 2, 1)
	ax1.set_title('Zonal Wind', fontsize = LS)
	plt.xlabel('WoFS U - Obs U (m/s)', fontsize = ls)
	ax1.set_yticks(heightticks)
	plt.ylabel('Altitude (m)', fontsize = ls, rotation = 'horizontal')
	ax1.axvline(x = 0, color = 'k')
	ax1.axhline(y = hWoFSM[PBLind], color = 'k')
	ax1.axhline(y = hWoFSM[ysuind], color = 'r')
	ax1.axhline(y = hWoFSM[myjind], color = 'b')
	ax1.axhline(y = hWoFSM[mynnind], color = 'm')
	ax1.set_xlim(-10,10)
	ax1.set_ylim(-100,2000)
	ax2 = fig.add_subplot(1 ,2 ,2)
	ax2.set_title('Meridional Wind', fontsize = LS)
	ax2.axvline(x = 0, color = 'k')
	ax2.axhline(y = hWoFSM[PBLind], color = 'k')
	plt.xlabel('WoFS V - Obs V (m/s)', fontsize = ls)
	ax2.set_yticks(heightticks)
	ax2.axhline(y = hWoFSM[PBLind], color = 'k')
	ax2.axhline(y = hWoFSM[ysuind], color = 'r')
	ax2.axhline(y = hWoFSM[myjind], color = 'b')
	ax2.axhline(y = hWoFSM[mynnind], color = 'm')
	ax2.set_xlim(-10,10)
	ax2.set_ylim(-100,2000)
	#Plotting each member color coded by pbl scheme
	for i in ysu:
		ax1.plot(diffU[i,:],vertical_levels, linewidth = lw, alpha = a, color = 'r')
		ax2.plot(diffV[i,:],vertical_levels, linewidth = lw, alpha = a, color = 'r')
	for i in myj:
		ax1.plot(diffU[i,:],vertical_levels, linewidth = lw, alpha = a, color = 'b')
		ax2.plot(diffV[i,:],vertical_levels, linewidth = lw, alpha = a, color = 'b')
	for i in mynn:
		ax1.plot(diffU[i,:],vertical_levels, linewidth = lw, alpha = a, color = 'm')
		ax2.plot(diffV[i,:],vertical_levels, linewidth = lw, alpha = a, color = 'm')
	#plotting the mean U by pbl scheme
	p1, = ax1.plot(diffUysu[:],vertical_levels, linewidth = 4.0, label = 'YSU', color = 'r')
	p1, = ax1.plot(diffUmyj[:],vertical_levels, linewidth = 4.0, label = 'MYJ', color = 'b')
	p1, = ax1.plot(diffUmynn[:],vertical_levels, linewidth = 4.0, label = 'MYNN', color = 'm')
 	#plotting the mean V by pbl scheme
	p1, = ax2.plot(diffVysu[:],vertical_levels, linewidth = 4.0, label = 'YSU', color = 'r')
	p1, = ax2.plot(diffVmyj[:],vertical_levels, linewidth = 4.0, label = 'MYJ', color = 'b')
	p1, = ax2.plot(diffVmynn[:],vertical_levels, linewidth = 4.0, label = 'MYNN', color = 'm')
 	#plotting the ensemble average
	p1, = ax1.plot(diffUM[:],vertical_levels, color = 'k', linewidth = 5.0, label = 'WoFS Mean')
	p1, = ax2.plot(diffVM[:],vertical_levels, color = 'k', linewidth = 5.0)
	
	box = ax1.get_position()
	ax1.set_position([box.x0, box.y0, box.width, box.height])
	legend = ax2.legend(loc='lower left', bbox_to_anchor=(1,0), shadow=True)
	plt.savefig('plots/Wind_Plots/windPBLdiffplt_' + str(date) + '_' + str(time) + '.png')
# 	
	
	# 	Thermo Differences plot
	fig = pylab.figure(figsize = (16,13))  
	fig.suptitle('Observation ' + time + 'UTC')
	ax1 = fig.add_subplot(1, 2, 1)
	ax1.set_title('WoFS \u03B8 - Obs \u03B8 versus vertical level')
	plt.xlabel('Difference from Observation (K), positive difference implies warm bias in model')
	ax1.set_yticks(heightticks)
	plt.ylabel('Vertical level in WoFS')
	ax1.axvline(x = 0, color= 'k')
	ax1.axhline(y = hWoFSM[PBLind], color = 'k')
	ax1.axhline(y = hWoFSM[PBLind], color = 'k')
	ax1.axhline(y = hWoFSM[ysuind], color = 'r')
	ax1.axhline(y = hWoFSM[myjind], color = 'b')
	ax1.axhline(y = hWoFSM[mynnind], color = 'm')
	ax1.set_xlim(-10,10)
	ax1.set_ylim(-100,2000)
	ax2 = fig.add_subplot(1 ,2 ,2)
	ax2.set_title('WoFS qv - Obs qv versus vertical level')
	ax2.axvline(x = 0, color = 'k')
	ax2.axhline(y = hWoFSM[PBLind], color = 'k')
	ax2.set_yticks(heightticks)
	ax2.set_xlim(-10,10)
	ax2.axhline(y = hWoFSM[PBLind], color = 'k')
	ax2.axhline(y = hWoFSM[ysuind], color = 'r')
	ax2.axhline(y = hWoFSM[myjind], color = 'b')
	ax2.axhline(y = hWoFSM[mynnind], color = 'm')
	plt.xlabel('Difference from Observation (Kg/Kg), positive difference implies positive mixing ratio bias in model')
	ax2.set_xlim(-0.01,0.01)
	ax2.set_ylim(-100,2000)
	#Plotting each member color coded by pbl scheme
	for i in ysu:
		ax1.plot(diffT[i,:],vertical_levels, linewidth = lw, alpha = a, color = 'r')
		ax2.plot(diffq[i,:],vertical_levels, linewidth = lw, alpha = a, color = 'r')
	for i in myj:
		ax1.plot(diffT[i,:],vertical_levels, linewidth = lw, alpha = a, color = 'b')
		ax2.plot(diffq[i,:],vertical_levels, linewidth = lw, alpha = a, color = 'b')
	for i in mynn:
		ax1.plot(diffT[i,:],vertical_levels, linewidth = lw, alpha = a, color = 'm')
		ax2.plot(diffq[i,:],vertical_levels, linewidth = lw, alpha = a, color = 'm')
	#plotting the mean U by pbl schemes
	p1, = ax1.plot(diffTysu[:],vertical_levels, linewidth = 4.0, label = 'YSU', color = 'r')
	p1, = ax1.plot(diffTmyj[:],vertical_levels, linewidth = 4.0, label = 'MYJ', color = 'b')
	p1, = ax1.plot(diffTmynn[:],vertical_levels, linewidth = 4.0, label = 'MYNN', color = 'm')
 	#plotting the mean V by pbl scheme
	p1, = ax2.plot(diffqysu[:],vertical_levels, linewidth = 4.0, label = 'YSU', color = 'r')
	p1, = ax2.plot(diffqmyj[:],vertical_levels, linewidth = 4.0, label = 'MYJ', color = 'b')
	p1, = ax2.plot(diffqmynn[:],vertical_levels, linewidth = 4.0, label = 'MYNN', color = 'm')
 	#plotting the ensemble average
	p1, = ax1.plot(diffTM[:],vertical_levels, color = 'k', linewidth = 5.0, label = 'WoFS Mean')
	p1, = ax2.plot(diffqM[:],vertical_levels, color = 'k', linewidth = 5.0)
	
	box = ax1.get_position()
	ax1.set_position([box.x0, box.y0, box.width, box.height])
	legend = ax2.legend(loc='lower left', bbox_to_anchor=(1,0), shadow=True)

	plt.savefig('plots/Thermo_Plots/thermoPBLdiffplt_' + str(date) + '_' + str(time) + '.png')

	plt.rcParams.update({'figure.max_open_warning': 0})
	
	
