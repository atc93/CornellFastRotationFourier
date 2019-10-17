#Josh Fagin
#Monte Carlo Fast Rotation Signal

from math import exp,erfc,sqrt,pi
import time
from random import uniform,randint
from multiprocessing import Pool
import numpy as np
import ROOT as r
from array import array
import os
from shutil import copyfile,rmtree
import json
import sys
from glob import glob1

#open config file
try:
	with open(sys.argv[1]) as config_file:
		data = json.load(config_file)
except:
	print("Config file \""+str(sys.argv[1])+"\" not found.\nYou must use the correct path to a config json file.") 
	sys.exit(1) #aborts

#welcome message :)
print('')
print(' ------------------------------')
print(' |                            |')
print(' |           CORNELL          |')
print(' |        FAST ROTATION       |')
print(' |         MONTE CARLO        |')
print(' |                            |')
print(' | contact:jsf243@cornell.edu |')
print(' |                            |')
print(' ------------------------------')
print('')

##get config data##

#can use numba

try:
	use_numba = data["use_numba"]
	if use_numba:
		try:
			from numba import jit
		except:
			print("Could not import numba! Not using numba.")
			use_numba = False
	numfrs = data["num_frs"]
	startnum = data["start_num"]
	numProcess = data["num_processors"]
	multi = data["multiprocess"]
	tm = data["tm"]
	binwidth = data["binwidth"]
	name = data["name"]
	fmin = data["fmin"] #[MHz]
	fmax = data["fmax"]

	use_directory = data["put_in_directory"]
	directory = data["name_directory"]

	#frequency
	frequency_random_num_gauss = data["frequency_random_num_gauss"]
	frequency_num_gauss_lower = data["frequency_num_gauss_lower"]
	frequency_num_gauss_upper = data["frequency_num_gauss_upper"]
	frequency_num_gauss = data["frequency_num_gauss"]

	frequency_random_center = data["frequency_random_center"]        
	frequency_center_lower = data["frequency_center_lower"]
	frequency_center_upper = data["frequency_center_upper"]
	frequency_center = data["frequency_center"]

	frequency_random_width = data["frequency_random_width"]
	frequency_width_lower = data["frequency_width_lower"]
	frequency_width_upper = data["frequency_width_upper"]
	frequency_width = data["frequency_width"]

	frequency_random_skew = data["frequency_random_skew"]
	frequency_skew_lower = data["frequency_skew_lower"]
	frequency_skew_upper = data["frequency_skew_upper"]
	frequency_skew = data["frequency_skew"]

	frequency_random_height = data["frequency_random_height"]
	frequency_height_lower = data["frequency_height_lower"]
	frequency_height_upper = data["frequency_height_upper"]
	frequency_height = data["frequency_height"]

	#longitude
	longitude_random_num_gauss = data["longitude_random_num_gauss"]
	longitude_num_gauss_lower = data["longitude_num_gauss_lower"]
	longitude_num_gauss_upper = data["longitude_num_gauss_upper"]
	longitude_num_gauss = data["longitude_num_gauss"]

	longitude_random_center = data["longitude_random_center"]
	longitude_center_lower = data["longitude_center_lower"]
	longitude_center_upper = data["longitude_center_upper"]        
	longitude_center = data["longitude_center"]

	longitude_random_width = data["longitude_random_width"]
	longitude_width_lower = data["longitude_width_lower"]
	longitude_width_upper = data["longitude_width_upper"]
	longitude_width = data["longitude_width"]

	longitude_random_skew = data["longitude_random_skew"]
	longitude_skew_lower = data["longitude_skew_lower"]
	longitude_skew_upper = data["longitude_skew_upper"]
	longitude_skew = data["longitude_skew"]

	longitude_random_height = data["longitude_random_height"]
	longitude_height_lower = data["longitude_height_lower"]
	longitude_height_upper = data["longitude_height_upper"]
	longitude_height = data["longitude_height"]

	height_multigauss_factor_freq = data["height_multigauss_factor_freq"] 
	height_multigauss_factor_long = data["height_multigauss_factor_long"] 

	random_t0 = data["random_t0"]
	t0_range_lower = data["t0_range_lower"]
	t0_range_upper = data["t0_range_upper"]
	t0_val = data["t0_val"]
	real_t0_optimization = data["real_t0_optimization"]
	ts = data["ts_for_optimization"] #micro seconds
	create_Fourier_config = data["create_Fourier_config"]
	Fourier_config_name = data["Fourier_config_name"]

	realistic = data["add_noise"]
	random_signal_to_noise_ratio = data["random_signal_to_noise_ratio"]
	signal_to_noise_ratio_lower = data["signal_to_noise_ratio_lower"]
	signal_to_noise_ratio_upper = data["signal_to_noise_ratio_upper"]
	signal_to_noise_ratio = data["signal_to_noise_ratio"]
	noise_type = data["noise_type"]
	exponential_increasing_noise = data["exponential_increasing_noise"]
	
	add_momentum_time_correlation = data["add_momentum_time_correlation"]
	random_momentum_time_correlation = data["random_momentum_time_correlation"]
	momentum_time_correlation = data["momentum_time_correlation"]
	momentum_time_correlation_lower = data["momentum_time_correlation_lower"]
	momentum_time_correlation_upper = data["momentum_time_correlation_upper"]
	
	momentum_time_correlation_scan = data["momentum_time_correlation_scan"]
	momentum_time_scan_dr = data["momentum_time_scan_dr"]
	
	print_graphs = data["print_graphs"]
	print_time = data["print_time"]
	display_counter = data["display_counter"]
	save_result_text = data["save_result_text"]
	n = data["field_index"]

	n_t0_opt = data["n_t0_opt"]

	background_fit = data["background_fit"]
	remove_background = data["remove_background"]
	background_frequencies = data["background_frequencies"]
	background_removal_threshold = data["background_removal_threshold"]
	t0_background_threshold = data["t0_background_threshold"]
	tS = data["tS"]
	freq_step_size = data["freq_step_size"]
except:
	print("Error! Missing parameters in config file. Exiting.")
	sys.exit(1) #aborts

if momentum_time_correlation_scan:
	numfrs = int((momentum_time_correlation_upper-momentum_time_correlation_lower)/momentum_time_scan_dr)+1

#defines important constants
T = 0.149126 
f0 = 1.0/T   
vacuumMin = 6.66279793407 #[MHz]
vacuumMax = 6.74765032038
decayRate = 64.46 #micro seconds
r.gStyle.SetOptStat(0)
r.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")

#CE calculation
muonMass = 0.105658
magicP = 3.094
c = 299792458.
E = sqrt(muonMass*muonMass+magicP*magicP)
gamma = E / muonMass
beta = magicP/(gamma*muonMass)
speed = beta*c
magicradius = 7.112 #[m]

#assertion statements#
assert tm > ts or not real_t0_optimization,"Cannot optimize t0 unless tm > "+str(ts)+" micro seconds"
assert noise_type == "gauss" or noise_type == "uniform","Parameter noise_type must be \"gauss\" or \"uniform\" but got:"+"\""+str(noise_type)+"\""
#makes sure the boolean parameters are boolean
boolList = [real_t0_optimization,use_directory,frequency_random_num_gauss,frequency_random_center,
		frequency_random_width,frequency_random_skew,frequency_random_height,longitude_random_num_gauss,longitude_random_center,
		longitude_random_width,longitude_random_skew,longitude_random_height,realistic,random_signal_to_noise_ratio,
		exponential_increasing_noise,print_time,display_counter,print_graphs,create_Fourier_config,remove_background,use_numba]
for boolPar in boolList:
	assert type(boolPar) == bool,"\""+str(boolPar)+"\""+" is not a boolean!"
#makes sure the integer parameters are integers
intList = [numfrs,startnum,numProcess,frequency_num_gauss_lower,frequency_num_gauss_upper,frequency_num_gauss,
		longitude_num_gauss_lower,longitude_num_gauss_upper,longitude_num_gauss,n_t0_opt]
for intPar in intList:
	assert type(intPar) == int and intPar >= 0,"\""+str(intPar)+"\""+" is not a positive integer!"

#Function to calculates the E-field correction
def getCE(f,flist):
    radius,amp = [],[]
    for i in range(len(f)):
        binCenter = flist[i]*1e6
        distance = speed/(2*pi*binCenter)
        if speed/(2*pi*vacuumMax*1e6) <= distance <= speed/(2*pi*vacuumMin*1e6):
            radius.append(distance)
            amp.append(f[i])
    var,mean,add = 0.0,0.0,0.0
    for i,j in zip(radius,amp):
        add += j
        var += j*(i-magicradius)*(i-magicradius)
        mean+= i*j
    var /= add
    mean/= add
    std = 0.0
    for i,j in zip(radius,amp):
        std += j*(i-mean)*(i-mean)
    std = np.sqrt(std/add)
    CE = -2*beta*beta*n*(1-n)*var/(magicradius*magicradius)*1e9
    return CE,mean,std
#functions to optimize the real value of t0

if use_numba:
	@jit(nopython = True)
	def cosineFFT(t0,bincontent,bincenter):  
		binnum = int((fmax-fmin)/binwidth) 
		frequency = []
		for i in range(binnum):
		    f = fmin + i*binwidth
		    integral = bincontent*np.cos(2*pi*f*(bincenter-t0))
		    frequency.append(np.sum(integral))
		return frequency
else:
	def cosineFFT(t0,bincontent,bincenter):  
		binnum = int((fmax-fmin)/binwidth) 
		frequency = []
		for i in range(binnum):
		    f = fmin + i*binwidth
		    integral = bincontent*np.cos(2*pi*f*(bincenter-t0))
		    frequency.append(np.sum(integral))
		return frequency		
		
def grapht0(numt0,mint0,stept0,fr):
    t0list = []
    for i in range(numt0):
        t0 = mint0 + stept0*i
        t0list.append(t0)
    freqList = []
    a,b = [],[]
    for i in range(int(ts/binwidth),int(tm/binwidth)):
    	a.append(fr[i])
    	b.append(i*binwidth)
    bincontent = np.asarray(a)
    bincenter = np.asarray(b)

    for t0 in t0list:
        freqList.append(cosineFFT(t0,bincontent,bincenter))
    return freqList,t0list

def findt0(freqList):
    sym = []
    for freq in freqList:
        magicBin = np.argmax(freq)
        if magicBin == 0 or magicBin == len(freq):
        	magicBin = int((f0-fmin)/binwidth)
        sym.append(abs(min(freq[0:magicBin])-min(freq[magicBin:])))
    index = sym.index(min(sym))
    return index,sym

     
#where the collimator appurture is
deltac = T*vacuumMax-1.0

#defines our distributions
def function(gausslist,numgauss,height_multigauss_factor):
	height,center,std,decay = gausslist[0]
	if numgauss == 1:
		if use_numba:
			@jit(nopython = True)
			def fun(x):
				return height*exp(-(x-center)**2/(2*std*std))*erfc(decay*(x-center)/std-1.0)
		else:
			def fun(x):
				return height*exp(-(x-center)**2/(2*std*std))*erfc(decay*(x-center)/std-1.0)
	elif numgauss == 2:
		height2,center2,std2,decay2 = gausslist[1]
		if use_numba:
			@jit(nopython = True)
			def fun(x):
				return height*exp(-(x-center)**2/(2*std*std))*erfc(decay*(x-center)/std-1.0)+(
					height2/height_multigauss_factor)*exp(-(x-center2-center)**2/(2*std2*std2))*erfc(decay2*(x-center2-center)/std2-1.0)
		else:
			def fun(x):
				return height*exp(-(x-center)**2/(2*std*std))*erfc(decay*(x-center)/std-1.0)+(
					height2/height_multigauss_factor)*exp(-(x-center2-center)**2/(2*std2*std2))*erfc(decay2*(x-center2-center)/std2-1.0)
	elif numgauss == 3:
		height2,center2,std2,decay2 = gausslist[1]
		height3,center3,std3,decay3 = gausslist[2]
		if use_numba:
			@jit(nopython = True)
			def fun(x):
				return height*exp(-(x-center)**2/(2*std*std))*erfc(decay*(x-center)/std-1.0)+(
					height2/height_multigauss_factor)*exp(-(x-center2)**2/(2*std2*std2))*erfc(decay2*(x-center2)/std2-1.0)+(
					height3/height_multigauss_factor)*exp(-(x-center3)**2/(2*std3*std3))*erfc(decay3*(x-center3)/std3-1.0)
		else:
			def fun(x):
				return height*exp(-(x-center)**2/(2*std*std))*erfc(decay*(x-center)/std-1.0)+(
					height2/height_multigauss_factor)*exp(-(x-center2)**2/(2*std2*std2))*erfc(decay2*(x-center2)/std2-1.0)+(
					height3/height_multigauss_factor)*exp(-(x-center3)**2/(2*std3*std3))*erfc(decay3*(x-center3)/std3-1.0)
	return fun

#computes the fast rotation signal
if use_numba:
	@jit(nopython = True)
	def time_bin(t):
		if t < 10.0: return 0.000010
		elif t < 20.0: return 0.000020
		elif t < 30.0: return 0.000025
		elif t < 40.0: return 0.000050
		elif t < 50.0: return 0.000200
		elif t < 70.0: return 0.000250
		else: return 0.000500
else:
	def time_bin(t):
		if t < 10.0: return 0.000010
		elif t < 20.0: return 0.000020
		elif t < 30.0: return 0.000025
		elif t < 40.0: return 0.000050
		elif t < 50.0: return 0.000200
		elif t < 70.0: return 0.000250
		else: return 0.000500

#Deletes graph files already made if they are not in a directory
if not use_directory:
	total_num = len(glob1("","*.root"))
	num_removed = 0
	for dist in range(startnum,startnum+numfrs):
		plot_dir = "Dist_"+str(dist+1)
		if os.path.exists(plot_dir):
			print("Directory \""+plot_dir+"\" already exists. Replacing.")
			rmtree(plot_dir)
			os.remove(name+str(dist+1)+".root")
			num_removed+=1
	num_extra_dist = total_num-num_removed		
	
def frs(distribution):
	#t0 is chosen
	t0 = uniform(t0_range_lower,t0_range_upper) if random_t0 else t0_val
	
	#initializes a frequency distribution
	
	numgauss = randint(frequency_num_gauss_lower,frequency_num_gauss_upper) if frequency_random_num_gauss else frequency_num_gauss
	gausslistFreq = []
	for j in range(numgauss):
		stdrho = uniform(frequency_width_lower,frequency_width_upper) if frequency_random_width else frequency_width 
		scew = uniform(frequency_skew_lower,frequency_skew_upper) if frequency_random_skew else frequency_skew
		center = uniform(frequency_center_lower,frequency_center_upper) if frequency_random_center else frequency_center
		height = uniform(frequency_height_lower,frequency_height_upper) if frequency_random_height else frequency_height
		gauss = (height,center,stdrho,scew)
		gausslistFreq.append(gauss)
	
	rho = function(gausslistFreq,numgauss,height_multigauss_factor_freq)
	
	#initializes a random longitudinal beam profile
	numgauss = randint(longitude_num_gauss_lower,longitude_num_gauss_upper) if longitude_random_num_gauss else longitude_num_gauss
	gausslistLong = []
	for j in range(numgauss):
		stdxi = uniform(longitude_width_lower,longitude_width_upper) if longitude_random_width else longitude_width 
		scew = uniform(longitude_skew_lower,longitude_skew_upper) if longitude_random_skew else longitude_skew
		center = uniform(longitude_center_lower,longitude_center_upper) if longitude_random_center else longitude_center
		height = uniform(longitude_height_lower,longitude_height_upper) if longitude_random_height else longitude_height
		gauss = (height,center,stdxi,scew)
		gausslistLong.append(gauss)
	xi = function(gausslistLong,numgauss,height_multigauss_factor_long)
	
	#sets the momentum time correlation
	if add_momentum_time_correlation:
		#for momentum time scan
		if momentum_time_correlation_scan:
			correlation = momentum_time_correlation_lower+distribution*momentum_time_scan_dr
		else:
			correlation = uniform(momentum_time_correlation_lower,momentum_time_correlation_upper) if random_momentum_time_correlation else momentum_time_correlation 
	else:
		correlation = 0.0
		
	#gets the mean of the beam profile distribution for use in the correlation
	axis = np.linspace(-3*T,3*T,num = 5000) #the longitudinal beam profile should be fully confined within 3 period
	weights = []
	for tprime in axis:
		weights.append(xi(tprime))
	mean_beam_profile = np.average(axis,weights = weights)
	
	#defines the function to outputs frs at a given time	
	if use_numba:
		@jit(nopython = True)
		def S(t):
			content = 0.0
			begin = int(t*vacuumMin)
			if begin <= 10:
				begin=0
			else:
				begin-=8
			dt = time_bin(t)
			for n in range(begin,int(t*vacuumMax)+8):
				periods = n*T+t0
				for j in range(int((t-(1+deltac)*periods)/dt),int((t-(1-deltac)*periods)/dt)):
					tprime = j*dt
					content += dt*xi(tprime)*rho((t-tprime)/(periods)-1.0+(correlation*(tprime-mean_beam_profile)/T))/(periods)
			return content	

	else:
		def S(t):
			content = 0.0
			begin = int(t*vacuumMin)
			if begin <= 10:
				begin=0
			else:
				begin-=8
			dt = time_bin(t)
			for n in range(begin,int(t*vacuumMax)+8):
				periods = n*T+t0
				for j in range(int((t-(1+deltac)*periods)/dt),int((t-(1-deltac)*periods)/dt)):
					tprime = j*dt
					content += dt*xi(tprime)*rho((t-tprime)/(periods)-1.0+(correlation*(tprime-mean_beam_profile)/T))/(periods)
			return content	

	f = r.TH1D('f', 'f', int((fmax-fmin)/binwidth), fmin, fmax)
	frs = r.TH1D('frs', 'frs', int(tm/binwidth), 0.0, tm)
	
	#creates the frs
	numbins = int(tm/binwidth)
	tlist = []
	intensity = []
	for i in range(numbins):
		t = frs.GetBinCenter(i+1)
		tlist.append(t)
		intensity.append(S(t))
	
	#normalizes the fast rotation signal
	intensity = np.asarray(intensity)
	norm = numbins/np.sum(intensity)
	intensity = norm*intensity
	
	#adds if you want a more realistic distribution
	if realistic:
		#finds the noise amplitude
		max_after = 4.0 if tm > 4.0 else 0.0
		signalAmp = np.max(list(intensity)[int(max_after/binwidth):])-1.0 #max value after 4 micros and subtracts mean value of 1
		if random_signal_to_noise_ratio:
			signalToNoiseRation = uniform(signal_to_noise_ratio_lower,signal_to_noise_ratio_upper)
		else:
			signalToNoiseRation = signal_to_noise_ratio
		noiseAmp = signalAmp/signalToNoiseRation
		
		if noise_type == "gauss": #Gaussian Noise
			if exponential_increasing_noise:
				intensity = intensity + np.exp(np.asarray(tlist)/decayRate)*np.random.normal(0.0,noiseAmp/3.0,numbins)
			else:
				intensity = intensity + np.random.normal(0.0,noiseAmp/3.0,numbins)
		elif noise_type == "uniform": #Uniform Noise
			if exponential_increasing_noise:
				intensity = intensity + np.exp(np.asarray(tlist)/decayRate)*np.random.uniform(-noiseAmp,noiseAmp,numbins)
			else:
				intensity = intensity + np.random.uniform(-noiseAmp,noiseAmp,numbins)
		
		#normalizes again		
		norm = numbins/np.sum(intensity)
		intensity = norm*intensity
	
	intensity = list(intensity)
	
	#creates real frequency distribution
	binnum = int((fmax-fmin)/binwidth)
	freq = []
	if not add_momentum_time_correlation:
		for i in range(binnum):
			delta = -1*(f.GetBinCenter(i+1)/f0-1.0)
			if -deltac <= delta <= deltac:
				freq.append(rho(delta))
			else:
				freq.append(0.0)
	else:		
		correlation_freq = np.zeros(binnum)
		flist = []
		for i in range(binnum):
			flist.append(f.GetBinCenter(i+1))
		for tprime in axis:
			delta = -1*(np.array(flist)/f0-1.0)-correlation*(mean_beam_profile-tprime)/T
			rho_vals = []
			for delta_val in delta: 
				if -deltac <= delta_val <= deltac: #only include physical values (in the ring)
					rho_vals.append(rho(delta_val))
				else:
					rho_vals.append(rho(delta_val))
			correlation_freq += np.array(rho_vals)*xi(tprime)
		for i,f_val in enumerate(flist):
			if vacuumMin <= f_val <= vacuumMax:
				freq.append(correlation_freq[i+1])
			else:
			 	freq.append(0.0)
			
	longR = r.TH1D('longR', 'longR', 4*binnum,-2.0*binnum*binwidth,2.0*binnum*binwidth)
	#creates longitudinal beam profile
	longitude,longlist = [],[]
	for i in range(4*binnum):
		x = longR.GetBinCenter(i+1)
		longitude.append(xi(x))
		longlist.append(x)
		
	#normalizes real distributions
	norm = max(freq)
	freq = [el/norm for el in freq]
	norm = max(longitude)
	longitude = [el/norm for el in longitude]
	flist = []
	
	#Creates ROOT file :)
	for i,freqVal in enumerate(freq):
		f.SetBinContent(i+1,freqVal) #kHz
		flist.append(f.GetBinCenter(i+1))
	for i,frsVal in enumerate(intensity):
		frs.SetBinContent(i+1,frsVal)
	for i,longVal in enumerate(longitude):
		longR.SetBinContent(i+1,longVal)
	radius = array("d")
	amp = array("d")
	for i in range(binnum):
		binCenter = f.GetBinCenter(i+1)*1e6
		distance = speed/(2*pi*binCenter)
		radius.append(1000*distance) #mm
		amp.append(f.GetBinContent(i+1))
	radiusGraph = r.TGraph(binnum,radius,amp)
	
	if momentum_time_correlation_scan:
		dist_str = "_r_{}".format(round(correlation,6))
	else:
		dist_str = str(distribution+1)	 
	#creates the ROOT file
	if use_directory:
		fr = r.TFile(directory+"/root_files/"+name+dist_str+'.root','RECREATE')
	else:
		fr = r.TFile("root_files/"+name+dist_str+'.root','RECREATE')
	
	f.SetBins(binnum,1000*fmin,1000*fmax) #sets the frequency to kHz
	f.Write('freq')
	frs.Write('frs')
	radiusGraph.Write('radius')
	longR.Write('longitude')
	fr.Close()
	
	#creates plots
	plot_dir = "Dist_"+dist_str 
	if print_graphs or save_result_text:
		os.mkdir(directory+"/Plots_Results/"+plot_dir) if use_directory else os.mkdir("Plots_Results/"+plot_dir)
	path = directory+"/Plots_Results/"+plot_dir+"/" if use_directory else "Plots_Results/"+plot_dir+"/"
	
	if print_graphs:
		c = r.TCanvas( 'c', 'Graph', 200, 10, 700, 500 )
		f.Draw('hist')
		f.SetTitle('Frequency Distribution')
		f.GetXaxis().SetTitle('frequency [MHz]')
		c.Print(path+'Frequency.png')
		
		radiusGraph.Draw('AC')
		radiusGraph.SetTitle('Radial Distribution')
		radiusGraph.GetXaxis().SetTitle('radius [mm]')
		c.Print(path+'Radius.png')
		
		longR.Draw('hist')
		longR.SetTitle('Longitudinal Beam Profile')
		longR.GetXaxis().SetTitle('t [#mu s]')
		c.Print(path+'Longitude.png')
		
		#Draws different ranges of the fast rotation signal
		frs.Draw('hist')
		frs.SetTitle('Fast Rotation Signal')
		frs.GetXaxis().SetTitle('t [#mu s]')
		c.Print(path+'FRSFull.png')
		if tm>=1:
			frs.Draw('hist')
			frs.SetTitle('Fast Rotation Signal')
			frs.GetXaxis().SetTitle('t [#mus]')
			frs.GetXaxis().SetRangeUser(0.0,1.0)
			c.Print(path+'FRSfirst_turn.png')
		if tm>=10:
			frs.Draw('hist')
			frs.SetTitle('Fast Rotation Signal')
			frs.GetXaxis().SetTitle('t [#mus]')
			frs.GetXaxis().SetRangeUser(0.0,10.0)
			c.Print(path+'FRSbeginning.png')
		if tm>=11:
			frs.Draw('hist')
			frs.SetTitle('Fast Rotation Signal')
			frs.GetXaxis().SetTitle('t [#mus]')
			frs.GetXaxis().SetRangeUser(10.0,11.0)
			c.Print(path+'FRStiny.png')
		if tm>=20:
			frs.Draw('hist')
			frs.SetTitle('Fast Rotation Signal')
			frs.GetXaxis().SetTitle('t [#mus]')
			frs.GetXaxis().SetRangeUser(10.0,20.0)
			c.Print(path+'FRSsmall.png')
		if tm>=40:
			frs.Draw('hist')
			frs.SetTitle('Fast Rotation Signal')
			frs.GetXaxis().SetTitle('t [#mus]')
			frs.GetXaxis().SetRangeUser(10.0,40.0)
			c.Print(path+'FRSmid.png')
		if tm>=100:
			frs.Draw('hist')
			frs.SetTitle('Fast Rotation Signal')
			frs.GetXaxis().SetTitle('t [#mus]')
			frs.GetXaxis().SetRangeUser(10.0,100.0)
			c.Print(path+'FRSlarge.png')
		if tm>=200:
			frs.Draw('hist')
			frs.SetTitle('Fast Rotation Signal')
			frs.GetXaxis().SetTitle('t [#mus]')
			frs.GetXaxis().SetRangeUser(0.0,200.0)
			c.Print(path+'FRS200.png')
		if tm>=400:
			frs.Draw('hist')
			frs.SetTitle('Fast Rotation Signal')
			frs.GetXaxis().SetTitle('t [#mus]')
			frs.GetXaxis().SetRangeUser(0.0,400.0)
			c.Print(path+'FRS400.png')
		#Zooms into FRS at large times to see the noise
		if realistic and tm>=50:
			frs.Draw('hist')
			frs.SetTitle('Fast Rotation Signal')
			frs.GetXaxis().SetTitle('t [#mus]')
			frs.GetXaxis().SetRangeUser(48.0,50.0)
			c.Print(path+'Noise_t_50.png')
		if realistic and tm>=100:
			frs.Draw('hist')
			frs.SetTitle('Fast Rotation Signal')
			frs.GetXaxis().SetTitle('t [#mus]')
			frs.GetXaxis().SetRangeUser(98.0,100.0)
			c.Print(path+'Noise_t_100.png')
		if realistic and tm>=200:
			frs.Draw('hist')
			frs.SetTitle('Fast Rotation Signal')
			frs.GetXaxis().SetTitle('t [#mus]')
			frs.GetXaxis().SetRangeUser(198.0,200.0)
			c.Print(path+'Noise_t_200.png')
		if realistic and tm>=300:
			frs.Draw('hist')
			frs.SetTitle('Fast Rotation Signal')
			frs.GetXaxis().SetTitle('t [#mus]')
			frs.GetXaxis().SetRangeUser(298.0,300.0)
			c.Print(path+'Noise_t_300.png')
		if realistic and tm>=400:
			frs.Draw('hist')
			frs.SetTitle('Fast Rotation Signal')
			frs.GetXaxis().SetTitle('t [#mus]')
			frs.GetXaxis().SetRangeUser(398.0,400.0)
			c.Print(path+'Noise_t_400.png')
	
	#optimizes the real value of t0		
	if real_t0_optimization:
		freqCalc = []
		t0Graph,minDifft0 = [],[]
		t0Graph2,minDifft02 = [],[]
		t0Graph3,minDifft03 = [],[]
		iterationFreq,iterationCorrection = [],[]
		frList = [],[]
		
		#look over a full period
		mint0 = 0.000
		maxt0 = 0.150
		numt0 = int((maxt0 - mint0)/binwidth)+1
		freqList,t0list = grapht0(numt0,mint0,0.001,intensity[:])
		index,minDiff = findt0(freqList)
		t0Graph = t0list
		minDifft0 = minDiff
		centert0 = t0list[index]
		#finds t0 to the nearest .1 ns
		freqList,t0list = grapht0(21,centert0-.001,0.0001,intensity[:])
		index,minDiff = findt0(freqList)
		t0Graph2 = t0list
		minDifft02 = minDiff
		centert0 = t0list[index]
		frequency = freqList[index]
		#finds t0 to the nearest .01 ns
		freqList,t0list = grapht0(21,centert0-.0001,0.00001,intensity[:])
		index,minDiff = findt0(freqList)
		t0Graph3 = t0list
		minDifft03 = minDiff
		t0_optimized = t0list[index]
		
		if print_graphs:
			t0_optimization_graph = r.TGraph(len(t0Graph),array("d",t0Graph),array("d",minDifft0))
			t0_optimization_graph2 = r.TGraph(len(t0Graph2),array("d",t0Graph2),array("d",minDifft02))
			t0_optimization_graph3 = r.TGraph(len(t0Graph3),array("d",t0Graph3),array("d",minDifft03))
			
			t0_optimization_graph.Draw('AC')
			t0_optimization_graph.SetTitle('t0 optimization t0 1 ns')
			t0_optimization_graph.GetXaxis().SetTitle('t0 [#mus]')
			c.Print(path+'t0_optimization1.png')
			
			t0_optimization_graph2.Draw('AC')
			t0_optimization_graph2.SetTitle('t0 optimization to 0.1 ns')
			t0_optimization_graph2.GetXaxis().SetTitle('t0 [#mus]')
			c.Print(path+'t0_optimization2.png')
			
			t0_optimization_graph3.Draw('AC')
			t0_optimization_graph3.SetTitle('t0 optimization to 0.01 ns')
			t0_optimization_graph3.GetXaxis().SetTitle('t0 [#mus]')
			c.Print(path+'t0_optimization3.png')
	else:
		t0_optimized = t0
	#creates a text file with all the values used
	if save_result_text:
		with open(path+"values.txt", "w") as f:
			CE,mean,stdev = getCE(freq,flist) #set to MHz
			f.write("Possibly Random Parameters")
			f.write("\n")
			f.write("\n")
			f.write("c_e: "+str(CE)+" ppb, eq_radius: "+str(1000.0*(mean-magicradius))+" mm , std: "+str(1000*stdev)+" mm")
			f.write("\n")
			if real_t0_optimization:
				f.write("t0 optimized to: "+str(t0_optimized)+" micro seconds")
				f.write("\n")
			if realistic:
				f.write("Signal to noise ratio: "+str(signalToNoiseRation))
				f.write("\n")
			if add_momentum_time_correlation:
				f.write("Correlation: "+str(correlation))
				f.write("\n")
			f.write("t0: " + str(t0)+" micro s")
			f.write("\n")
			f.write("\n")
			f.write("Fixed parameters")
			f.write("\n")
			f.write("\n")
			f.write("n: "+str(n)+", Height multigauss factor freq: "+str(height_multigauss_factor_freq)+", Height multigauss factor long: "+str(height_multigauss_factor_long)+", binwidth: "+str(binwidth)+" micro s")
			f.write("\n")
			f.write("Add noise: "+str(realistic)+", Noise type: "+str(noise_type)+", Exponential increasing noise: "+str(exponential_increasing_noise))
			f.write("\n")
			f.write("\n")
			f.write("Frequency and Longitude Parameters")
			f.write("\n")
			f.write("\n")
			f.write("Num Gauss Frequency: "+str(len(gausslistFreq)))
			f.write("\n")
			for i,gauss in enumerate(gausslistFreq):
				f.write("Frequency Gauss: "+str(i+1))
				f.write("\n")
				f.write("Height: " +str(gauss[0])+", Center: " +str(gauss[1])+", Stdev: " +str(gauss[2])+", Skew: " +str(gauss[3]))
				f.write("\n")	
			f.write("\n")
			f.write("Num Gauss Longitude: "+str(len(gausslistLong)))
			f.write("\n")
			for i,gauss in enumerate(gausslistLong):
				f.write("Longitude Gauss "+str(i+1))
				f.write("\n")
				f.write("Height: " +str(gauss[0])+", Center: " +str(gauss[1])+", Stdev: " +str(gauss[2])+", Skew: " +str(gauss[3]))
				f.write("\n")
	
	#creates json files to run Fourier analysis code
	if create_Fourier_config:
		if real_t0_optimization == False:
			print()
			print("Warning t0 was not optimized! Used the fixed value of t0 instead.")
			print()
		path = directory+"/root_files/"+name+dist_str+'.root' if use_directory else "root_files/"+name+dist_str+'.root'
		#config file for fourier analysis
		
			
		jsonOutput = {
		    "tag":"Output_dist"+dist_str,
		    "root_file":path,
			"histo_name":"frs",
		    "background_correction":"fit",
		    "background_fit":background_fit,
		    "background_frequencies":background_frequencies,
		    "remove_background":remove_background,        
		    "background_removal_threshold":background_removal_threshold,
		    "t0_background_threshold":t0_background_threshold,
			"lower_t0":t0_optimized-0.0020, #Looks at +/- 2.0 ns from the optimized t0
			"upper_t0":t0_optimized+0.0020,
			"t0_step_size":0.000025,
			"tS":tS,
			"tM":tm,
		    "n_t0_opt":n_t0_opt,
		    "n_fit_param":0,
		    "freq_step_size":freq_step_size,
		    "lower_freq":1000.0*fmin,
		    "upper_freq":1000.0*fmax,
		    "rebin_wiggle_factor":149,
		    "rebin_frs_factor":1,
		    "field_index":n,
		    "start_fit_time":30,
			"poly_order":4,
        	"cbo_freq":370,
			"print_plot":True,
			"calc_sine":False,
		    "verbose":1,
		    "append_results":False,

		    "compare_with_truth":False,
		    "truth_root_file":"",
			"truth_freq_histo_name":"",
			"truth_rad_histo_name":"",

		    "fix_t0":False,
		    "fixed_t0_value":0.1216,
		    "run_t0_scan":False,
		
		    "check_positron_hits":False,
		    "positron_hits_threshold":100000,
			
			"fix_fit_bound": False,
  	        "fit_lower_bound": 6671,
    	    "fit_upper_bound": 6725,
		    
		    "run_tS_scan":False,
			"lower_tS":4.000,
			"upper_tS":15.000,
		    "tS_step_size":0.1,

		    "run_tM_scan":False,
			"lower_tM":300.000,
			"upper_tM":300.300,
		    "tM_step_size":0.001,

		    "run_freq_step_scan":False,
		    "lower_freq_step_size":0.25,
		    "upper_freq_step_size":3.75,
		    "freq_step_size_increment":0.25,

		    "stat_fluctuation": False,
		    "n_stat_fluctuation": 0,

		    "run_background_threshold_scan":False,
		    "lower_background_threshold":1,
		    "upper_background_threshold":5,
		    "background_threshold_step_size":0.25,

		    "run_background_removal_threshold_scan":False,
		    "lower_background_removal_threshold":1,
		    "upper_background_removal_threshold":5,
		    "background_removal_threshold_step_size":0.25
		}
		path = directory+"/Fourier_config/" if use_directory else "Fourier_config/"
		with open(path+Fourier_config_name+dist_str+'.json', 'w') as json_file:
			json.dump(jsonOutput, json_file)
		
	
	#counter looks at how many ROOT files are made so far. This is a good way to get a counter to work when multiprocessing.
	if display_counter:
		num_files = len(glob1(directory+"/root_files","*.root")) if use_directory else len(glob1("root_files","*.root"))-num_extra_dist
		string = str(round(100*float(num_files)/float(numfrs),2))+ " % complete. " + str(num_files)+"/"+str(numfrs)
		sys.stdout.write('\r'+str(string))
		sys.stdout.flush()

if __name__ == '__main__':
	start = time.time()
	#copies the config file
	if use_directory:
		if os.path.exists(directory):
			print("Directory \"" + directory + "\" already exists. Replacing.")
			rmtree(directory)
		os.mkdir(directory)
		copyfile(sys.argv[1],directory + '/' + 'config.json')
	else:
		copyfile(sys.argv[1],'config.json')
	
	#creates the Fourier config files
	if create_Fourier_config:
		os.mkdir(directory+"/Fourier_config") if use_directory else os.mkdir("Fourier_config")
	#creates the file to put plots and results in
	if print_graphs or save_result_text:
		os.mkdir(directory+"/Plots_Results") if use_directory else os.mkdir("Plots_Results")
	#creates the file to put ROOT files in
	os.mkdir(directory+"/root_files") if use_directory else os.mkdir("root_files")
	print("")
	#counter
	if display_counter: 
		sys.stdout.write('\r'+"0.0 % complete. "+"0/"+str(numfrs))
		sys.stdout.flush()
	#Creates the FRS
	if not momentum_time_correlation_scan:
		if multi: Pool(numProcess).map(frs,range(startnum,startnum+numfrs))
		else: 
			for dist in range(startnum,startnum+numfrs): frs(dist)
	else:
		#for momentum time scan
		if multi: Pool(numProcess).map(frs,range(numfrs))
		else: 
			for dist in range(numfrs): frs(dist)	
	#prints how long it took
	print("\n")
	if print_time:
		dif = time.time()-start
		if dif < 100.0: print('This took ' + str(round(dif,3)) + ' seconds')
		elif  7000.0 >= dif >= 100.0: print('This took ' + str(round(dif/60.0,3)) + ' minutes')
		else: print('This took ' + str(round(dif/3600.0,3)) + ' hours')
		
