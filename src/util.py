import os
import math
import shutil
import ROOT as r
import numpy as np
import src.constants as constants

# compute the cosine Fourier transform
def calc_cosine_transform(t0, binContent, binCenter, freq_step_size, n_freq_step, lower_freq):

    a = []
    b = []

    for i in range(0, n_freq_step):
        frequency   = ( lower_freq/1000 + freq_step_size/1000/2) + i*freq_step_size/1000 # in MHz
        a.append( frequency*1000 ) # back to kHz
        b.append( np.sum( binContent*np.cos(2*math.pi*frequency*(binCenter-t0))*0.001 ) )# time hist in micro-sec, dt is 1 ns

    return a, b


# compute the sine Fourier transform
def calc_sine_transform(t0, binContent, binCenter):

    a = []
    b = []

    for i in range(0, constants.nFreq):
        frequency   = ( constants.lower_freq/1000 + constants.freq_step/1000/2) + i*constants.freq_step/1000 # in MHz
        a.append( frequency*1000 ) # back to kHz
        b.append( np.sum( binContent*np.sin(2*math.pi*frequency*(binCenter-t0))*0.001 ) )# time hist in micro-sec, dt is 1 ns

    return a, b

# convert Frequency to Radius ==#
# assume velocity corresponding to magic momentum for all the muons. Velocity depends very little on momentum ==#
def convert_freq_to_radius( freqHist, radius, intensity, n_freq_step ):
    for i in range( 1, n_freq_step+1 ):
        radius.append   ( constants.speedOfLight * constants.magicBeta / ( 2*math.pi*freqHist.GetBinCenter(i) ) )
        intensity.append( freqHist.GetBinContent(i))


# compute the parabola correction distribution
def calc_parabola(t0, tS, firstApprox, parabola):
    for i in range(0, constants.nFreq):
        frq = ( constants.lower_freq + constants.freq_step/2) + i*constants.freq_step # in MHz
        integral = 0
        for j in range(1, constants.nFreq+1):
            if (firstApprox.GetBinCenter(j)-frq) != 0:
                integral += firstApprox.GetBinContent(j)*np.sin(2*math.pi*(frq-firstApprox.GetBinCenter(j))*(tS-t0)/1000)/(1000*(frq-firstApprox.GetBinCenter(j))  )
            else:
                integral += 2*math.pi*firstApprox.GetBinContent(j)*(tS-t0)/1000000
        parabola.SetBinContent(i+1,integral)


# perform the background minimization
def minimization(parabola, cosine):

    x = []
    y = []

    for i in range(1, constants.nFreq+1):
        if ( cosine.GetBinCenter(i) < constants.lowerCollimatorFreq or cosine.GetBinCenter(i) > constants.upperCollimatorFreq ):
            x.append(parabola.GetBinContent(i))
            y.append(cosine.GetBinContent(i))
            x.append(parabola.GetBinContent(i))
            y.append(cosine.GetBinContent(i))

    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A, y,rcond=-1)[0]

    return m, c


# manage results directory
def manage_results_directory(dir_path):

    print(' ### Step 1/4: manage results directory\n')

    if os.path.isdir(dir_path) and not os.listdir(dir_path):
        print('    output directoy already exists and is empty --> will save results there\n')
        os .mkdir(dir_path + '/t0_optimization')
    elif os.path.isdir(dir_path) and os.listdir(dir_path):
        print('    output directoy already exists and is NOT empty --> deleting/recreating directory (previous results will be lost)\n')
        shutil.rmtree(dir_path)
        os.mkdir(dir_path)
        dir_path = dir_path + '/t0_optimization'
        os.mkdir(dir_path)
    else:
        print('    output directory does not exist --> creating it\n')
        os      .mkdir(dir_path)
        dir_path = dir_path + '/t0_optimization'
        os.mkdir(dir_path)

# print to terminal the welcome message
def print_welcome_message():

    print('')
    print(' ------------------------------')
    print(' |                            |')
    print(' |           CORNELL          |')
    print(' |        FAST ROTATION       |')
    print(' |      FOURIER ANALYSIS      |')
    print(' |                            |')
    print(' | contact: atc93@cornell.edu |')
    print(' |                            |')
    print(' ------------------------------')
    print('')

#== Convert ROOT histogram to Numpy array ==#
def rootHistToNumpArray( hist, tS, tM ):

    startBin = hist.FindBin(tS)
    endBin   = hist.FindBin(tM)

    binCenter   = np.empty( int(endBin-startBin+1), dtype=float )
    binContent  = np.empty( int(endBin-startBin+1), dtype=float )

    for j in range(startBin, endBin+1):
            binContent[j-startBin]  = hist.GetBinContent(j)
            binCenter[j-startBin]   = hist.GetBinCenter(j)

    return binCenter, binContent

#== Convert from ring global radial coordinate to beam local radial coordinate ==#
def globalToLocalRadialCoordinate( graph ):

    nPoint = graph.GetN()

    for i in range( 0, nPoint ):

        x, y = r.Double(), r.Double()
        graph.GetPoint( i, x, y )
        graph.SetPoint( i, x-constants.magicR, y )

#== Compute Radial Mean within collimator aperture ==#
def computeRadialMean( radius, intensity ):

    mean = 0
    sumI = 0

    for x,y in zip( radius, intensity ):

        #== Discard data point if radius outside of collimator aperture ==#
        if ( x < constants.lowerCollimatorRad or x > constants.upperCollimatorRad ):
            continue

        #== Else compute the Mean ==#
        mean += x*y
        sumI += y

    mean /= sumI

    return mean

#== Compute Radial Standard Deviation within collimator aperture ==#
def computeRadialSTD( radius, intensity, meanRad ):

    std  = 0
    sumI = 0

    for x,y in zip( radius, intensity ):
        
        #== Discard data point if radius outside of collimator aperture ==#
        if ( x < constants.lowerCollimatorRad or x > constants.upperCollimatorRad ):
            continue

        #== Else compute the STD ==#
        sumI += y
        std  += y * ( x-meanRad ) * ( x-meanRad )

    std /= sumI
    std = math.sqrt(std) 

    return std

#== Compute the E-field correction ==#

def computeEfieldCorrection( n, mean, std ):

    return ( - 2* math.pow( constants.magicBeta, 2 ) * n * ( 1-n ) * ( math.pow ( mean, 2 ) + math.pow ( std, 2 ) ) / ( math.pow ( constants.magicR, 2 ) )* 1e9 )

#== Extract the two minima of the Cosine Fourier transform ==#
def extractMinima( hist ):
    hist.GetXaxis().SetRangeUser( constants.lower_freq, constants.magicFreq )
    min1 = hist.GetMinimum()
    minBinIdx1 = hist.GetMinimumBin()
    hist.GetXaxis().SetRangeUser( constants.magicFreq, constants.upper_freq )
    min2 = hist.GetMinimum()
    minBinIdx2 = hist.GetMinimumBin()
    hist.GetXaxis().SetRangeUser( constants.lower_freq, constants.upper_freq )
    return min1, min2, abs( min1-min2 ), minBinIdx1, minBinIdx2


import numba as nb
@nb.njit(fastmath=True,parallel=True)
def compute_numba(t0, a, b, it):
    res=np.empty((it,a.shape[0]))
    #res=np.empty((a.shape[0]))
    #res=np.empty((it))
    ita=np.arange(0,it)

    #ntegral=np.zero(0, it).astype(np.float32)
    integral = np.zeros((it))
    #print('size: ', integral.size)
    for i in nb.prange(ita.shape[0]):
    #for i in nb.prange(ita.shape[0]):
    #for i in nb.prange(1):
    #for i in nb.prange(ita.shape[0]):
        #cpt = i
        #print(i)
        frequency   = ( constants.lower_freq/1000 + constants.freq_step/1000/2) + i*constants.freq_step/1000 # in MHz
        #print(i, frequency)
        #t=ita[i]
        #print(cpt, res[cpt])
        lala = 0.
        for j in range(a.shape[0]):
            #res[i,j]=a[j] * np.cos( 2. * np.pi * it * b[j])
            #res[j] = a[j] * np.cos( 2. * np.pi * it * b[j])
            #integral[i] += a[j] * np.cos( 2. * np.pi * frequency * (b[j]-t0)) * 0.001
            lala += a[j] * np.cos( 2. * np.pi * frequency * (b[j]-t0)) * 0.001
            #res[j] = a[j] * np.cos( 2. * np.pi * frequency * (b[j]-t0)) * 0.001
            #print(integral[i])
            #print( a[j] * np.cos( 2. * np.pi * frequency * (b[j]-t0)) * 0.001 )
        #print(lala)
        #print(i, integral[i])
        #print(np.sum(res))
        integral[i] = lala
        #integral[i] = np.sum(res[i])
        #print(integral[i])
        #print(sum)
            #np.sum(res)
            #np.sum( a[j] * np.cos( 2. * np.pi * frequency * (b[j]-t0) )*0.001 )

    return integral
