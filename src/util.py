from scipy.optimize import curve_fit
import scipy.special as special
import os
import sys
import math
import shutil
import ROOT as r
import numpy as np
import src.constants as constants


def fit_sinc(a, b, err, scale_factor, verbose):

    def func(x, a, b, c, d):
        return a * np.sin(np.asarray(x-c)/d)/((x-c)/d) + b

    try:
        popt, pcov = curve_fit(func, a, b, bounds=(
            [-10*scale_factor, -10, 6695, 5], [0, 0, 6715, 50]), sigma=err, maxfev=1000)
    except:
        print("Sinc background fit failed: return chi2 of -1")
        fit_status = -1
        return -1, -1, -1, -1, -1

    if (verbose > 1):
        print('    bkg fit param: ', popt)

    residuals = b - func(a, *popt)
    chi2 = sum((residuals / err) ** 2)/len(residuals)

    return func, chi2, np.std(residuals), popt, pcov

def fit_triangle(a, b, err, ts, verbose):

    def func(f,amp, time, a, b, c):
        return -amp*(((f-a)/(b-a))*(special.sici(time*2*math.pi*(f-a))[0]-special.sici(time*2*math.pi*(f-b))[0])+(np.cos(time*2*math.pi*(f-a))-np.cos(time*2*math.pi*(f-b)))/(time*2*math.pi*(b-a))+
                 ((c-f)/(b-c))*(special.sici(time*2*math.pi*(f-c))[0]-special.sici(time*2*math.pi*(f-b))[0])-(np.cos(time*2*math.pi*(f-c))-np.cos(time*2*math.pi*(f-b)))/(time*2*math.pi*(b-c)))

    initial_guess = [10.0, ts/1000, constants.magicFreq-10, constants.magicFreq+10, constants.magicFreq] # tS converted to ms because frequency in kHz
    try:
        popt, pcov = curve_fit(func, a, b, p0=initial_guess, maxfev=100000)
    except:
        print("    Triangle background fit failed: return chi2 of -1")
        fit_status = -1
        return -1, -1, -1, -1, -1        

    if (verbose > 1):
        print('    bkg fit param: ', popt)

    residuals = b - func(a, *popt)
    chi2 = sum((residuals / err) ** 2)/len(residuals)

    return func, chi2, np.std(residuals), popt, pcov

def fit_erfi(a, b, err, ts, verbose):

    def func(f, amp, f_mean, f_sigma, time):
        # the 2.pi goes with 'time' because we are in the frequency domain and not in the angular frequency domain
        # time is 'ts-t0'
        return -amp*np.exp(-np.asarray(f-f_mean)**2/(2*abs(f_sigma)**2))*np.imag(special.erfi((np.asarray(f-f_mean)+1j*(abs(f_sigma)**2)*time*2*math.pi)/np.sqrt(2*abs(f_sigma)**2)))

    initial_guess = [10.0, 6705, 10, ts/1000] # tS converted to ms because frequency in kHz
    popt, pcov = curve_fit(func, a, b, p0=initial_guess, maxfev=10000)

    if (verbose > 1):
        print('    bkg fit param: ', popt)

    residuals = b - func(a, *popt)
    chi2 = sum((residuals / err) ** 2)/len(residuals)

    return func, chi2, np.std(residuals), popt, pcov


def fit_pol(a, b, err, poly_order, noise_sigma):

    try:
        fit = np.polyfit(a, b, poly_order, w=err)
    except:
        print("Poly background fit failed: return chi2 of -1")
        return -1, -1, -1

    func = np.poly1d(fit)
    chi2 = np.sum((np.polyval(fit, a) - b) ** 2 /
                  noise_sigma ** 2)/(len(a)-poly_order)
    residuals = []
    for i in range(len(a)):
        residuals.append(np.polyval(fit, a[i])-b[i])

    return func, chi2, np.std(residuals)


def calc_cosine_transform(t0, binContent, binCenter, freq_step_size, n_freq_step, lower_freq):

    a = []
    b = []

    dt = binCenter[1]-binCenter[0]
    for i in range(0, n_freq_step):
        frequency = (lower_freq + freq_step_size/2) + i*freq_step_size
        a.append(frequency)
        frequency /= 1000  # convert from kHz to MHz because time is in micro-sec
        # time hist in micro-sec, dt is 1 ns
        b.append(np.sum(binContent*np.cos(2*math.pi*frequency*(binCenter-t0))*dt))

    return a, b


# compute the sine Fourier transform
def calc_sine_transform(t0, binContent, binCenter, freq_step_size, n_freq_step, lower_freq):

    a = []
    b = []

    dt = binCenter[1]-binCenter[0]
    for i in range(0, n_freq_step):
        frequency = (lower_freq + freq_step_size/2) + i*freq_step_size
        a.append(frequency)
        frequency /= 1000  # convert from kHz to MHz because time is in micro-sec
        # time hist in micro-sec, dt is 1 ns
        b.append(np.sum(binContent*np.sin(2*math.pi*frequency*(binCenter-t0))*dt))

    return a, b


# convert Frequency to Radius ==#
# assume velocity corresponding to magic momentum for all the muons. Velocity depends very little on momentum ==#

def convert_freq_to_radius(freqHist, radius, intensity, n_freq_step):
    for i in range(1, n_freq_step+1):
        radius.append(constants.speedOfLight * constants.magicBeta /
                      (2*math.pi*freqHist.GetBinCenter(i)))
        intensity.append(freqHist.GetBinContent(i))


# compute the parabola correction distribution
def calc_parabola(t0, tS, firstApprox, parabola, n_freq_step, freq_step_size, lower_freq):
    for i in range(0, n_freq_step):
        frq = (lower_freq + freq_step_size/2) + i*freq_step_size  # in MHz
        integral = 0
        for j in range(1, n_freq_step+1):
            if (firstApprox.GetBinCenter(j)-frq) != 0:
                integral += firstApprox.GetBinContent(j)*np.sin(2*math.pi*(
                    frq-firstApprox.GetBinCenter(j))*(tS-t0)/1000)/(1000*(frq-firstApprox.GetBinCenter(j)))
            else:
                integral += 2*math.pi * \
                    firstApprox.GetBinContent(j)*(tS-t0)/1000000
        parabola.SetBinContent(i+1, integral)


# perform the background minimization
def minimization(parabola, cosine, n_freq_step, fit_boundary1, fit_boundary2):

    x = []
    y = []

    for bin_idx in range(1, n_freq_step+1):
        if ((cosine.GetBinCenter(bin_idx) < fit_boundary1 and cosine.GetBinCenter(bin_idx) > constants.lowerCollimatorFreq) or
                (cosine.GetBinCenter(bin_idx) > fit_boundary2 and cosine.GetBinCenter(bin_idx) < constants.upperCollimatorFreq)):
            # if (cosine.GetBinCenter(bin_idx) < constants.lowerCollimatorFreq or cosine.GetBinCenter(bin_idx) > constants.upperCollimatorFreq):
            x.append(parabola.GetBinContent(bin_idx))
            y.append(cosine.GetBinContent(bin_idx))
            x.append(parabola.GetBinContent(bin_idx))
            y.append(cosine.GetBinContent(bin_idx))

    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A, y, rcond=-1)[0]

    return m, c


# manage results directory
def manage_results_directory(dir_path):

    print(' ### Step 1/4: manage results directory\n')

    if os.path.isdir(dir_path) and not os.listdir(dir_path):
        print(
            '    output directoy already exists and is empty --> will save results there\n')
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


def rootHistToNumpArray(hist, tS, tM):

    startBin = hist.FindBin(tS)
    endBin = hist.FindBin(tM)

    binCenter = np.empty(int(endBin-startBin+1), dtype=float)
    binContent = np.empty(int(endBin-startBin+1), dtype=float)

    for j in range(startBin, endBin+1):
        binContent[j-startBin] = hist.GetBinContent(j)
        binCenter[j-startBin] = hist.GetBinCenter(j)

    return binCenter, binContent

#== Convert from ring global radial coordinate to beam local radial coordinate ==#


def global_to_local_radial_coordinate(graph):

    nPoint = graph.GetN()

    for i in range(0, nPoint):

        x, y = r.Double(), r.Double()
        graph.GetPoint(i, x, y)
        graph.SetPoint(i, x-constants.magic_r, y)

#== Compute Radial Mean within collimator aperture ==#


def compute_radial_mean(radius, intensity):

    mean = 0
    sumI = 0

    for x, y in zip(radius, intensity):

        #== Discard data point if radius outside of collimator aperture ==#
        if (x < constants.lowerCollimatorRad or x > constants.upperCollimatorRad):
            continue

        #== Else compute the Mean ==#
        mean += x*y
        sumI += y

    mean /= sumI

    return mean

#== Compute Radial Standard Deviation within collimator aperture ==#


def compute_radial_std(radius, intensity, meanRad, label):

    std = 0
    sumI = 0

    for x, y in zip(radius, intensity):

        #== Discard data point if radius outside of collimator aperture ==#
        if (label == 'ring'):
            if (x < constants.lowerCollimatorRad or x > constants.upperCollimatorRad):
                continue
        if (label == 'beam'):
            if (x < -45 or x > 45):
                continue

        #== Else compute the STD ==#
        sumI += y
        std += y * (x-meanRad) * (x-meanRad)

    std /= sumI
    std = math.sqrt(std)

    return std

#== Compute the E-field correction ==#


def compute_efield_correction(n, mean, std):

    return (- 2 * math.pow(constants.magicBeta, 2) * n * (1-n) * (math.pow(mean, 2) + math.pow(std, 2)) / (math.pow(constants.magic_r, 2)) * 1e9)

#== Extract the two minima of the Cosine Fourier transform ==#


def extractMinima(hist):
    hist.GetXaxis().SetRangeUser(constants.lower_freq, constants.magicFreq)
    min1 = hist.GetMinimum()
    minBinIdx1 = hist.GetMinimumBin()
    hist.GetXaxis().SetRangeUser(constants.magicFreq, constants.upper_freq)
    min2 = hist.GetMinimum()
    minBinIdx2 = hist.GetMinimumBin()
    hist.GetXaxis().SetRangeUser(constants.lower_freq, constants.upper_freq)
    return min1, min2, abs(min1-min2), minBinIdx1, minBinIdx2


'''
import numba as nb
@nb.njit(fastmath=True, parallel=True)
def compute_numba(t0, a, b, it):
    res = np.empty((it, a.shape[0]))
    # res=np.empty((a.shape[0]))
    # res=np.empty((it))
    ita = np.arange(0, it)

    #ntegral=np.zero(0, it).astype(np.float32)
    integral = np.zeros((it))
    #print('size: ', integral.size)
    for i in nb.prange(ita.shape[0]):
        # for i in nb.prange(ita.shape[0]):
        # for i in nb.prange(1):
        # for i in nb.prange(ita.shape[0]):
        #cpt = i
        # print(i)
        frequency = (constants.lower_freq/1000 + constants.freq_step /
                     1000/2) + i*constants.freq_step/1000  # in MHz
        #print(i, frequency)
        # t=ita[i]
        #print(cpt, res[cpt])
        lala = 0.
        for j in range(a.shape[0]):
            #res[i,j]=a[j] * np.cos( 2. * np.pi * it * b[j])
            #res[j] = a[j] * np.cos( 2. * np.pi * it * b[j])
            #integral[i] += a[j] * np.cos( 2. * np.pi * frequency * (b[j]-t0)) * 0.001
            lala += a[j] * np.cos(2. * np.pi * frequency * (b[j]-t0)) * 0.001
            #res[j] = a[j] * np.cos( 2. * np.pi * frequency * (b[j]-t0)) * 0.001
            # print(integral[i])
            #print( a[j] * np.cos( 2. * np.pi * frequency * (b[j]-t0)) * 0.001 )
        # print(lala)
        #print(i, integral[i])
        # print(np.sum(res))
        integral[i] = lala
        #integral[i] = np.sum(res[i])
        # print(integral[i])
        # print(sum)
        # np.sum(res)
        #np.sum( a[j] * np.cos( 2. * np.pi * frequency * (b[j]-t0) )*0.001 )

    return integral
'''
