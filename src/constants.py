#==============================================#
#== IMPORT FILE CONTAINING VARIOUS CONSTANTS ==#
#==============================================#

import math

#== Physics constants ==#
muonMass = .105658  # in GeV
magicR = 7112  # in mm
magicP = 3.094  # in GeV
magicE = math.sqrt(math.pow(muonMass, 2) + math.pow(magicP, 2))  # in GeV
magicGamma = magicE / muonMass  # dimensionless
magicBeta = magicP / (magicGamma * muonMass)  # dimensionless
speedOfLight = 299792458  # in m/s
magicFreq = magicBeta*speedOfLight / (2 * math.pi * (magicR))  # MHz

#== For plotting purposes ==#

##== Collimator frequency boundaries to draw line ==##
lowerCollimatorFreq = float(6662.799)  # in kHz
upperCollimatorFreq = float(6747.651)  # in kHz

##== Collimator frequency boundaries to print out text ==##
lowerCollimatorTextFreq1 = float(6650)  # in kHz
lowerCollimatorTextFreq2 = float(6674)  # in kHz
upperCollimatorTextFreq1 = float(6737)  # in kHz
upperCollimatorTextFreq2 = float(6759)  # in kHz

##== Collimator radius boundaries to draw line ==##
lowerCollimatorRad = float(7112-45)  # in mm
upperCollimatorRad = float(7112+45)  # in mm

##== Collimator radial boundaries to print out text ==##
lowerCollimatorTextRad1 = float(7060)  # in kHz
lowerCollimatorTextRad2 = float(7074)  # in kHz
upperCollimatorTextRad1 = float(7150)  # in kHz
upperCollimatorTextRad2 = float(7164)  # in kHz

#== For analysis purposes ==#

##== Noise threshold in unit of statistical standard deviation ==#
noiseThreshold = 2
