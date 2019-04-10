#=======================================================#
#== IMPORT FILE CONTAINING STAT FLUCTUATION FUNCTIONS ==#
#=======================================================#

import math
from random import randint
import numpy as np

#== Square root of N fluctuation ==#
def squareRootN(  h ):
    nBins = h.GetXaxis().GetNbins()
    for iBin in range ( 1, nBins+1 ):
        binContent = h.GetBinContent( iBin )
        fluc = math.sqrt( binContent )
        rand = randint(0, 1)
        if ( rand > 0.5 ):
            h.SetBinContent( iBin, binContent + fluc )
        else:
            h.SetBinContent( iBin, binContent - fluc )

    return h

#== Poisson fluctuation ==#
def poisson( h ):
    nBins = h.GetXaxis().GetNbins()
    for iBin in range ( 1, nBins+1 ):
        binContent = h.GetBinContent( iBin )
        fluc = math.sqrt( np.random.poisson( binContent ) )
        rand = randint(0, 1)
        if ( rand > 0.5 ):
            h.SetBinContent( iBin, binContent + fluc )
        else:
            h.SetBinContent( iBin, binContent - fluc )
    return h
