import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

plt.rcParams['figure.figsize'] = [9, 6]
plt.rcParams.update({'font.size': 17})


inFileName = '../results/60h_ethScan.txt' 
tag = '60h'
label='Run-1 \n60-hour'

eth = np.loadtxt(inFileName, usecols=(1,))
t0 = np.loadtxt(inFileName, usecols=(3,))*1000 # t0
chi2 = np.loadtxt(inFileName, usecols=(5,)) # chi2
noise = np.loadtxt(inFileName, usecols=(7,)) # noise
noise_th = np.loadtxt(inFileName, usecols=(9,)) # noise threshold
ts = np.loadtxt(inFileName, usecols=(13,)) # ts
tm = np.loadtxt(inFileName, usecols=(15,)) # tm
df = np.loadtxt(inFileName, usecols=(17,)) # df
findex = np.loadtxt(inFileName, usecols=(19,)) # field index
xe = np.loadtxt(inFileName, usecols=(21,)) # xe
width = np.loadtxt(inFileName, usecols=(23,)) # width
ce = np.loadtxt(inFileName, usecols=(25,)) # ce

for i in range(6):

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.text(-0.215, 0.99, label, transform=ax.transAxes, fontsize=12)
    plt.subplots_adjust(left=0.175, bottom=0.125, right=0.975, top=0.95, wspace=0, hspace=0)
    plt.xlabel('Positron energy threshold [MeV]')

    if (i == 0):
        ax.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(float(x), '.2f')))
        plt.plot(eth, t0, 'o')
        #ax = plt.gca()
        #ax.set_xlim( 450, 2150 )
        plt.ylabel('$\mathregular{t_{0}}$ [ns]')
        plt.savefig(tag + '_t0_vs_eth.eps', format='eps')
        plt.close()

    if (i == 1):
        plt.plot(eth, xe, 'o', ms=10)
        #ax = plt.gca()
        #ax.set_xlim( 1150, 1950 )
        #ax.set_ylim( 5.94, 6.11 )
        plt.ylabel('$\mathregular{x_{e}}$ [mm]')
        plt.savefig(tag + '_xe_vs_eth.eps', format='eps')
        plt.close()
        
    if (i == 2):
        plt.plot(eth, width, 'o', ms=10)
        #ax = plt.gca()
        #ax.set_xlim( 1150, 1950 )
        #ax.set_ylim( 8.76, 8.83 )
        plt.ylabel('$\mathregular{\sigma}$ [mm]')
        plt.savefig(tag + '_width_vs_eth.eps', format='eps')
        plt.close()
        
    if (i == 3):
        plt.plot(eth, ce, 'o', ms=10)
        #ax = plt.gca()
        #ax.set_xlim( 1150, 1950 )
        #ax.set_ylim( -441, -425 )
        plt.ylabel('$\mathregular{C_{E}}$ [ppb]')
        plt.savefig(tag + '_ce_vs_eth.eps', format='eps')
        plt.close()
        
    if (i == 4):
        plt.plot(eth, noise, 'o', ms=10)
        #ax = plt.gca()
        #ax.set_xlim( 1150, 1950 )
        plt.ylabel('Background fit residuals [a.u.]')
        plt.savefig(tag + '_fitredisuals_vs_eth.eps', format='eps')
        plt.close()
        
    if (i == 5):
        plt.plot(eth, chi2, 'o', ms=10)
        #ax = plt.gca()
        #ax.set_xlim( 1150, 1950 )
        plt.ylabel('$\mathregular{{\chi}^2}}$/d.o.f.')
        plt.savefig(tag + '_chi2_vs_eth.eps', format='eps')
        plt.close()

    del fig, ax
