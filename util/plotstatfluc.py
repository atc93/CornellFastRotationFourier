import sys 
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats.stats import pearsonr   

plt.rcParams['figure.figsize'] = [9, 6]
plt.rcParams.update({'font.size': 17})

inFileName = '/home/achapela/mountCornell/CornellFastRotationFourier/results/60h_statfluc_t0bkgth1/results.txt'
tag = '../results/60h_statfluc_t0bkgth1'
label='Run-1 \n60-hour'

#text_file = open(str(outFileName), "w")

t0 = np.loadtxt(inFileName, usecols=(1,))*1000 # t0
chi2 = np.loadtxt(inFileName, usecols=(3,)) # chi2
noise = np.loadtxt(inFileName, usecols=(5,)) # noise
noise_th = np.loadtxt(inFileName, usecols=(7,)) # noise threshold
ts = np.loadtxt(inFileName, usecols=(11,)) # ts
tm = np.loadtxt(inFileName, usecols=(13,)) # tm
df = np.loadtxt(inFileName, usecols=(15,)) # df
findex = np.loadtxt(inFileName, usecols=(17,)) # field index
xe = np.loadtxt(inFileName, usecols=(19,)) # xe
width = np.loadtxt(inFileName, usecols=(21,)) # width
ce = np.loadtxt(inFileName, usecols=(23,)) # ce


mean_t0 = 0
std_t0  = 0
mean_ce = 0
std_ce  = 0
mean_xe = 0
std_xe  = 0
mean_width = 0
std_width  = 0

for i in range(0, len(ce)):
    mean_ce += ce[i]
    mean_xe += xe[i]
    mean_width += width[i]
    mean_t0 += t0[i]

mean_ce /= len(ce)
mean_xe /= len(ce)
mean_width /= len(ce)
mean_t0 /= len(ce)

for i in range(0, len(ce)):
    std_ce += ( ce[i] - mean_ce ) * ( ce[i] - mean_ce )
    std_xe += ( xe[i] - mean_xe ) * ( xe[i] - mean_xe )
    std_width += ( width[i] - mean_width ) * ( width[i] - mean_width )
    std_t0 += ( t0[i] - mean_t0 ) * ( t0[i] - mean_t0 )

std_ce /= len(ce)-1
std_ce = math.sqrt(std_ce)
std_xe /= len(ce)-1
std_xe = math.sqrt(std_xe)
std_width /= len(ce)-1
std_width = math.sqrt(std_width)
std_t0 /= len(ce)-1
std_t0 = math.sqrt(std_t0)

print ( 'x_e = ' + str(mean_xe) + ' +- ' + str(std_xe) )
print ( 'width ' + str(mean_width) + ' +- ' + str(std_width) )
print ( 'C_E = ' + str(mean_ce) + ' +- ' + str(std_ce) )

#text_file.write('xe %f %f width %f %f ce %f %f t0 %f %f \n' %
#                (mean_xe, std_xe, mean_width, std_width, mean_ce, std_ce, mean_t0, std_t0) )
#text_file.close()


print('correlation: ', np.corrcoef(xe,width))
print('covariance: ', np.cov(xe,width))

#mask = (xe < 6.06)
mask = (xe < 1000)
ce = ce[mask]

for i in range(11):

    fig = plt.figure(1)
    plt.subplots_adjust(left=0.175, bottom=0.125, right=0.975, top=0.95, wspace=0, hspace=0)
    if (i < 7):
        ax = fig.add_subplot(111)
        ax.text(-0.215, 0.99, label, transform=ax.transAxes, fontsize=12)
    else:
        ax = fig.add_subplot(111)
        ax.text(-0.27, 0.99, label, transform=ax.transAxes, fontsize=12)

    if (i == 0):
        plt.xlabel('$\mathregular{C_{E}}$ [ppb]')
        plt.ylabel('#')
        plt.suptitle('$\mathregular{<C_{E}>}$ = ' + '{0:.1f}'.format(mean_ce) + ' $\pm$ ' + '{0:.1f}'.format(std_ce) + ' ppb', y=0.99, x=0.575, fontsize=16)
        plt.hist(ce, bins=int(math.sqrt(len(ce))))
        plt.savefig(tag + '_ce.eps', format='eps')
        plt.close()

    if (i == 1):
        plt.xlabel('$\mathregular{x_{e}}$ [mm]')
        plt.ylabel('#')
        plt.suptitle('$\mathregular{<x_{e}>}$=' + '{0:.2f}'.format(mean_xe) + ' +- ' + '{0:.2f}'.format(std_xe) + ' mm', y=0.99, x=0.575, fontsize=16)
        plt.hist(xe, bins=int(math.sqrt(len(ce))))
        plt.savefig(tag + '_xe.eps', format='eps')
        plt.close()
        
    if (i == 2):
        plt.xlabel('$\mathregular{\sigma}$ [mm]')
        plt.ylabel('#')
        plt.suptitle('$\mathregular{<\sigma>}$=' + '{0:.2f}'.format(mean_width) + ' +- ' + '{0:.2f}'.format(std_width) + ' mm', y=0.99, x=0.575, fontsize=16)
        plt.hist(width, bins=int(math.sqrt(len(ce))))
        plt.savefig(tag + '_width.eps', format='eps')
        plt.close()
        
    if (i == 3):
        plt.xlabel('$\mathregular{t_{0}}$ [ns]')
        plt.ylabel('#')
        plt.suptitle('$\mathregular{<t_{0}>}$=' + '{0:.2f}'.format(mean_t0) + ' +- ' + '{0:.2f}'.format(std_t0) + ' ns', y=0.99, x=0.575, fontsize=16)
        plt.hist(t0, bins=int(math.sqrt(len(ce))))
        plt.savefig(tag + '_t0.eps', format='eps')
        plt.close()
        
    if (i == 4):
        plt.xlabel('$\mathregular{x_{e}}$ [mm]')
        plt.ylabel('$\mathregular{\sigma}$ [mm]')
        plt.suptitle('Correlation coefficient = {0:.3f}'.format(np.corrcoef(xe,width)[1][0]), y=0.99, x=0.575, fontsize=16)
        plt.plot(xe, width, 'o')
        plt.savefig(tag + '_widthVSxe.eps', format='eps')
        plt.close()

    if (i == 5):
        plt.xlabel('$\mathregular{t_{0}}$ [ns]')
        plt.ylabel('$\mathregular{\sigma}$ [mm]')
        plt.suptitle('Correlation coefficient = {0:.3f}'.format(np.corrcoef(t0, width)[1][0]), y=0.99, x=0.575, fontsize=16)
        plt.plot(t0, width, 'o')
        plt.savefig(tag + '_widthVSt0.eps', format='eps')
        plt.close()
        
    if (i == 6):
        plt.xlabel('$\mathregular{t_{0}}$ [ns]')
        plt.ylabel('$\mathregular{x_{e}}$ [mm]')
        plt.suptitle('Correlation coefficient = {0:.3f}'.format(np.corrcoef(t0, xe)[1][0]), y=0.99, x=0.575, fontsize=16)
        plt.plot(t0, xe, 'o')
        plt.savefig(tag + '_xeVSt0.eps', format='eps')
        plt.close()

    if (i == 7):
        #ax = plt.gca()
        #ax.set_ylim( 8.76, 8.89 )
        #ax.set_xlim( 5.6, 5.8 )
        line =ax.scatter(xe[mask], width[mask], c=abs(ce), s=100, marker='o', edgecolors='none', cmap='gist_rainbow')
        cb = plt.colorbar(line)
        cb.set_label('$|\mathregular{C_{e}}|$ [ppb]', rotation=90)
        ax.set_xlabel('$\mathregular{x_{e}}$ [mm]')
        ax.set_ylabel('$\mathregular{\sigma}$ [mm]')
        plt.savefig(tag + '_ceVSwithVSxe.eps', format='eps')
        plt.close()

    if (i == 8):
        line =ax.scatter(xe, width, c=t0-mean_t0, s=100, marker='o', edgecolors='none', cmap='gist_rainbow')
        cb = plt.colorbar(line)
        cb.set_label('$\mathregular{\Delta t_{0}}$ [ns]', rotation=90)
        ax.set_xlabel('$\mathregular{x_{e}}$ [mm]')
        ax.set_ylabel('$\mathregular{\sigma}$ [mm]')
        plt.savefig(tag + '_t0VSwithVSxe.eps', format='eps')
        plt.close()
        
    if (i == 9):
        plt.ticklabel_format(axis='z', style='sci')
        plt.subplots_adjust(left=0.175, bottom=0.125, right=0.94, top=0.95, wspace=0, hspace=0)
        line =ax.scatter(xe, width, c=noise, s=100, marker='o', edgecolors='none', cmap='gist_rainbow')
        cb = plt.colorbar(line, format='%.1e')
        cb.set_label('Background fit residuals [a.u.]', rotation=90)
        ax.set_xlabel('$\mathregular{x_{e}}$ [mm]')
        ax.set_ylabel('$\mathregular{\sigma}$ [mm]')
        plt.savefig(tag + '_residualsVSwithVSxe.eps', format='eps')
        plt.close()
        
    if (i == 10):
        #ax = plt.gca()
        #ax.set_ylim( 8.76, 8.89 )
        #ax.set_xlim( 5.6, 5.8 )
        line =ax.scatter(xe, width, c=chi2, s=100, marker='o', edgecolors='none', cmap='gist_rainbow')
        cb = plt.colorbar(line)
        cb.set_label('$\mathregular{\chi^2}}$/d.o.f.', rotation=90)
        ax.set_xlabel('$\mathregular{x_{e}}$ [mm]')
        ax.set_ylabel('$\mathregular{\sigma}$ [mm]')
        plt.savefig(tag + '_chi2VSwithVSxe.eps', format='eps')
        plt.close()

    del fig, ax
