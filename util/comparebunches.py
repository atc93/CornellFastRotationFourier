import sys
import math
import ROOT as r
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('../')
import src.style as style

plt.rcParams['figure.figsize'] = [9, 6]
plt.rcParams.update({'font.size': 17})

tag='60h'
label='Run-1 \n60-hour'

name = '../results/60h_perBunch.txt'

fileList    = []
histoList   = []
colorList = [1, 2, 4, 6, 7, 8, 9, 28]

statCeMean = []
statXeMean = []
statWMean  = []
statT0Mean = []
statCeStd = []
statXeStd = []
statWStd  = []
statT0Std = []

c = r.TCanvas('c','c',900,600)
style.setTCanvasStyle( c )

leg = r.TLegend(0.8535,0.5,0.95,0.925)

for i in range(0, 8):
    fileList.append( r.TFile('../results/60h_nominal_bunch' + str(i) + '/results.root') )
    histoList.append( fileList[i].Get('rad') )
    style.setTH1Style( histoList[i], '', 'Radius [mm]', 'Arbitrary' )
    histoList[i].SetLineColor ( colorList[i] )
    histoList[i].SetMarkerColor ( colorList[i] )
    leg.AddEntry(histoList[i], 'bunch ' + str(i), "l")
    if ( i == 0 ):
        histoList[i].Draw()
    else:
        histoList[i].Draw("samelp")

graph_min = -0.05
graph_max = 1.05

leg.Draw("same")

c.SaveAs('../results/' + tag + '_radial_allBunches.eps')

bunchNum = np.loadtxt(name, usecols=(0,))
t0      = np.loadtxt(name, usecols=(2,))*1000
xe      = np.loadtxt(name, usecols=(20,))
width   = np.loadtxt(name, usecols=(22,))
ce      = np.loadtxt(name, usecols=(24,))

avg_ce = 0
avg_xe = 0
avg_w = 0
avg_t0 = 0
for i in ce:
    avg_ce += i
for i in xe:
    avg_xe += i
for i in width:
    avg_w += i
for i in t0:
    avg_t0 += i
avg_ce /= 8
avg_xe /= 8
avg_w /= 8
avg_t0 /= 8

std_ce = 0
sum = 0
for i in ce:
    sum += 1
    std_ce += (i-avg_ce) * (i-avg_ce)
std_ce /= sum-1
std_ce = math.sqrt(std_ce)

std_xe = 0
sum = 0
for i in xe:
    sum += 1
    std_xe += (i-avg_xe) * (i-avg_xe)
std_xe /= sum-1
std_xe = math.sqrt(std_xe)

std_w = 0
sum = 0
for i in width:
    sum += 1
    std_w += (i-avg_w) * (i-avg_w)
std_w /= sum-1
std_w = math.sqrt(std_w)

std_t0 = 0
sum = 0
for i in t0:
    sum += 1
    std_t0 += (i-avg_t0) * (i-avg_t0)
std_t0 /= sum-1
std_t0 = math.sqrt(std_t0)

print( 'C_E = ', avg_ce, ' +- ', std_ce, ' standrd deviation' )
print( 'C_E = ', avg_ce, ' +- ', std_ce/math.sqrt(8), ' standard error' )

print( 'x_e = ', avg_xe, ' +- ', std_xe, ' standard deviation' )
print( 'x_e = ', avg_xe, ' +- ', std_xe/math.sqrt(8), ' standard error' )

print( 'width = ', avg_w, ' +- ', std_w, ' standard deviation' )
print( 'width = ', avg_w, ' +- ', std_w/math.sqrt(8), ' standard error' )

print( 't0 = ', avg_t0, ' +- ', std_t0, ' standard deviation'  )
print( 't0 = ', avg_t0, ' +- ', std_t0/math.sqrt(8), ' standard error'  )


for i in range(4):

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.text(-0.215, 0.99, label, transform=ax.transAxes, fontsize=12)
    plt.subplots_adjust(left=0.175, bottom=0.125, right=0.975, top=0.95, wspace=0, hspace=0)
    plt.xlabel('Bunch #')

    if (i == 0):
        plt.suptitle('$\mathregular{<t_{0}>}=' + '{0:.1f}'.format(avg_t0) + ' \pm$ ' + '{0:.1f}'.format(std_t0/math.sqrt(8)) + ' ns', y=0.99, x=0.575, fontsize=16)
        plt.plot(bunchNum,t0, 'o', ms=12.5)
        plt.ylabel('$\mathregular{t_{0}}$ [ns]')
        plt.savefig('../results/' + tag + '_t0_vs_bunchNum.eps', format='eps')
        plt.close()

    if (i == 1):
        plt.suptitle('$\mathregular{<x_{e}>}=' + '{0:.1f}'.format(avg_xe) + ' \pm$ ' + '{0:.1f}'.format(std_xe/math.sqrt(8)) + ' mm', y=0.99, x=0.575, fontsize=16)
        plt.plot(bunchNum, xe, 'o', ms=12.5)
        plt.ylabel('$\mathregular{x_{e}}$ [mm]')
        plt.savefig('../results/' + tag + '_xe_vs_bunchNum.eps', format='eps')
        plt.close()

    if (i == 2):
        plt.plot(bunchNum, width, 'o', ms=12.5)
        plt.suptitle('$\mathregular{<\sigma>}=' + '{0:.1f}'.format(avg_w) + ' \pm$ ' + '{0:.1f}'.format(std_w/math.sqrt(8)) + ' mm', y=0.99, x=0.575, fontsize=16)
        plt.ylabel('$\mathregular{\sigma}$ [mm]')
        plt.savefig('../results/' + tag + '_width_vs_bunchNum.eps', format='eps')
        plt.close()

    if (i == 3):
        plt.plot(bunchNum, ce, 'o', ms=12.5)
        plt.suptitle('$\mathregular{<C_{E}>}=' + '{0:.1f}'.format(avg_ce) + ' \pm$ ' + '{0:.1f}'.format(std_ce/math.sqrt(8)) + ' ppb', y=0.99, x=0.575, fontsize=16)
        plt.ylabel('$\mathregular{C_{E}}$ [ppb]')
        plt.savefig('../results/' + tag + '_ce_vs_bunchNum.eps', format='eps')
        plt.close()

    del fig, ax
