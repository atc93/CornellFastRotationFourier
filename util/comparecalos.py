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

name = '../results/60h_perCalo.txt'

fileList    = []
histoList   = []
colorList   = [600, 601, 602, 603, 604, 599, 632, 633, 634, 635, 636, 631, 416, 417, 418, 419, 420, 415, 616, 617, 618, 619, 620, 615]

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

#leg = r.TLegend(0.895,0.28,0.95,0.925)
leg = r.TLegend(0.15,0.25,0.23,0.925)

for i in range(0, 24):
    fileList.append( r.TFile('../results/60h_nominal_calo' + str(i+1) + '/results.root') )
    histoList.append( fileList[i].Get('rad') )
    style.setTH1Style( histoList[i], '', 'Radius [mm]', 'Arbitrary units' )
    histoList[i].SetLineColor ( colorList[i] )
    histoList[i].SetMarkerColor ( colorList[i] )
    leg.AddEntry(histoList[i], '  calo ' + str(i+1), "l")
    if ( i == 0 ):
        histoList[i].Draw()
    else:
        histoList[i].Draw("samelp")

graph_min = -0.05
graph_max = 1.05

leg.Draw("same")

c.SaveAs('../results/' + tag + '_radial_allCalos.eps')

print(name)
caloNum = np.loadtxt(name, usecols=(0,))
t0      = np.loadtxt(name, usecols=(2,))*1000
xe      = np.loadtxt(name, usecols=(20,))
width   = np.loadtxt(name, usecols=(22,))
ce      = np.loadtxt(name, usecols=(24,))

avg_ce = 0
avg_xe = 0
avg_w = 0
for i in ce:
    avg_ce += i
for i in xe:
    avg_xe += i
for i in width:
    avg_w += i
avg_ce /= 24
avg_xe /= 24
avg_w /= 24

std_ce = 0
sum = 0
for i in ce:
    sum += 1
    std_ce += (i-avg_ce) * (i-avg_ce)
std_ce /= sum-1
std_ce = math.sqrt(std_ce)

print( 'C_E = ', avg_ce, ' +- ', std_ce )

std_xe = 0
sum = 0
for i in xe:
    sum += 1
    std_xe += (i-avg_xe) * (i-avg_xe)
std_xe /= sum-1
std_xe = math.sqrt(std_xe)

print( 'x_e = ', avg_xe, ' +- ', std_xe )

std_w = 0
sum = 0
for i in width:
    sum += 1
    std_w += (i-avg_w) * (i-avg_w)
std_w /= sum-1
std_w = math.sqrt(std_w)

print( 'width = ', avg_w, ' +- ', std_w )


for i in range(4):

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.text(-0.215, 0.99, label, transform=ax.transAxes, fontsize=12)
    plt.subplots_adjust(left=0.175, bottom=0.125, right=0.975, top=0.95, wspace=0, hspace=0)
    plt.xlabel('Calo #')

    if (i == 5):
        fit = np.polyfit(caloNum,t0,1)
        fit_fn = np.poly1d(fit) 
        plt.plot(caloNum,t0, 'o', caloNum, fit_fn(caloNum), 'k', ms=10, linewidth=2)
        plt.suptitle('$\mathregular{t_{0}}$ = ' + '{0:.2f} * calo # + {1:.2f} ns'.format(fit_fn.c[0], fit_fn.c[1]), y=0.99, x=0.575, fontsize=16)
        plt.ylabel('$\mathregular{t_{0}}$ [ns]')
        plt.savefig('../results/' + tag + '_t0_vs_caloNum.eps', format='eps')
        plt.close()

    if (i == 1):
        fit = np.polyfit(caloNum,t0,1)
        fit_fn = np.poly1d(fit) 
        plt.plot(caloNum,t0-fit_fn(caloNum), 'o', ms=10)
        plt.ylabel('$\mathregular{t_{0}}$ fit residuals [ns]')
        plt.savefig('../results/' + tag + '_t0residuals_vs_caloNum.eps', format='eps')
        plt.close()

    if (i == 2):
        plt.plot(caloNum, xe, 'o', ms=10, label='data', zorder=2)
        plt.subplots_adjust(left=0.125, bottom=0.125, right=0.95, top=0.95, wspace=0, hspace=0)
        plt.suptitle('$\mathregular{<x_{e}>} $= ' + '{0:.2f} $\pm$ {1:.2f} mm'.format(avg_xe, std_xe/math.sqrt(24)), y=0.99, x=0.575, fontsize=16)
        plt.ylabel('$\mathregular{x_{e}}$ [mm]')
        plt.savefig('../results/' + tag + '_xe_vs_caloNum.eps', format='eps')
        plt.close()

    if (i == 3):
        plt.plot(caloNum, width, 'o', ms=10, label='data', zorder=2)
        plt.suptitle('$\mathregular{\sigma} $ = ' + '{0:.2f} $\pm$ {1:.2f} mm'.format(avg_w, std_w/math.sqrt(24)), y=0.99, x=0.575, fontsize=16)
        plt.ylabel('$\mathregular{\sigma}$ [mm]')
        plt.savefig('../results/' + tag + '_width_vs_caloNum.eps', format='eps')
        plt.close()

    if (i == 4):
        plt.plot(caloNum, ce, 'o', ms=10, label='data', zorder=2)
        plt.suptitle('$\mathregular{<C_E>} $ = ' + '{0:.0f} $\pm$ {1:.0f} ppb'.format(avg_ce, std_ce/math.sqrt(24)), y=0.99, x=0.575, fontsize=16)
        plt.ylabel('$\mathregular{C_{E}}$ [ppb]')
        plt.savefig('../results/' + tag + '_ce_vs_caloNum.eps', format='eps')
        plt.close()

    del fig, ax
