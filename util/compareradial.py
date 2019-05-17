import matplotlib.pyplot as plt
import numpy as np
import ROOT as r
import sys
import os

sys.path.append('../')

import src.style as style

plt.rcParams['figure.figsize'] = [9, 6]
plt.rcParams.update({'font.size': 17})

n_file = 7

tag = []
tag.append('24376')
tag.append('24448')
tag.append('24582')
tag.append('24835')
tag.append('25026')
tag.append('25940')
tag.append('26150')

file_name = []
file_name.append('../results/run2/run2_' + tag[0] + '/results.root')
file_name.append('../results/run2/run2_' + tag[1] + '/results.root')
file_name.append('../results/run2/run2_' + tag[2] + '/results.root')
file_name.append('../results/run2/run2_' + tag[3] + '/results.root')
file_name.append('../results/run2/run2_' + tag[4] + '/results.root')
file_name.append('../results/run2/run2_' + tag[5] + '/results.root')
file_name.append('../results/run2/run2_' + tag[6] + '/results.root')

x_e = []
x_e.append(5.3)
x_e.append(5.3)
x_e.append(5.1)
x_e.append(5.5)
x_e.append(5.6)
x_e.append(5.6)
x_e.append(5.8)

width = []
width.append(8.5)
width.append(8.6)
width.append(8.5)
width.append(8.7)
width.append(8.7)
width.append(8.7)
width.append(8.7)

c_e =[]
c_e.append(-381)
c_e.append(-390)
c_e.append(-376)
c_e.append(-402)
c_e.append(-408)
c_e.append(-405)
c_e.append(-415)

c = r.TCanvas('c','c',900,600)
style.setTCanvasStyle(c)

leg = r.TLegend(0.15,0.6,0.45,0.925)

color = []
color.append(1)
color.append(2)
color.append(4)
color.append(6)
color.append(8)
color.append(9)
color.append(28)

for i in range(n_file):

    root_file = r.TFile(file_name[i])
    graph = root_file.Get('rad')
    graph.SetLineColor(color[i])
    graph.SetMarkerColor(color[i])
    leg.AddEntry(graph, tag[i] + ', x_{e}=' + str(x_e[i]) + ' mm, #sigma=' + str(width[i]) + ' mm, C_{E}=' + str(c_e[i]) + ' ppb     ', 'pl')

    if (i == 0):

        graph.SetTitle('')
        graph.GetXaxis().SetTitle("Radius [mm]")
        graph.GetYaxis().SetTitle("Arbitrary units")
        graph.GetXaxis().CenterTitle()
        graph.GetYaxis().CenterTitle()
        graph.GetXaxis().SetTitleOffset(1.4)
        graph.GetXaxis().SetTitleSize(0.055);
        graph.GetXaxis().SetLabelSize(0.05);
        graph.GetYaxis().SetTitleOffset(1.4)
        graph.GetYaxis().SetTitleSize(0.055);
        graph.GetYaxis().SetLabelSize(0.05)    
        graph.Draw()

    else:

        graph.Draw('samepl')

leg.Draw('same')
c.SaveAs('test.eps')     
