#==============================================#
#== IMPORT FILE CONTAINING STYLING FUNCTIONS ==#
#==============================================#

import ROOT as r
import src.constants as constants


#== Set TCanvas style ==#
def setTCanvasStyle(c):

    r.gStyle.SetOptStat(0)
    r.gStyle.SetOptFit(0)
    r.gStyle.SetOptTitle(1)
    r.gStyle.SetStatX(1)
    r.gStyle.SetStatY(1)
    r.gStyle.SetStatH(0.1)
    r.gStyle.SetStatW(0.15)
    r.gPad.SetTicks(1)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.075)
    c.SetBottomMargin(0.15)


#== Set TH1 style ==#
def setTH1Style(h, title, xAxisTitle, yAxisTitle, *args):

    h.SetTitle(title)
    h.GetXaxis().CenterTitle()
    h.GetXaxis().SetTitle(xAxisTitle)
    h.GetYaxis().CenterTitle()
    h.GetYaxis().SetTitle(yAxisTitle)
    h.GetXaxis().SetTitleOffset(1.4)
    h.GetYaxis().SetTitleOffset(1.4)
    h.GetXaxis().SetTitleSize(0.055)
    h.GetXaxis().SetLabelSize(0.05)
    h.GetYaxis().SetTitleSize(0.055)
    h.GetYaxis().SetLabelSize(0.05)
    h.SetLineColor(1)
    h.SetMarkerColor(1)
    h.SetLineWidth(2)

    # == The argument list args contains
    # == args[0]: scale factor for histogram minimum
    # == args[1]: scale factor for histogram maximum
    if (len(args) == 2):
        h.SetMinimum(h.GetMinimum()*args[0])
        h.SetMaximum(h.GetMaximum()*args[1])


#== Set TGraph style ==#
def setTGraphStyle(graph, title, xAxisTitle, yAxisTitle, graphMin, graphMax):

    graph.SetTitle(title)
    graph.GetXaxis().SetTitle(xAxisTitle)
    graph.GetYaxis().SetTitle(yAxisTitle)
    graph.GetXaxis().CenterTitle()
    graph.GetYaxis().CenterTitle()
    graph.GetXaxis().SetTitleOffset(1.4)
    graph.GetXaxis().SetTitleSize(0.055)
    graph.GetXaxis().SetLabelSize(0.05)
    graph.GetYaxis().SetTitleOffset(1.4)
    graph.GetYaxis().SetTitleSize(0.055)
    graph.GetYaxis().SetLabelSize(0.05)
    graph.GetXaxis().SetRangeUser(7052, 7172)
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(0.9)
    graph.SetMarkerColor(1)
    graph.SetLineColor(1)
    graph.SetMaximum(graphMax)
    graph.SetMinimum(graphMin)


#== Create and set vertical TLines style showing the collimator frequency aperture ==#
def set_collimator_aperture_tline(yMin, yMax, label):

    if (label == 'frequency'):
        innerLine = r.TLine(constants.lowerCollimatorFreq,
                            yMin, constants.lowerCollimatorFreq, yMax)
        outerLine = r.TLine(constants.upperCollimatorFreq,
                            yMin, constants.upperCollimatorFreq, yMax)
    if (label == 'radial'):
        innerLine = r.TLine(constants.lowerCollimatorRad,
                            yMin, constants.lowerCollimatorRad, yMax)
        outerLine = r.TLine(constants.upperCollimatorRad,
                            yMin, constants.upperCollimatorRad, yMax)

    innerLine.SetLineWidth(3)
    outerLine.SetLineWidth(3)

    return innerLine, outerLine


#== Create and set TPaveText style for radial collimator aperture ==#
def setCollimatorAperturePaveText(yMin, yMax, label):

    if (label == 'frequency'):
        pt = r.TPaveText(constants.lowerCollimatorTextFreq1,
                         yMin, constants.lowerCollimatorTextFreq2, yMax)
        pt2 = r.TPaveText(constants.upperCollimatorTextFreq1,
                          yMin, constants.upperCollimatorTextFreq2, yMax)
    if (label == 'radial'):
        pt = r.TPaveText(constants.lowerCollimatorTextRad1,
                         yMin, constants.lowerCollimatorTextRad2, yMax)
        pt2 = r.TPaveText(constants.upperCollimatorTextRad1,
                          yMin, constants.upperCollimatorTextRad2, yMax)

    pt.AddText("collimators")
    pt.AddText("aperture")
    pt.SetShadowColor(0)
    pt.SetBorderSize(1)
    pt.SetFillColor(0)
    pt.SetLineWidth(1)
    pt.SetLineColor(1)
    pt.SetTextAngle(90)

    pt2.AddText("collimators")
    pt2.AddText("aperture")
    pt2.SetShadowColor(0)
    pt2.SetBorderSize(1)
    pt2.SetFillColor(0)
    pt2.SetLineWidth(1)
    pt2.SetLineColor(1)
    pt2.SetTextAngle(90)

    return pt, pt2

#== Create and set vertical TLine style showing the magic radius ==#


def setMagicRadiusTLine(yMin, yMax):

    magicLine = r.TLine(constants.magic_r, yMin, constants.magic_r, yMax)
    magicLine.SetLineWidth(1)
    magicLine.SetLineStyle(7)

    return magicLine

#== Create and set TPaveText style for magic radius vertical TLine ==#


def setMagicRadiusPaveText(yMin, yMax):

    pt = r.TPaveText(7113, yMin, 7121, yMax)

    pt.AddText("magic")
    pt.AddText("radius")
    pt.SetShadowColor(0)
    pt.SetBorderSize(1)
    pt.SetFillColor(0)
    pt.SetLineWidth(1)
    pt.SetLineColor(1)

    return pt

#== Create and set TPaveText style containing results: x_e, width, c_e ==#


def set_radial_results_pave_text(eq_radius, std, c_e, yMin, yMax, radMin, radMax, coord):

    pt = r.TPaveText(radMin, yMax, radMax, yMin)

    if (coord == 'beam'):
        pt.AddText('x_{e} = ' + '{0:.2f}'.format(eq_radius) + ' mm')
        pt.AddText('#sigma = ' + '{0:.2f}'.format(std) + ' mm')
        pt.AddText('  C_{E} = ' + '{0:.0f}'.format(c_e) + ' ppb  ');
    elif (coord == 'ring'):
        pt.AddText('x_{e} = ' + '{0:.2f}'.format(eq_radius) + ' mm')
        pt.AddText('#sigma = ' + '{0:.2f}'.format(std) + ' mm   ')

    pt.SetShadowColor(0)
    pt.SetBorderSize(1)
    pt.SetFillColor(0)
    pt.SetLineWidth(1)
    pt.SetLineColor(1)
    pt.SetTextAngle(90)

    if (coord == 'truth'):
        pt.AddText('Truth level');
        pt.AddText('x_{e} = ' + '{0:.2f}'.format(eq_radius) + ' mm');
        pt.AddText(' #sigma = ' + '{0:.2f}'.format(std) + ' mm ');
        pt.AddText('      C_{E} = ' + '{0:.0f}'.format(c_e) + ' ppb     ');
        pt.SetLineColor(2);
        pt.SetTextColor(2);

    return pt

#== Create and set TPaveText style to compare rec/truth-level mean frequencies ==#


def setRecTruthFrequenciesPaveText(recFreq, trueFreq, rec_width, true_width, yMin, yMax):

    pt = r.TPaveText(6664, yMax, 6686, yMin)
    pt.AddText('<reco freq> = ' + '{0:.2f}'.format(recFreq) + ' kHz')
    pt.AddText('  <true freq> = ' + '{0:.2f}'.format(trueFreq) + ' kHz ')
    pt.GetListOfLines().Last().SetTextColor(2)
    pt.AddText('  reco width = ' + '{0:.2f}'.format(rec_width) + ' kHz    ')
    pt.AddText(' true width = ' + '{0:.2f}'.format(true_width) + ' kHz   ')
    pt.GetListOfLines().Last().SetTextColor(2)
    pt.SetShadowColor(0)
    pt.SetBorderSize(1)
    pt.SetFillColor(0)
    pt.SetLineWidth(1)
    pt.SetLineColor(1)
    pt.SetTextAngle(90)

    return pt

def setRecFrequenciesPaveText(freq, width, yMin, yMax):

    pt = r.TPaveText(6664, yMax, 6686, yMin)
    pt.AddText('<rfreq> = ' + '{0:.2f}'.format(freq) + ' kHz')
    pt.AddText('  width = ' + '{0:.2f}'.format(width) + ' kHz    ')
    pt.GetListOfLines().Last().SetTextColor(2)
    pt.SetShadowColor(0)
    pt.SetBorderSize(1)
    pt.SetFillColor(0)
    pt.SetLineWidth(1)
    pt.SetLineColor(1)
    pt.SetTextAngle(90)

    return pt
