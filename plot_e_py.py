# Implementation of the plotting step of the analysis
#
# The plotting combines the histograms to plots which allow us to study the
# inital dataset based on observables motivated through physics.

import ROOT
ROOT.gROOT.SetBatch(True)

# Declare a human-readable label for each variable
labels = {
        "pt_l": "Gen-level Electron p_{T} / GeV",
        "pt_n": "Gen-level Neutrino / GeV",
        "eta_l": "Gen-level Electron #eta",
        "eta_n": "Gen-level Neutrino #eta",
        "phi_l": "Gen-level Electron #phi",
        "phi_n": "Gen-level Neutrino #phi",
        "pt_genmet": "GenMET (p_{T}) / GeV",
        "phi_genmet": "GenMET (#phi)",
        "pdgId_l": "pdgId",
        "pdgId_n": "pdgId",
        #"mass_l": "Electron mass / GeV",
        #"mass_n": "Neutrino mass / GeV",
        "mt": "Gen-level M_{T} / GeV",
        "mt_met": "Gen-level Transverse mass / GeV",
        "m_inv": "Gen-level invariant M(l#nu) / GeV",
        "scalePDF": "invariant M(l#nu) / GeV",
        }

# Specify the color for each process
colors = {
        "enu_100to200": ROOT.TColor.GetColor(100, 192, 232),
        "enu_200to500": ROOT.TColor.GetColor(248, 206, 104),
        "enu_500to1000": ROOT.TColor.GetColor(200, 106, 100),
        "enu_1000to2000": ROOT.TColor.GetColor("#BF2229"),
        "enu_2000to3000": ROOT.TColor.GetColor("#00A88F"),
        "enu_3000to4000": ROOT.TColor.GetColor(155, 152, 204),
        "enu_4000to5000": ROOT.TColor.GetColor(222, 90, 106),
        "enu_5000to6000":  ROOT.TColor.GetColor(250, 202, 255),
        "enu_6000":  ROOT.TColor.GetColor(190, 110, 200),
        }


# Retrieve a histogram from the input file based on the process and the variable
# name
def getHistogram(tfile, name, variable, tag=""):
    name = "{}_{}{}".format(name, variable, tag)
    h = tfile.Get(name)
    if not h:
        raise Exception("Failed to load histogram {}.".format(name))
    return h


# Main function of the plotting step
#
# The major part of the code below is dedicated to define a nice-looking layout.
# The interesting part is the combination of the histograms to the QCD estimation.
# There, we take the data histogram from the control region and subtract all known
# processes defined in simulation and define the remaining part as QCD. Then,
# this shape is extrapolated into the signal region with a scale factor.
def main(variable):
    tfile = ROOT.TFile("histograms_e_pythia.root", "READ")

    # Styles
    ROOT.gStyle.SetOptStat(0)

    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
    ROOT.gStyle.SetCanvasDefH(600)
    ROOT.gStyle.SetCanvasDefW(600)
    ROOT.gStyle.SetCanvasDefX(0)
    ROOT.gStyle.SetCanvasDefY(0)

    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetPadBottomMargin(0.13)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.05)

    ROOT.gStyle.SetHistLineColor(1)
    ROOT.gStyle.SetHistLineStyle(0)
    ROOT.gStyle.SetHistLineWidth(1)
    ROOT.gStyle.SetEndErrorSize(2)
    ROOT.gStyle.SetMarkerStyle(20)

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetTitleFont(42)
    ROOT.gStyle.SetTitleColor(1)
    ROOT.gStyle.SetTitleTextColor(1)
    ROOT.gStyle.SetTitleFillColor(10)
    ROOT.gStyle.SetTitleFontSize(0.05)

    ROOT.gStyle.SetTitleColor(1, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleXOffset(1.00)
    ROOT.gStyle.SetTitleYOffset(1.60)

    ROOT.gStyle.SetLabelColor(1, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetLabelOffset(0.007, "XYZ")
    ROOT.gStyle.SetLabelSize(0.04, "XYZ")

    ROOT.gStyle.SetAxisColor(1, "XYZ")
    ROOT.gStyle.SetStripDecimals(True)
    ROOT.gStyle.SetTickLength(0.03, "XYZ")
    ROOT.gStyle.SetNdivisions(510, "XYZ")
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    ROOT.gStyle.SetPaperSize(20., 20.)
    ROOT.gStyle.SetHatchesLineWidth(5)
    ROOT.gStyle.SetHatchesSpacing(0.05)

    ROOT.TGaxis.SetExponentOffset(-0.08, 0.01, "Y")

    # Simulation
    h100 = getHistogram(tfile, "enu_100to200" , variable)
    h200 = getHistogram(tfile, "enu_200to500" , variable)
    h500 = getHistogram(tfile, "enu_500to1000" , variable)
    h1000 = getHistogram(tfile, "enu_1000to2000" , variable)
    h2000 = getHistogram(tfile, "enu_2000to3000", variable)
    h3000 = getHistogram(tfile, "enu_3000to4000", variable)
    h4000 = getHistogram(tfile, "enu_4000to5000", variable)
    h5000 = getHistogram(tfile, "enu_5000to6000", variable)
    h6000 = getHistogram(tfile, "enu_6000", variable)

    # Draw histograms
    for x, l in [ (h100, "enu_100to200"), (h200, "enu_200to500"), (h500, "enu_500to1000"), (h1000, "enu_1000to2000"), (h2000, "enu_2000to3000"), (h3000, "enu_3000to4000"), (h4000, "enu_4000to5000"), (h5000, "enu_5000to6000"), (h6000, "enu_6000")]:
        x.SetLineWidth(0)
        x.SetFillColor(colors[l])

    stack = ROOT.THStack("", "")
    for x in [h100, h200, h500, h1000, h2000, h3000, h4000, h5000, h6000]:
        stack.Add(x)

    c = ROOT.TCanvas("", "", 600, 600)
    c.SetLogy()
    stack.Draw("hist")
    name = h1000.GetTitle()
    if name in labels:
        title = labels[name]
    else:
        title = name
    stack.GetXaxis().SetTitle(title)
    stack.GetYaxis().SetTitle("N_{Events}")
    stack.SetMaximum(stack.GetMaximum() * 50)
    stack.SetMinimum(0.00000000001)

    #h1000.Draw("HIST SAME")
    #h2000.Draw("HIST SAME")
    #h3000.Draw("HIST SAME")
    #h4000.Draw("HIST SAME")

    # Add legend
    legend = ROOT.TLegend(0.4, 0.73, 0.90, 0.88)
    legend.SetNColumns(2)
    legend.AddEntry(h100  , "W#rightarrowe#nu(M100-200)"  , "f")
    legend.AddEntry(h200  , "W#rightarrowe#nu(M200-500)"  , "f")
    legend.AddEntry(h500  , "W#rightarrowe#nu(M500-1000)" , "f")
    legend.AddEntry(h1000 , "W#rightarrowe#nu(M1000-2000)", "f")
    legend.AddEntry(h2000 , "W#rightarrowe#nu(M2000-3000)", "f")
    legend.AddEntry(h3000 , "W#rightarrowe#nu(M3000-4000)", "f")
    legend.AddEntry(h4000 , "W#rightarrowe#nu(M4000-5000)", "f")
    legend.AddEntry(h5000 , "W#rightarrowe#nu(M5000-6000)", "f")
    legend.AddEntry(h6000 , "W#rightarrowe#nu(M6000-inf)" , "f")
    legend.SetBorderSize(0)
    legend.Draw()

    # Add title
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.SetTextFont(42)
    latex.DrawLatex(0.65, 0.935,  "(13.6 TeV)") #"11.5 fb^{-1} (2022, 13.6 TeV)")
    latex.DrawLatex(0.16, 0.935, "#bf{CMS Simulation}")

    # Save
    c.SaveAs("{}_e_pythia.pdf".format(variable))
    #c.SaveAs("{}.png".format(variable))


# Loop over all variable names and make a plot for each
if __name__ == "__main__":
    for variable in labels.keys():
        main(variable)
