# Implementation of the plotting step of the analysis
#
# The plotting combines the histograms to plots which allow us to study the
# inital dataset based on observables motivated through physics.

import ROOT
ROOT.gROOT.SetBatch(True)

# Declare a human-readable label for each variable
labels = {
        "pt_l": "LHE-level Electron p_{T} / GeV",
        "pt_n": "LHE-level Neutrino / GeV",
        "eta_l": "LHE-level Electron #eta",
        "eta_n": "LHE-level Neutrino #eta",
        "phi_l": "LHE-level Electron #phi",
        "phi_n": "LHE-level Neutrino #phi",
        "pt_genmet": "GenMET (p_{T}) / GeV",
        "phi_genmet": "GenMET (#phi)",
        "pdgId_l": "pdgId",
        "pdgId_n": "pdgId",
        #"m_l": "Electron mass / GeV",
        #"m_n": "Neutrino mass / GeV",
        "mt": "LHE-level M_{T} / GeV",
        "mt_met": "Gen-level Transverse mass / GeV",
        "m_inv": "LHE-level invariant M(l#nu) / GeV",
        }

# Specify the color for each process
colors = {
        "enu_120-200": ROOT.TColor.GetColor(100, 192, 232),
        "enu_200-400": ROOT.TColor.GetColor(248, 206, 104),
        #"enu_400-800": ROOT.TColor.GetColor(200, 106, 100),
        "enu_800-1500": ROOT.TColor.GetColor("#BF2229"),
        "enu_1500-2500": ROOT.TColor.GetColor("#00A88F"),
        "enu_2500-4000": ROOT.TColor.GetColor(155, 152, 204),
        "enu_4000-6000": ROOT.TColor.GetColor(222, 90, 106),
        #"enu_6000-inf":  ROOT.TColor.GetColor(250, 202, 255),
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
    tfile = ROOT.TFile("histograms_e.root", "READ")

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
    h120 = getHistogram(tfile, "enu_120-200" , variable)
    h200 = getHistogram(tfile, "enu_200-400" , variable)
    #h400 = getHistogram(tfile, "enu_400-800" , variable)
    h800 = getHistogram(tfile, "enu_800-1500" , variable)
    h1500 = getHistogram(tfile, "enu_1500-2500", variable)
    h2500 = getHistogram(tfile, "enu_2500-4000", variable)
    h4000 = getHistogram(tfile, "enu_4000-6000", variable)
    #h6000 = getHistogram(tfile, "enu_6000-inf", variable)

    # Draw histograms
    #for x, l in [ (h120, "enu_120-200"), (h200, "enu_200-400"), (h400, "enu_400-800"), (h800, "enu_800-1500"), (h1500, "enu_1500-2500"), (h2500, "enu_2500-4000"), (h4000, "enu_4000-6000")], (h6000, "enu_6000-inf")]:
    for x, l in [ (h120, "enu_120-200"), (h200, "enu_200-400"), (h800, "enu_800-1500"), (h1500, "enu_1500-2500"), (h2500, "enu_2500-4000"), (h4000, "enu_4000-6000")]:
        x.SetLineWidth(0)
        x.SetFillColor(colors[l])

    stack = ROOT.THStack("", "")
    #for x in [h120, h200, h400, h800, h1500, h2500, h4000, h6000]:
    for x in [h120, h200, h800, h1500, h2500, h4000]:
        stack.Add(x)

    c = ROOT.TCanvas("", "", 600, 600)
    c.SetLogy()
    stack.Draw("hist")
    name = h800.GetTitle()
    if name in labels:
        title = labels[name]
    else:
        title = name
    stack.GetXaxis().SetTitle(title)
    stack.GetYaxis().SetTitle("N_{Events}")
    stack.SetMaximum(stack.GetMaximum() * 50)
    stack.SetMinimum(0.00000000001)

    #h800.Draw("HIST SAME")
    #h1500.Draw("HIST SAME")
    #h2500.Draw("HIST SAME")
    #h4000.Draw("HIST SAME")

    # Add legend
    legend = ROOT.TLegend(0.4, 0.73, 0.90, 0.88)
    legend.SetNColumns(2)
    legend.AddEntry(h120  , "W#rightarrowe#nu(M120-200)", "f")
    legend.AddEntry(h200  , "W#rightarrowe#nu(M200-400)", "f")
    #legend.AddEntry(h400  , "W#rightarrowe#nu(M400-800)", "f")
    legend.AddEntry(h800  , "W#rightarrowe#nu(M800-1500)", "f")
    legend.AddEntry(h1500 , "W#rightarrowe#nu(M1500-2500)", "f")
    legend.AddEntry(h2500 , "W#rightarrowe#nu(M2500-4000)", "f")
    legend.AddEntry(h4000 , "W#rightarrowe#nu(M4000-6000)", "f")
    #legend.AddEntry(h6000 , "W#rightarrowe#nu(M6000-inf)", "f")
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
    c.SaveAs("{}_e.pdf".format(variable))
    #c.SaveAs("{}.png".format(variable))


# Loop over all variable names and make a plot for each
if __name__ == "__main__":
    for variable in labels.keys():
        main(variable)
