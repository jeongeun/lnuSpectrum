# Implementation of the histogramming step of the analysis
#
# The histogramming step produces histograms for each variable in the dataset
# and for each physics process resulting into the final state with a muon and a
# tau. Then, the resulting histograms are passed to the plotting step, which
# combines the histograms so that we can study the physics of the decay.

import ROOT
ROOT.gROOT.SetBatch(True)

# Declare the range of the histogram for each variable
#
# Each entry in the dictionary contains of the variable name as key and a tuple
# specifying the histogram layout as value. The tuple sets the number of bins,
# the lower edge and the upper edge of the histogram.
default_nbins = 30

ranges = {
        "pt_l": (100, 0, 5000),
        "pt_n": (100, 0, 5000),
        "eta_l": (default_nbins, -2.5, 2.5),
        "eta_n": (default_nbins, -2.5, 2.5),
        "phi_l": (default_nbins, -3.14, 3.14),
        "phi_n": (default_nbins, -3.14, 3.14),
        "pdgId_l": (40, -20, 20),
        "pdgId_n": (40, -20, 20),
        "pt_genmet": (100, 0, 5000),
        "phi_genmet": (default_nbins, -3.14, 3.14),
        "mass_l": (default_nbins, 0, 2),
        "mass_n": (default_nbins, 0, 2),
        "mt": (400, 0, 8000),
        "mt_met": (400, 0, 8000),
        "m_inv": (400, 0, 8000),
        "scalePDF": (400, 0, 8000),
        }

# Book a histogram for a specific variable
def bookHistogram(df, variable, range_):
    return df.Histo1D(ROOT.ROOT.RDF.TH1DModel(variable, variable, range_[0], range_[1], range_[2]),\
                      variable, "weight")


# Write a histogram with a given name to the output ROOT file
def writeHistogram(h, name):
    h.SetName(name)
    h.Write()

# Main function of the histogramming step
#
# The function loops over the outputs from the skimming step and produces the
# required histograms for the final plotting.
def main():
    # Set up multi-threading capability of ROOT
    ROOT.ROOT.EnableImplicitMT()
    poolSize = ROOT.ROOT.GetImplicitMTPoolSize()
    print("Pool size: {}".format(poolSize))

    # Create output file
    tfile = ROOT.TFile("histograms_m_pythia.root", "RECREATE")
    variables = ranges.keys()

    # Loop through skimmed datasets and produce histograms of variables
    for name, label in [
            ("100to200"  , "mnu_100to200"),
            ("200to500"  , "mnu_200to500"),
            ("500to1000"  , "mnu_500to1000"),
            ("1000to2000"  , "mnu_1000to2000"),
            ("2000to3000" , "mnu_2000to3000"),
            ("3000to4000" , "mnu_3000to4000"),
            ("4000to5000" , "mnu_4000to5000"),
            ("5000to6000" , "mnu_5000to6000"),
            ("6000"  , "mnu_6000"),
        ]:
        print(">>> Process skimmed sample {} for process {}".format(name, label))

        # Load skimmed dataset
        df = ROOT.ROOT.RDataFrame("Events", name + "_Skim_m_pythia.root")

        # Book histograms
        hists = {}
        for variable in variables:
            hists[variable] = bookHistogram(df, variable, ranges[variable])
        report1 = df.Report()

        # Write histograms to output file
        for variable in variables:
            writeHistogram(hists[variable], "{}_{}".format(label, variable))

        # Print cut-flow report
        print("Cut-flow report (signal region):")
        report1.Print()
    tfile.Close()

if __name__ == "__main__":
    main()
