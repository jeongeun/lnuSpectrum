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
        "pt_l": (500, 0, 5000),
        "pt_n": (500, 0, 5000),
        "eta_l": (default_nbins, -2.5, 2.5),
        "eta_n": (default_nbins, -2.5, 2.5),
        "phi_l": (default_nbins, -3.14, 3.14),
        "phi_n": (default_nbins, -3.14, 3.14),
        "pdgId_l": (40, -20, 20),
        "pdgId_n": (40, -20, 20),
        "pt_genmet": (500, 0, 5000),
        "phi_genmet": (default_nbins, -3.14, 3.14),
        "m_l": (default_nbins, 0, 2),
        "m_n": (default_nbins, 0, 2),
        "mt": (500, 0, 10000),
        "mt_met": (500, 0, 10000),
        "m_inv": (500, 0, 10000),
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
    tfile = ROOT.TFile("histograms_t.root", "RECREATE")
    variables = ranges.keys()

    # Loop through skimmed datasets and produce histograms of variables
    for name, label in [
            ("120-200"  , "tnu_120-200"),
            ("200-400"  , "tnu_200-400"),
            #("400-800"  , "tnu_400-800"),
            ("800-1500"  , "tnu_800-1500"),
            ("1500-2500" , "tnu_1500-2500"),
            ("2500-4000" , "tnu_2500-4000"),
            ("4000-6000" , "tnu_4000-6000"),
            ("6000-inf"  , "tnu_6000-inf"),
        ]:
        print(">>> Process skimmed sample {} for process {}".format(name, label))

        # Load skimmed dataset
        df = ROOT.ROOT.RDataFrame("Events", name + "_Skim_t.root")

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
