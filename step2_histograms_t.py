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
#    "lhe_pt_l"  ,"lhe_eta_l" ,"lhe_phi_l" ,"lhe_m_l"    ,"lhe_pdgId_l"   ,
#    "lhe_pt_n"  ,"lhe_eta_n" ,"lhe_phi_n" ,"lhe_m_n"    ,"lhe_pdgId_n"   ,
#    "pt_genmet" ,"phi_genmet","lhe_p4_l"  ,"lhe_p4_n"   ,"lhe_p4"        ,
#    "lhe_m_inv" ,"scalePDF"  ,"lhe_mt"    ,"lhe_mt_met" ,"weight" 

ranges = {
        "lhe_pt_l": (100, 0, 5000),
        "lhe_pt_n": (100, 0, 5000),
        "lhe_eta_l": (default_nbins, -3.0, 3.0),
        "lhe_eta_n": (default_nbins, -3.0, 3.0),
        "lhe_phi_l": (default_nbins, -3.14, 3.14),
        "lhe_phi_n": (default_nbins, -3.14, 3.14),
        "lhe_pdgId_l": (40, -20, 20),
        "lhe_pdgId_n": (40, -20, 20),
        "pt_genmet": (100, 0, 5000),
        "phi_genmet": (default_nbins, -3.14, 3.14),
        "lhe_m_l": (default_nbins, 0, 2),
        "lhe_m_n": (default_nbins, 0, 2),
        "lhe_mt": (200, 0, 8000),
        "lhe_mt_met": (200, 0, 8000),
        "lhe_m_inv": (200, 0, 8000),
        "scalePDF": (200, 0, 8000),
       } 


# Book a histogram for a specific variable
def bookHistogram(df, variable, range_):
    return df.Histo1D(ROOT.ROOT.RDF.TH1DModel(variable, variable, range_[0], range_[1], range_[2]), variable, "weight")


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
    #poolSize = ROOT.ROOT.GetImplicitMTPoolSize()
    #print("Pool size: {}".format(poolSize))

    # Create output file
    tfile = ROOT.TFile("histograms_t_mgmlm.root", "RECREATE")
    variables = ranges.keys()

    # Loop through skimmed datasets and produce histograms of variables
    for name, label in [
            #("120to200"   , "tnu_120to200"   ),
            #("200to400"   , "tnu_200to400"   ),
            ("400to800"   , "tnu_400to800"   ),
            ("800to1500"  , "tnu_800to1500"  ),
            ("1500to2500" , "tnu_1500to2500" ),
            ("2500to4000" , "tnu_2500to4000" ),
            ("4000to6000" , "tnu_4000to6000" ),
            ("6000"       , "tnu_6000"       ),
        ]:
        print(">>> Process skimmed sample {} for process {}".format(name, label))

        # Load skimmed dataset
        df = ROOT.ROOT.RDataFrame("Events", name + "_Skim_t_mgmlm.root")

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
