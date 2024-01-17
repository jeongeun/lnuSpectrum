#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
#include "TStopwatch.h"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include "ROOT/RDF/RInterface.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "Math/Vector4Dfwd.h"
#include "TStyle.h"
#include "TLorentzVector.h"

using namespace std;
using namespace ROOT::VecOps;
using floats = RVec<float>;
using bools = RVec<bool>;
using ints = RVec<int>;
using FourVector = ROOT::Math::PtEtaPhiMVector;
using FourVectorVec = std::vector<FourVector>;

//GenPart_statusFlags : //| UShort_t gen status flags stored bitwise, //bits are: 
//0 : isPrompt, 
//1 : isDecayedLeptonHadron, 
//2 : isTauDecayProduct, 
//3 : isPromptTauDecayProduct, 
//4 : isDirectTauDecayProduct, 
//5 : isDirectPromptTauDecayProduct, 
//6 : isDirectHadronDecayProduct, 
//7 : isHardProcess, 
//8 : fromHardProcess, 
//9 : isHardProcessTauDecayProduct, 
//10 : isDirectHardProcessTauDecayProduct, 
//11 : fromHardProcessBeforeFSR, 
//12 : isFirstCopy, 
//13 : isLastCopy, 
//14 : isLastCopyBeforeFSR,                           

/*
 *  * Base path to local filesystem or to EOS containing the datasets
 *   */
const std::string samplesBasePath = "/pnfs/knu.ac.kr/data/cms/store/user/jelee/WtoLNu-4Jets_mgmlm/EE/";
/*
 *  * Names of the datasets to be found in the base path and processed for the analysis
*/  
const std::vector<std::string> sampleNames = {
//"120to200",//"200to400",
"400to800",
"800to1500",
"1500to2500",
"2500to4000",
"4000to6000",
"6000"
};

const float integratedLuminosity = 1;//27.0072 * 1000.0; // Run2022EFG (PostEE)
std::map<std::string, float> eventWeights = {
     //{"120to200",   1.672e+02     / 123093091.0 * integratedLuminosity},//     {"200to400",  ??    / 1460015.0 * integratedLuminosity},
     {"400to800",   1.600e+00  / 3598198.0 * integratedLuminosity},
     {"800to1500",  1.091e-01  / 3295801.0 * integratedLuminosity},
     {"1500to2500", 6.536e-03  / 3302117.0 * integratedLuminosity},
     {"2500to4000", 3.484e-04  / 3351427.0 * integratedLuminosity},
     {"4000to6000", 1.077e-05  / 3222531.0 * integratedLuminosity},
     {"6000",       4.209e-07  / 3668638.0 * integratedLuminosity}
};
/*
 * Helper function to compute the difference in the azimuth coordinate taking
 * the boundary conditions at 2 * pi into account.
*/
namespace Helper {
template <typename T>
float DeltaPhi(T v1, T v2, const T c = M_PI)
{
    auto r = std::fmod(v2 - v1, 2.0 * c);
    if (r < -c) {
        r += 2.0 * c;
    }
    else if (r > c) {
        r -= 2.0 * c;
    }
    return r;
}
}

template <typename T>
auto FindGoodGenTaus(T &df) {
    return df.Define("goodGenTaus", "abs(LHEPart_pdgId) == 15 && LHEPart_status == 1");
}
template <typename T>
auto FindGoodGenNeutrinos(T &df) {
    return df.Define("goodGenNeutrinos", "abs(LHEPart_pdgId) == 16 && LHEPart_status == 1");
}

/*
 *  * Reduce the dataset to the interesting events containing at least one interesting
 *   * muon and tau candidate.
 *   */
template <typename T>
auto FilterGoodGenEvents(T &df) {
    return df.Filter("Sum(goodGenTaus) > 0", "Event has good genTaus")
             .Filter("Sum(goodGenNeutrinos) > 0", "Event has good genNeutrinos");
}
template <typename T>
auto DeclareLHEVariables(T &df) {
    auto add_p4 = [](float pt, float eta, float phi, float mass)
    {
        return ROOT::Math::PtEtaPhiMVector(pt, eta, phi, mass);
    };
    auto compute_mt = [](float pt_1, float phi_1, float pt_met, float phi_met)
    {
        const auto dphi = Helper::DeltaPhi(phi_1, phi_met);
        return std::sqrt(2.0 * pt_1 * pt_met * (1.0 - std::cos(dphi)));
    };

    return df.Define("lhe_pt_l"    , "LHEPart_pt[2]")
             .Define("lhe_eta_l"   , "LHEPart_eta[2]")
             .Define("lhe_phi_l"   , "LHEPart_phi[2]")
             .Define("lhe_m_l"     , "LHEPart_mass[2]")
             .Define("lhe_pdgId_l" , "LHEPart_pdgId[2]")
             .Define("lhe_pt_n"    , "LHEPart_pt[3]")
             .Define("lhe_eta_n"   , "LHEPart_eta[3]")
             .Define("lhe_phi_n"   , "LHEPart_phi[3]")
             .Define("lhe_m_n"     , "LHEPart_mass[3]")
             .Define("lhe_pdgId_n" , "LHEPart_pdgId[3]")
             .Define("pt_genmet"   , "GenMET_pt")
             .Define("phi_genmet"  , "GenMET_phi")
             .Define("lhe_p4_l"    , add_p4, {"lhe_pt_l", "lhe_eta_l", "lhe_phi_l", "lhe_m_l"})
             .Define("lhe_p4_n"    , add_p4, {"lhe_pt_n", "lhe_eta_n", "lhe_phi_n", "lhe_m_n"})
             .Define("lhe_p4"      , "lhe_p4_l + lhe_p4_n")
             .Define("lhe_m_inv"   , "float(lhe_p4.M())")
             .Define("lhe_mt"      , compute_mt, {"lhe_pt_l", "lhe_phi_l", "lhe_pt_n", "lhe_phi_n"})
             .Define("lhe_mt_met"  , compute_mt, {"lhe_pt_l", "lhe_phi_l", "pt_genmet", "phi_genmet"})
             .Define("scalePDF"    , "Generator_scalePDF");
}

/*
 *  * Add the event weight to the dataset as the column "weight"
*/
template <typename T>
auto AddEventWeight(T &df, const std::string& sample) {
    const auto weight = eventWeights[sample];
    return df.Define("weight", [weight](){ return weight; });
}

/*
 * Declare all variables which shall end up in the final reduced dataset
 */
const std::vector<std::string> finalVariables = {
    "lhe_pt_l"  ,"lhe_eta_l" ,"lhe_phi_l" ,"lhe_m_l"    ,"lhe_pdgId_l"   ,
    "lhe_pt_n"  ,"lhe_eta_n" ,"lhe_phi_n" ,"lhe_m_n"    ,"lhe_pdgId_n"   ,
    "pt_genmet" ,"phi_genmet","lhe_p4_l"  ,"lhe_p4_n"   ,"lhe_p4"        ,
    "lhe_m_inv" ,"scalePDF"  ,"lhe_mt"    ,"lhe_mt_met" ,"weight" 
};

/*
 * Main function of the skimming step of the analysis
 * The function loops over all required samples, reduces the content to the
 * interesting events and writes them to new files.
 */
int main() {
    ROOT::EnableImplicitMT();
    //const auto poolSize = ROOT::GetImplicitMTPoolSize();
    //std::cout << "Pool size: " << poolSize << std::endl;

    for (const auto &sample : sampleNames) {
        std::cout << ">>> Process sample " << sample << ":" << std::endl;
        TStopwatch time;
        time.Start();

        ROOT::RDataFrame df("Events", samplesBasePath + sample + "/*.root");
        std::cout << "RDataFrame " << std::endl;
        std::cout << "Number of events: " << *df.Count() << std::endl;

        auto df1 = FindGoodGenTaus(df);
        auto df2 = FindGoodGenNeutrinos(df1);
        auto df3 = FilterGoodGenEvents(df2);
        auto df4 = DeclareLHEVariables(df3);
        //auto df5 = FindENuPair(df3);
        //auto df6 = DeclareGENVariables(df5);
        auto df5 = AddEventWeight(df4, sample);
        auto dfFinal = df5;
        auto report = dfFinal.Report();
        dfFinal.Snapshot("Events", sample + "_Skim_t_mgmlm.root", finalVariables);
        time.Stop();

        report->Print();
        time.Print();
    }
}
