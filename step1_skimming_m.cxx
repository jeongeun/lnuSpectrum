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
"120to200",//"200to400",
"400to800",
"800to1500",
"1500to2500",
"2500to4000",
"4000to6000",
"6000"
};

const float integratedLuminosity = 1;//27.0072 * 1000.0; // Run2022EFG (PostEE)
std::map<std::string, float> eventWeights = {
     {"120to200",   471.507     / 123093091.0 * integratedLuminosity},//     {"200to400",   74.0872     / 1460015.0 * integratedLuminosity},
     {"400to800",  7.78643    / 3598198.0 * integratedLuminosity},
     {"800to1500", 0.690298   / 3295801.0 * integratedLuminosity},
     {"1500to2500", 0.0482412  / 3302117.0 * integratedLuminosity},
     {"2500to4000", 0.00298857  / 3351427.0 * integratedLuminosity},
     {"4000to6000", 0.00011405  / 3222531.0 * integratedLuminosity},
     {"6000", 0.0000046658  / 3668638.0 * integratedLuminosity}
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
auto FindENuPair(T &df) {
using namespace ROOT::VecOps;
using floats = RVec<float>;
using ints = RVec<int>;
  auto LEPorigins = [](int nGenPart, ints GenPart_pdgId, ints GenPart_status, floats GenPart_pt)
      {
        ints out;
        int mu_idx = -1;
        int nu_idx  = -1;
        for( int i = 0; i < nGenPart; i++){
            //bool isLast = (GenPart_statusFlags[i] & (int) std::pow(2, 13)) != 0 ;// 13 = isLastCopy
            //bool isPrompt = (GenPart_statusFlags[i] & 1) != 0; //  0 = isPrompt
            bool isLastMuon = ( abs(GenPart_pdgId[i]) == 13  && GenPart_status[i] == 1); // && isPrompt);    
            if(isLastMuon){ mu_idx = i ; break; }
        }

        for( int j = 0; j < nGenPart; j++){
            //bool isLast = (GenPart_statusFlags[i] & (int) std::pow(2, 13)) != 0 ;// 13 = isLastCopy
            //bool isPrompt = (GenPart_statusFlags[i] & 1) != 0; //  0 = isPrompt
            bool isLastNeutrino = ( abs(GenPart_pdgId[j]) == 14  && GenPart_status[j] == 1 );     
            if(isLastNeutrino){ nu_idx  = j ; break; }
        }

        //if( mu_idx != -1 && nu_idx != -1 && (GenPart_pdgId[mu_idx] * GenPart_pdgId[nu_idx]) < 0){
        //    cout << "*** e " << GenPart_pdgId[mu_idx] << " nu " << GenPart_pdgId[nu_idx] 
        //         << " ept " << GenPart_pt[mu_idx] << ", npt " << GenPart_pt[nu_idx] << endl;
        //}
        out.emplace_back(mu_idx);
        out.emplace_back(nu_idx);

        return out;
    };

    return df.Define("LEPorigins", LEPorigins, 
             {"nGenPart","GenPart_pdgId","GenPart_status","GenPart_pt"})
             .Define("mu_idx"   , "LEPorigins[0]")
             .Define("nu_idx"    , "LEPorigins[1]")
             .Define("pid1"      , "GenPart_pdgId[mu_idx]")
             .Define("pid2"      , "GenPart_pdgId[nu_idx]")
             .Filter("mu_idx != -1", "Valid ele in selected idx")
             .Filter("nu_idx  != -1", "Valid nu in selected idx")
             .Filter("pid1 * pid2 < 0", "opposite charge pid1 * pid2 < 0");
}

/*
 *  * Declare all variables which we want to study in the analysis
*/
template <typename T>
auto DeclareVariables(T &df) {
    auto add_p4 = [](float pt, float eta, float phi, float mass)
    {
        return ROOT::Math::PtEtaPhiMVector(pt, eta, phi, mass);
    };
    auto compute_mt = [](float pt_1, float phi_1, float pt_met, float phi_met)
    {
        const auto dphi = Helper::DeltaPhi(phi_1, phi_met);
        return std::sqrt(2.0 * pt_1 * pt_met * (1.0 - std::cos(dphi)));
    };

    return df.Define("pt_l"      , "GenPart_pt[mu_idx]")
             .Define("eta_l"     , "GenPart_eta[mu_idx]")
             .Define("phi_l"     , "GenPart_phi[mu_idx]")
             .Define("mass_l"    , "GenPart_mass[mu_idx]")
             .Define("pdgId_l"   , "GenPart_pdgId[mu_idx]")
             .Define("pt_n"      , "GenPart_pt[nu_idx]")
             .Define("eta_n"     , "GenPart_eta[nu_idx]")
             .Define("phi_n"     , "GenPart_phi[nu_idx]")
             .Define("mass_n"    , "GenPart_mass[nu_idx]")
             .Define("pdgId_n"   , "GenPart_pdgId[nu_idx]")
             .Define("pt_genmet" , "GenMET_pt")
             .Define("phi_genmet", "GenMET_phi")
             .Define("p4_l"      , add_p4, {"pt_l", "eta_l", "phi_l", "mass_l"})
             .Define("p4_n"      , add_p4, {"pt_n", "eta_n", "phi_n", "mass_n"})
             .Define("p4"        , "p4_l + p4_n")
             .Define("m_inv"     , "float(p4.M())")
             .Define("scalePDF"  , "Generator_scalePDF")
             .Define("mt"        , compute_mt, {"pt_l", "phi_l", "pt_n", "phi_n"})
             .Define("mt_met"    , compute_mt, {"pt_l", "phi_l", "pt_genmet", "phi_genmet"});
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
    "pt_l"      ,"eta_l"     ,"phi_l"     ,"mass_l"    ,"pdgId_l"   ,
    "pt_n"      ,"eta_n"     ,"phi_n"     ,"mass_n"    ,"pdgId_n"   ,
    "pt_genmet" ,"phi_genmet","p4_l"      ,"p4_n"      ,"p4"        ,
    "m_inv"     ,"scalePDF"  ,"mt"        ,"mt_met"    ,"weight"
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

        auto df2 = FindENuPair(df);
        auto df3 = DeclareVariables(df2);
        auto df4 = AddEventWeight(df3, sample);

        auto dfFinal = df4;
        auto report = dfFinal.Report();
        dfFinal.Snapshot("Events", sample + "_Skim_m_mgmlm.root", finalVariables);
        time.Stop();

        report->Print();
        time.Print();
    }
}
