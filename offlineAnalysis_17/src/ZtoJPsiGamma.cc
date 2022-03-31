#define ZtoJPsiGamma_cxx
#include "ZtoJPsiGamma.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include "TDirectory.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TString.h>
#include <assert.h>
#include <TFile.h>
#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1K.h"
#include "TObject.h"
#include "TH1D.h"
#include "TSystem.h"
#include "configana.h"

using namespace std;


//*******************************
/*
namespace ZtoJPsiGamma::AnaUtil {
  void tokenize(const string& str, vector<string>& tokens, const string& delimiters) {
    // Skip delimiters at beginning.                                                                                                         
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);

-    // Find first "non-delimiter".                                                                                                            
    string::size_type pos = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)  {
      // Found a token, add it to the vector.                                                                                                 
      tokens.push_back(str.substr(lastPos, pos - lastPos));

      // Skip delimiters.  Note the "not_of"                                                                                                  
      lastPos = str.find_first_not_of(delimiters, pos);

      // Find next "non-delimiter"                                                                                                            
      pos = str.find_first_of(delimiters, lastPos);
    }
  }

  double cutValue(const map<string, double>& m, const string& mkey) {
    if (m.find(mkey) == m.end()) {
      cerr << ">>> key: " << mkey << " not found in the map!" << endl;
      for (auto const& el: m)
        cerr << el.first << ": " << setw(7) << el.second << endl;
      return -999;
    }
    //assert(m.find(cname) != m.end());                                                                                                       
    return m.find(mkey)->second;
  }
  const map<string, double>& cutMap(const map<string, map<string, double>>& hmap, const string& mkey) {
    assert(hmap.find(mkey) != hmap.end());
    return hmap.find(mkey)->second;
  }
  void buildList(const vector<string>& tokens, vector<string>& list) {
    for (vector<string>::const_iterator it  = tokens.begin()+1;
	 it != tokens.end(); ++it) {
      list.push_back(*it);
    }
  }
  void buildMap(const vector<string>& tokens, map<string, int>& hmap) {
    string key = tokens.at(1) + "-" + tokens.at(2) + "-" + tokens.at(3);
    hmap.insert({key, 1});
  }
  void buildMap(const vector<string>& tokens, unordered_map<string, int>& hmap) {
    string key = tokens.at(1) + "-" + tokens.at(2) + "-" + tokens.at(3);
    hmap.insert({key, 1});
  }
  void storeCuts(const vector<string>& tokens, map<string, map<string, double>>& hmap) {
    const string& key = tokens.at(0);
    auto pos = hmap.find(key);
    if (pos == hmap.end()) {
      map<string, double> m;
      for (auto it = tokens.begin()+1; it != tokens.end(); ++it) {
        // Split the line into words                                                                                                          
        vector<string> cutstr;
        tokenize(*it, cutstr, "=");
        if (cutstr.size() < 2) continue;
        m.insert({cutstr.at(0), atof(cutstr.at(1).c_str())});
      }
      if (!m.empty()) hmap.insert({key, m});
    }
  }
  void showCuts(const map<string, map<string, double>>& hmap, ostream& os) {
    for (auto const& el: hmap) {
      os << ">>> " << el.first << endl;
      auto const& cutm = el.second;
      os << std::setprecision(2);
      for (auto const& il: cutm)
        os << setw(16) << il.first << ": "
           << setw(7)  << il.second << endl;
    }
  }
}
*/
//*********************************************


void ZtoJPsiGamma::tokenize(const string& str, vector<string>& tokens, const string& delimiters) {
  // Skip delimiters at beginning.                                                                                                         
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  
  // Find first "non-delimiter".                                                                                                            
  string::size_type pos = str.find_first_of(delimiters, lastPos);
  
  while (string::npos != pos || string::npos != lastPos)  {
    // Found a token, add it to the vector.                                                                                                 
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    
    // Skip delimiters.  Note the "not_of"                                                                                                  
    lastPos = str.find_first_not_of(delimiters, pos);
    
    // Find next "non-delimiter"                                                                                                            
    pos = str.find_first_of(delimiters, lastPos);
  }
}

void ZtoJPsiGamma::buildList(const vector<string>& tokens, vector<string>& list) {
  for (vector<string>::const_iterator it  = tokens.begin()+1; it != tokens.end(); ++it){
    list.push_back(*it);
  }
}

int ZtoJPsiGamma::setInputFile(const string& fname){
  auto found = fname.find("root:");
  if (found == string::npos && gSystem->AccessPathName(fname.c_str())){
    cerr << ">>> Warning: File <<" << fname << ">> was not found!!" << endl;
    return static_cast<int>(chain->GetEntries());
  }
  chain->Add(fname.c_str(), -1);
  return static_cast<int>(chain->GetEntries());

}

void ZtoJPsiGamma::storeCuts(const std::vector<string>& tokens, map<string, map<string, double>>& hmap) {
  const string& key = tokens.at(0);
  auto pos = hmap.find(key);
  if (pos == hmap.end()) {
    map<string, double> m;
    for (auto it = tokens.begin()+1; it != tokens.end(); ++it) {
      // Split the line into words
      vector<string> cutstr;
      tokenize(*it, cutstr, "=");
      if (cutstr.size() < 2) continue;
      m.insert({cutstr.at(0), atof(cutstr.at(1).c_str())});
    }
    if (!m.empty()) hmap.insert({key, m});
  }
}

TH1* ZtoJPsiGamma::getHist1D(const char* hname) {
  TObject *obj = gDirectory->GetList()->FindObject(hname);
  if (obj == nullptr) {
    cerr << "**** getHist1D: Histogram for <" << hname
	 << "> not found! ("
	 << __FILE__ << ":" << __LINE__ << ")"
	 << endl;
    return nullptr;
  }
  TH1* h = nullptr;
  if (obj->InheritsFrom("TH1D"))
    h = dynamic_cast<TH1D*>(obj);
  else if (obj->InheritsFrom("TH1C"))
    h = dynamic_cast<TH1C*>(obj);
  else if (obj->InheritsFrom("TH1K"))
    h = dynamic_cast<TH1K*>(obj);
  else if (obj->InheritsFrom("TH1S"))
    h = dynamic_cast<TH1S*>(obj);
  else if (obj->InheritsFrom("TH1I"))
    h = dynamic_cast<TH1I*>(obj);
  else
    h = dynamic_cast<TH1F*>(obj);
  if (h == nullptr) {
    cerr << "**** getHist1D: <" << hname
	 << "> may not be a 1D Histogram! ("
	 << __FILE__ << ":" << __LINE__ << ")"
	 << endl;
  }
  return h;
}

TH1* ZtoJPsiGamma::getHist1D(const string& hname) {
  return getHist1D(hname.c_str());
}



TH2* ZtoJPsiGamma::getHist2D(const char* hname) {
  TObject *obj = gDirectory->GetList()->FindObject(hname);
  if (obj == nullptr) {
    cerr << "**** getHist2D: Histogram for <" << hname
	 << "> not found! ("
	 << __FILE__ << ":" << __LINE__ << ")"
	 << endl;
    return nullptr;
  }
  TH2* h = nullptr;
  if (obj->InheritsFrom("TH2D"))
    h = dynamic_cast<TH2D*>(obj);
  else if (obj->InheritsFrom("TH2C"))
    h = dynamic_cast<TH2C*>(obj);
  else if (obj->InheritsFrom("TH2S"))
    h = dynamic_cast<TH2S*>(obj);
  else if (obj->InheritsFrom("TH2I"))
    h = dynamic_cast<TH2I*>(obj);
  else
    h = dynamic_cast<TH2F*>(obj);
  if (h == nullptr) {
    cerr << "**** getHist2D: <" << hname
	 << "> may not be a 2D Histogram! ("
	 << __FILE__ << ":" << __LINE__ << ")"
	 << endl;
  }
  return h;
}

TH2* getHist2D(const string& hname) {
  return getHist2D(hname.c_str());
}

double ZtoJPsiGamma::cutValue(const map<string, double>& m, const string& mkey) {
  if (m.find(mkey) == m.end()) {
    cerr << ">>> key: " << mkey << " not found in the map!" << endl;
    cout << ">>> key: " << mkey << " not found in the map!" << endl;
    for (auto const& el: m)
      cerr << el.first << ": " << setw(7) << el.second << endl;
    return -999;
  }
  //assert(m.find(cname) != m.end());                                                                                                       
  return m.find(mkey)->second;
}
const map<string, double>& ZtoJPsiGamma::cutMap(const map<string, map<string, double>>& hmap, const string& mkey) {
  assert(hmap.find(mkey) != hmap.end());
  return hmap.find(mkey)->second;
}

const std::map<std::string, double>& ZtoJPsiGamma::lumiWtMap() {return ZtoJPsiGamma::cutMap(hmap_, "lumiWtList");}
const std::map<std::string, double>& ZtoJPsiGamma::ratioCutMap() {return ZtoJPsiGamma::cutMap(hmap_, "ratioCutList");}
const std::map<std::string, double>& ZtoJPsiGamma::SRCutMap() {return ZtoJPsiGamma::cutMap(hmap_, "SRCutList");}

/*
double ZtoJPsiGamma::dopuweight(int& num_pu_vt) const {
  
  const char* fname = gSystem->ExpandPathName(pufilename);
  if (gSystem->AccessPathName(fname)) {
    cerr << ">>> Warning: File <<" << pufilename << ">> not found!!" << endl;
    return false;
  }
  
  
  TFile file(fname);
  //TH1F *h = dynamic_cast<TH1F*>(file.Get(puhistname.c_str()));
  //TFile* fPU = new TFile("mcPileup2017.root");
  TH1D *h = dynamic_cast<TH1D*>(fPU->Get(puhistname));
  //TH1F *h = dynamic_cast<TH1F*>(fPU->Get(puhistname));
  if (!h) {
    cerr << ">>> Warning: Histogram <<" << puhistname << ">> not found!!" << endl;
    return false;
  }

  //if (num_pu_vt < 0) continue;                                                                                                             
  h->Scale(1/h->Integral());
  int nx = h->GetXaxis()->FindBin(num_pu_vt);
  std::cout.precision(2);
  //for (int i = 0; i < nx; ++i) {                                                                                                            
  double wt = h->GetBinContent(nx);
  //if (verbose)
    std::cout << setw(8) << nx
      << setw(8) << wt
      << std::endl;
  
  //puWtList_.push_back(wt);                                                                                                             
  //}                                                                                                                                       
  return wt;

  
}
*/
/*
double ZtoJPsiGamma::cutValue(const map<string, double>& m, const string& mkey) {
  if (m.find(mkey) == m.end()) {
    cerr << ">>> key: " << mkey << " not found in the map!" << endl;
    for (auto& el: m)
      cerr << el.first << ": " << setw(7) << el.second << endl;
    return -999;
  }
  //assert(m.find(cname) != m.end());
    return m.find(mkey)->second;
}


map<string, double>& ZtoJPsiGamma::cutMap(const map<string, map<string, double>>& hmap, const string& mkey) {
  assert(hmap.find(mkey) != hmap.end());
  return hmap.find(mkey)->second;
}
*/
//double ZtoJPsiGamma::lumiWt(double evtWeightSum, int nent) const
double ZtoJPsiGamma::lumiWt(double evtWeightSum)
{
  //    double nevt = (evtWeightSum > -1) ? evtWeightSum : nent;
    double nevt = (evtWeightSum > -1) ? evtWeightSum : ZtoJPsiGamma::cutValue(ZtoJPsiGamma::lumiWtMap(), "nevents"); 
  //  if (!verbose)
  //  xsec = xsec_;
  
  /*  std::cout << "-- intLumi: " << intLumi
	    << " xsec: " << xsec
	    << " nevt: " << nevt << std::endl;*/
  
  std::cout << "-- intLumi: " << ZtoJPsiGamma::cutValue(ZtoJPsiGamma::lumiWtMap(), "intLumi")
	    << " xsec: " << ZtoJPsiGamma::cutValue(ZtoJPsiGamma::lumiWtMap(), "xsec")
	    << " nevt: " << nevt << std::endl;
  
  //    return (intLumi * xsec / nevt);
  return (ZtoJPsiGamma::cutValue(ZtoJPsiGamma::lumiWtMap(), "intLumi") * ZtoJPsiGamma::cutValue(ZtoJPsiGamma::lumiWtMap(), "xsec") / nevt);
}

// taken from HZZ4lUtil to show cutflow efficiency
void ZtoJPsiGamma::showEfficiency(const string& hname,const std::vector<std::string>& slist,const string& header,const string& tag,std::ostream& os)
{
  os << ">>> " << header << " Efficiency" << endl;
  TH1 *h = ZtoJPsiGamma::getHist1D(hname);
  if (h != nullptr) {
    os << setw(64) << "CutFlow"
       << setw(20) << tag
       << setw(20) << "AbsEff"
       << setw(20) << "RelEff"
       << endl;
    os.precision(3);
    int nbins = h->GetNbinsX();
    for (int i = 1; i <= nbins; ++i) {
      double cont  = h->GetBinContent(1);
      double conti = h->GetBinContent(i);
      double contj = h->GetBinContent(i-1);
      os << setw(64) << slist[i-1]
	 << setprecision(8)
	 << setw(20) << conti
	 << setprecision(5)
	 << setw(20) << ((conti > 0) ? conti/cont : 0.0)
	 << setw(20) << ( i == 1 ? 1.0 :(contj > 0) ? conti/contj : 0.0)
	 << endl;
    }
  }
}

void ZtoJPsiGamma::objectEfficiency() {
  // Object selection Efficiency
  // outfile->cd();
  // outfile->cd("ObjectSelection");
  
  // Muon  
  vector<string> muLabels {
    "All",
      "looseid",
      "tightid"
      };
  ZtoJPsiGamma::showEfficiency("muCutFlow", muLabels, "Muon Selection","Muons");

  // Photon
  vector<string> phoLabels {
    "All",
      "mvaid",
      "pt",
      "ScEta"
      };
  ZtoJPsiGamma::showEfficiency("photonCutFlow", phoLabels, "Photon Selection","Photons");
}

void ZtoJPsiGamma::eventEfficiency(){
  //outfile->cd();
  //outfile->cd("ZtoJPsiGamma");
  std::vector<std::string> evlabels {
    "Events processed",
      //      "Gen level filter : Minimum 2 gen Muons & Minimum 1 gen photon",
      //"Gen level filter after detector cuts, eta =< 2.5",
      "Events with > 0 good vertex",
      "Events passing trigger",
      "Atleast 2 muons after muon selection",
      "Relative Isolation of leading mu",
      "leading mu pt",
      "subleading mu pt",
      "Atleast 1 photon after photon selection",
      "deltaR between mu1 & pho",
      "deltaR between mu2 & pho",
      "deltaR between dimu & pho",
      "deltaphi between dimu & pho",
      "JPsi mass cut",
      "Z mass cut",
      "ratio cut for pt of dimu",
      "ratio cut for Et of pho"
      };
  
  ZtoJPsiGamma::showEfficiency("evtCutFlow",evlabels, "Event Selection","Eventflow");

  double lumiFac = 1.0;
  lumiFac = lumiWt(evtWeightSum_);
  cout << endl
       << "evtWeightSum: " << setw(10) << setprecision(0) << evtWeightSum_ << endl
       << "      lumiWt: " << setw(10) << setprecision(5) << lumiFac
       << endl;

  //  ZtoJPsiGamma::scaleHistogram("evtCutFlowWt", lumiFac);

  TH1 *h = ZtoJPsiGamma::getHist1D("evtCutFlowWt");
  if (h != nullptr) {
    int nbins = h->GetNbinsX();
    for (int i = 1; i <= nbins; ++i) {
      double cont = h->GetBinContent(i) * lumiFac;
      double err  = h->GetBinError(i) * lumiFac;
      h->SetBinContent(i, cont);
      h->SetBinError(i, err);
    }
  }

  ZtoJPsiGamma::showEfficiency("evtCutFlowWt", evlabels, "Event Selection (Weighted)", "EventFlow");

      std::vector<std::string> catlabels {
      "total final yield",
	"cat1",
	"cat2",
	"cat3"
	//"cat4"
      };

        ZtoJPsiGamma::showEfficiency("nPhoCat",catlabels, "category of events", "category");

  double lumiFac_cat = 1.0;
  lumiFac_cat = lumiWt(evtWeightSum_);
  cout << endl
       << "evtWeightSum: " << setw(10) << setprecision(0) << evtWeightSum_ << endl
       << "      lumiWt: " << setw(10) << setprecision(5) << lumiFac_cat
       << endl;
  
  TH1 *hcat = ZtoJPsiGamma::getHist1D("nPhoCatWt");
  if (hcat != nullptr) {
    int nbins_cat = hcat->GetNbinsX();
    for (int j = 1; j <= nbins_cat; ++j) {
      double cont_cat = hcat->GetBinContent(j) * lumiFac_cat;
      double err_cat  = hcat->GetBinError(j) * lumiFac_cat;
      hcat->SetBinContent(j, cont_cat);
      hcat->SetBinError(j, err_cat);
    }
  }

  ZtoJPsiGamma::showEfficiency("nPhoCatWt",catlabels, "category of weighted events", "Category");

}
      
void ZtoJPsiGamma::showCategoryYield(){
  // just for test, please ignore

}

// Helicity angle between J/Psi and mu+ 
float ZtoJPsiGamma::costheta1(TLorentzVector particle_orig, TLorentzVector parent_orig){
  TLorentzVector particle = particle_orig; // make a dummy
  TLorentzVector parent = parent_orig;
  TVector3 boosttoparent = -(parent.BoostVector());
  particle.Boost(boosttoparent);
  parent.Boost(boosttoparent);
 
  TVector3 particle3 = particle.Vect();
  TVector3 parent3 = parent.Vect();
  Float_t numerator = particle3.Dot(parent3);
  Float_t denominator = (particle3.Mag())*(parent3.Mag());
  Float_t temp = numerator/denominator;

  return temp;
}

float ZtoJPsiGamma::costhetastar(TLorentzVector particle_orig, TLorentzVector parent_orig){
  TLorentzVector particle = particle_orig; // make a dummy
  TLorentzVector parent = parent_orig;
  TVector3 boosttoparent = -(parent.BoostVector());
  particle.Boost(boosttoparent);
  parent.Boost(boosttoparent);
 
  //  TVector3 particle3 = particle.Vect();
  //  TVector3 parent3 = parent.Vect();
  Float_t temp = particle.CosTheta();

  return temp;
}

void ZtoJPsiGamma::bookHistograms() {
  //TString of =  "ZtoJPsiGamma_hist.root";
  //outfile = new TFile(of,"RECREATE");
  outfile->cd();
  outfile->mkdir("ObjectSelection");
  outfile->cd("ObjectSelection");
  new TH1D("muCutFlow", "Muon Cut Flow", 3, -0.5, 2.5);
  new TH1D("photonCutFlow", "Photon Cut Flow", 4, -0.5, 3.5);
  
  new TH1F("muPt", "Muon p_{T}", 200, 0., 100.);//mu loose id
  new TH1F("muEta", "Muon #eta", 200, -4., 4.);
  new TH1F("muDxy", "Muon d_{xy}", 50, -10., 10.);
  new TH1F("muDz", "Muon d_z", 50, -10., 10.);
  new TH1F("muSIP3d", "Muon sip3d", 50, -10., 10.);
  new TH1F("muSIP3dincm", "Muon sip3d in cm", 50, -10., 10.);
  new TH1F("muLooseid", "Muon looseid", 2, 0., 1.);
  new TH1F("muIspfcand", "Muon isPFcand", 2, 0., 1.);//mutightid
  new TH1F("muHighptid", "Muon isHighptid", 2, 0., 1.);
  new TH1F("muTightid", "Muon tightId", 2, 0., 1.);

  new TH1F("phoMvaidwp90", "Photon is Mvaidwp90", 2, 0., 1.);//photon
  new TH1F("phoPt", "Photon p_{T}", 400, 0., 200.);
  new TH1F("Photon_isScEtaEB", "Photon isScEtaEB", 2, 0., 1.);
  new TH1F("Photon_isScEtaEE", "Photon isScEtaEE", 2, 0., 1.);
  new TH1D("nPhoScEta", "#Photons in ScEta regions", 2, -0.5, 1.5);
  new TH1D("phoR9", "Photon R9", 40, 0.4, 1.0);
  new TH1D("phoR9_ScEB", "Photon R9 for ScEta in EB", 40, 0.4, 1.0);
  new TH1D("phoR9_ScEE", "Photon R9 for ScEta in EE", 40, 0.4, 1.0);

  outfile->cd();
  outfile->mkdir("ZtoJPsiGamma");
  outfile->cd("ZtoJPsiGamma");
  new TH1F("generatorweight", "MC generator weight", 3, -1.5, 1.5);
  new TH1F("puweight", "PU weight", 200, 0., 2.);
  new TH2F("pudist","Pu weight distribution", 80, 0., 80., 100., 0., 2. );
  /*
  new TH1F("genweight", "generator weight", 3, -1.5, 1.5);
  new TH1F("genMu1Pt","leading gen mu p_T ", 50, 0.,100.);
  new TH1F("genMu2Pt","subleading gen mu p_T ", 50,0., 100.);
  new TH1F("genPhoPt"," gen photon p_T ", 40, 0., 200.);
  new TH1F("genMu1Eta","leading gen mu eta ", 50, -2.5, 2.5);
  new TH1F("genMu2Eta","subleading gen mu eta ", 50, -2.5 , 2.5);
  new TH1F("genPhoEta"," gen photon eta ", 50, -2.5, 2.5);
  new TH1F("gendRMu1Mu2","dR between gen mu1 & mu2", 50, 0., 1.);
  new TH1F("gendRdiMuPho","dR between gen diMu & Pho", 50, 0., 5.);
  new TH1F("gendRMu1Pho","dR between gen diMu1 & Pho", 50, 0., 5.);
  new TH1F("gendRMu2Pho","dR between gen diMu2 & Pho", 50, 0., 5.);
  */
  //Preselection
  new TH1F("nPV", "No. of reconstructed primary vertices", 15, -0.5, 14.5);
  new TH1F("nPVgood", "No. of good reconstructed primary vertices, passing all the criteria", 15, -0.5, 14.5);

  //Helicity
  new TH1F("costheta1","Helicity angle between JPsi and mu",80, -1.0, 1.0);
  new TH1F("costheta1_rewt","Helicity angle between JPsi and mu (reweighted)",80, -1.0, 1.0);
  new TH1F("costhetas","Helicity angle between JPsi and beam",80, -1.0, 1.0);
  new TH1F("costhetas_rewt","Helicity angle between JPsi and beam (reweighted)",80, -1.0, 1.0);

  new TH1F("nLoosemuons", "No. of loose muons", 10, -0.5, 9.5);
  new TH1F("nTightmuons", "No. of tight muons", 10, -0.5, 9.5);
  new TH1F("nPhotons", "No. of photons", 5, -0.5, 4.5);

  new TH1F("muPfRelIso03all", "muPfRelIso03all", 100, 0., 1.);
  new TH1F("muPfRelIso03chg", "muPfRelIso03chg", 100, 0., 1.);
  new TH1F("muPfRelIso04all", "muPfRelIso04all", 100, 0., 1.);
  new TH1F("muPfRelIso03allByMuSubleadpt", "muPfRelIso03all", 80, 0., 2.);
  new TH1F("muLeadpt", "Leading Mu pt", 20, 0., 100.);
  new TH1F("muSubleadpt", "Subleading Mu pt", 20, 0., 100.);
  
  new TH1F("dRmu1pho", "dR between mu1 & photon", 100, 0., 5.);
  new TH1F("dRmu2pho", "dR between mu2 & photon", 100, 0., 5.);
  new TH1F("dRdimupho", "dR between diMuon & photon", 120, 0., 8.);
  new TH1F("dPhidimupho", "dPhi between diMuon & photon", 120, -4., 4.);
  new TH1F("ratiocutdimu", "ratio cut for diMuon", 100, 0., 1.);
  new TH1F("ratiocutpho", "ratio cut for photon", 100, 0., 1.);


  new TH1F("jpsimass_bmasscut","jpsimass before mass cut", 60, 2.8, 3.4);
  new TH1F("jpsimass_amasscut","jpsimass after mass cut", 60, 2.8, 3.4);
  new TH1F("jpsimass_afinalcut","jpsimass after final cut", 60, 2.8, 3.4);
  new TH1F("zmass_bmasscut","zmass before Z mass cut", 64, 40., 200.);
  new TH1F("zmass_amasscut","zmass after Z mass cut", 64, 40., 200.);
  new TH1F("zmass_afinalcut","zmass after final cut", 64, 40., 200.);

  new TH1F("vtxdist_x","DiMuon vertex position (X)", 40, -0.2, 0.2);
  new TH1F("vtxdist_y","DiMuon vertex position (Y)", 40, -0.2, 0.2);
  new TH1F("vtxdist_z","DiMuon vertex position (Z)", 40, -0.2, 0.2);
  
  new TH1D("evtCutFlow", "Event CutFlow", 16, -0.5, 15.5); // event cutflow
  new TH1D("evtCutFlowWt", "Event CutFlow(Weighted)", 16, -0.5, 15.5); // event cutflow weighted
  
  // Event category for Z signal
 
  new TH1D("nPhoCat", "Event Category", 4, -0.5, 3.5);
  new TH1D("nPhoCatWt", "Event Category weighted", 4, -0.5, 3.5);

  /*
  new TH1F("patdRmu1mu2_bHLT","dR between pat mu1 & mu2 before HLT", 50, 0., 1.);
  new TH1F("patMmu1mu2_bHLT","Invariant Mass of pat mu1 & mu2 before HLT", 30, 2.8, 3.4);
  new TH1F("patPTmu1mu2_bHLT","Pt of pat mu1 & mu2 before HLT", 45, 10., 100.);
  new TH1F("patdRmu1mu2_aHLT","dR between pat mu1 & mu2 after HLT", 50, 0., 1.);
  new TH1F("patMmu1mu2_aHLT","Invariant Mass of pat mu1 & mu2 after HLT", 30, 2.8, 3.4);
  new TH1F("patPTmu1mu2_aHLT","Pt of pat mu1 & mu2 after HLT", 45, 10., 100.);
  */
  // categorywise final plot 
  new TH1F("muLeadpt_cat1", "Leading Mu pt for cat1 ", 32, 0., 160.);
  new TH1F("muSubleadpt_cat1", "Subleading Mu pt for cat1", 32, 0., 160.);
  new TH1F("phoPt_cat1", "Leading Mu pt for cat1 ", 32, 0., 160.);
  new TH1F("muLeadeta_cat1","leading Mu eta for cat1 ", 25, -2.5, 2.5);
  new TH1F("muSubleadeta_cat1","Subleading Mu eta for cat1 ", 25, -2.5, 2.5);
  new TH1F("phoEta_cat1","leading photon eta for cat1 ", 25, -2.5, 2.5);
  new TH1F("dRmu1pho_cat1", "dR between mu1 & photon for cat1", 50, 0., 5.);
  new TH1F("dRmu2pho_cat1", "dR between mu2 & photon for cat1", 50, 0., 5.);
  new TH1F("dRdimupho_cat1", "dR between diMuon & photon for cat1", 80, 0., 8.);
  new TH1F("dPhidimupho_cat1", "dPhi between diMuon & photon for cat1", 80, -4., 4.);
  new TH1F("dRmu1mu2_cat1", "dR between mu1 & mu2 for cat1", 100, 0., 1.);
  new TH1F("jpsimass_cat1","jpsimass for cat1", 40, 2.9, 3.3);  
  new TH1F("diMuPt_cat1","DiMuon Pt for cat1", 32, 0. , 160.);  
  new TH1F("finalZPt_cat1","Pt for final dimuon+gamma for cat1", 25, 0., 100.);  
  new TH1F("zmass_cat1","Zmass for cat1", 64, 40., 200.);  

  new TH1F("muLeadpt_cat2", "Leading Mu pt for cat2 ", 32, 0., 160.);
  new TH1F("muSubleadpt_cat2", "Subleading Mu pt for cat2", 32, 0., 160.);
  new TH1F("phoPt_cat2", "Leading Mu pt for cat2 ", 32, 0., 160.);
  new TH1F("muLeadeta_cat2","leading Mu eta for cat2 ", 25, -2.5, 2.5);
  new TH1F("muSubleadeta_cat2","Subleading Mu eta for cat2 ", 25, -2.5, 2.5);
  new TH1F("phoEta_cat2","leading photon eta for cat2 ", 25, -2.5, 2.5);
  new TH1F("dRmu1pho_cat2", "dR between mu1 & photon for cat2", 50, 0., 5.);
  new TH1F("dRmu2pho_cat2", "dR between mu2 & photon for cat2", 50, 0., 5.);
  new TH1F("dRdimupho_cat2", "dR between diMuon & photon for cat2", 80, 0., 8.);
  new TH1F("dPhidimupho_cat2", "dPhi between diMuon & photon for cat2", 80, -4., 4.);
  new TH1F("dRmu1mu2_cat2", "dR between mu1 & mu2 for cat2", 100, 0., 1.);
  new TH1F("jpsimass_cat2","jpsimass for cat2", 40, 2.9, 3.3);
  new TH1F("diMuPt_cat2","DiMuon Pt for cat2", 32, 0. , 160.);
  new TH1F("finalZPt_cat2","Pt for final dimuon+gamma for cat2", 25, 0., 100.);
  new TH1F("zmass_cat2","Zmass for cat2", 64, 40., 200.);

  new TH1F("muLeadpt_cat3", "Leading Mu pt for cat3 ", 32, 0., 160.);
  new TH1F("muSubleadpt_cat3", "Subleading Mu pt for cat3", 32, 0., 160.);
  new TH1F("phoPt_cat3", "Leading Mu pt for cat3 ", 32, 0., 160.);
  new TH1F("muLeadeta_cat3","leading Mu eta for cat3 ", 25, -2.5, 2.5);
  new TH1F("muSubleadeta_cat3","Subleading Mu eta for cat3 ", 25, -2.5, 2.5);
  new TH1F("phoEta_cat3","leading photon eta for cat3 ", 25, -2.5, 2.5);
  new TH1F("dRmu1pho_cat3", "dR between mu1 & photon for cat3", 50, 0., 5.);
  new TH1F("dRmu2pho_cat3", "dR between mu2 & photon for cat3", 50, 0., 5.);
  new TH1F("dRdimupho_cat3", "dR between diMuon & photon for cat3", 80, 0., 8.);
  new TH1F("dPhidimupho_cat3", "dPhi between diMuon & photon for cat3", 80, -4., 4.);
  new TH1F("dRmu1mu2_cat3", "dR between mu1 & mu2 for cat3", 100, 0., 1.);
  new TH1F("jpsimass_cat3","jpsimass for cat3", 40, 2.9, 3.3);
  new TH1F("diMuPt_cat3","DiMuon Pt for cat3", 32, 0. , 160.);
  new TH1F("finalZPt_cat3","Pt for final dimuon+gamma for cat3", 25, 0., 100.);
  new TH1F("zmass_cat3","Zmass for cat3", 64, 40., 200.);
  

  // To select which cut to choose
  //  new TH1F("dRMu1Cut", "# of Events for different dR(mu1, #gamma)", )
    
  bookedHistograms = true;
  
}

void ZtoJPsiGamma::saveHistograms() {
  cout << "About to close the histogram file" << endl;
  if (bookedHistograms) {
    //    outfile->cd();
    //outfile->cd("ZtoJPsiGamma");
    outfile->Write();  
    //    outfile->Close();
  }
  cout << "histogram file not closed" << endl;

}
// ***********************************************


bool ZtoJPsiGamma::readJob(const std::string& jobFile, int& nFiles){
  
  // Open the file containing the datacards
  ifstream fin(jobFile.c_str(), std::ios::in);
  if (!fin) {
    cerr << "==> Input File: <<" << jobFile << ">> could not be opened " << endl;
    return false;
  }
  // eventFileList_.clear();

  static constexpr int BUF_SIZE = 256;
  char buf[BUF_SIZE];
  while (fin.getline(buf, BUF_SIZE, '\n')) {
    // Pops off the newline character
    string line(buf);
    if (line.empty() || line == "START") continue;

    // enable '#' and '//' style comments                           
    if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
    if (line == "END") break;

    // Split the line into words

    vector<string> tokens;
    ZtoJPsiGamma::tokenize(line, tokens);
    int vsize = tokens.size();
    assert(vsize > 1);

    const string& key   = tokens.at(0);
    const string& value = tokens.at(1);

    if (key == "dataType") {
      string vtmp(value);
      std::transform(vtmp.begin(), vtmp.end(), vtmp.begin(), ::toupper);
      vector<string> dt;
      tokenize(vtmp, dt, "#");
      if (dt.size()) {
        isMC_ = (dt.at(0) == "MC") ? true : false;
        if (isMC_ && dt.size() > 1) {
          isSignal_ = (dt.at(1) == "SIGNAL") ? true : false;
        }
      }
    }

    //else if (key == "readTrigObject")
    //readTrigObject_ = std::stoi(value.c_str()) > 0 ? true : false;
    //  else if (key == "readGenInfo")
    //readGenInfo_ = std::stoi(value.c_str()) > 0 ? true : false;
    //else if (key == "usePUWt")
    //	  usePUWt_ = std::stoi(value.c_str()) > 0 ? true : false;
    /*else if (key == "readPFObject")
      readPFObject_ = std::stoi(value.c_str()) > 0 ? true : false;
      else if (key == "useTrigger")
      useTrigger_ = std::stoi(value.c_str()) > 0 ? true : false;
      else if (key == "useLumiWt")
      useLumiWt_ = std::stoi(value.c_str()) > 0 ? true : false;
       else if (key == "usePUWt")
	  usePUWt_ = std::stoi(value.c_str()) > 0 ? true : false;
        else if (key == "maxEvent")
        maxEvt_ = std::stoi(value.c_str());
	else if (key == "startEvent")
      firstEvt_ = std::stoi(value.c_str());
      else if (key == "endEvent")
      lastEvt_ = std::stoi(value.c_str());
      else if (key == "bunchX")
      bunchCrossing_ = std::stoi(value.c_str());
    else if (key == "trigPathList")
    AnaUtil::buildList(tokens, trigPathList_);*/
    //else if (key == "xs")
    //xsec_ = std::stoi(value.c_str());
    else if (key == "histFile")
      histFile_ = value;
    else if (key == "hlt")
      hlt_ = std::stoi(value.c_str());
    else if (key == "maxEvent")
      maxEvt_ = std::stoi(value.c_str());
    //    else if (key == "puHistFile")
    // puHistFile_ = value;
    //else if (key == "puHistogram")
    // puHistogram_ = value;
    else if (key == "inputFile")
      ZtoJPsiGamma::buildList(tokens, fileList_);
    //else if (key == "eventId" && tokens.size() == 4) {
    //AnaUtil::buildMap(tokens, eventIdMap_);
    //}
    //else if (key == "useEventList")
    //useEventList_ = std::stoi(value.c_str()) > 0 ? true : false;
    //    else if (key == "inputEventFile")
    //eventFileList_.push_back(value);

    else {
      if (0) cout << "==> " << line << endl;
      ZtoJPsiGamma::storeCuts(tokens, hmap_);
    }
  }
  
  
  
  
  
  // Close the file                                                  
  fin.close();

  //if (!isMC_) usePUWt_ = false;
  //if (!isSignal_) readGenInfo_ = false;

  // Build the chain of root files
  for (const auto& fname: fileList_) {
    cout << ">>> INFO. Adding input file " << fname << " to TChain " << endl;
    ++nFiles;
    int nevt = setInputFile(fname);
    if (maxEvt_ > 0 && nevt >= maxEvt_) break;
  }

  if (!nFiles) {
    cerr << ">>> WARN. Input Root file list is empty! exiting ..." << endl;
    return false;
  }

  return true;



}

// ************************************

void ZtoJPsiGamma::Loop()
{
//   In a ROOT session, you can do:
//      root> .L ZtoJPsiGamma.C
//      root> ZtoJPsiGamma t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    chain->SetBranchStatus("*",0);  // disable all branches
//    chain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    chain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


   if (chain == 0) return;
   //Long64_t nentries = chain->GetEntriesFast();
   int nentries = static_cast<int>(chain->GetEntries());

   if (maxEvt_ > 0) nentries = std::min(nentries, maxEvt_);
   //std::cout << "max event: " << maxEvt_ << std::endl;
   std::cout << "Total number of entries: " << nentries << std::endl;
   
   outfile = TFile::Open(histFile_.c_str(), "RECREATE");
   
   Float_t finaljpsim,finalm, finalm1,finalm2,finalm3;
   TTree *T = new TTree("T","test final mass tree");
   T->Branch("finaljpsim",&finaljpsim,"finaljpsim/F");
   T->Branch("finalm",&finalm,"finalm/F");
   TTree *T1 = new TTree("T1","final mass tree cat1");
   T1->Branch("finalm1",&finalm1,"finalm1/F");
   TTree *T2 = new TTree("T2","final mass tree cat2");
   T2->Branch("finalm2",&finalm2,"finalm2/F");
   TTree *T3 = new TTree("T3","final mass tree cat3");
   T3->Branch("finalm3",&finalm3,"finalm3/F");
   //xsec = xsec_;
   
   bookHistograms();

   int num_pu_vtx = 0;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries; jentry++) {
   //      for (Long64_t jentry=0; jentry<100;jentry++) {
     //     if ( !(jentry == 43 || jentry == 44) ) continue;
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = chain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;
     //std::cout << "Analyzing entry: " << jentry << std::endl;
     //std::cout << "debug2 : " << std::endl;
     std::vector<float> vec_genmupt;
     std::vector<float> vec_genphopt;
     std::vector<float> vec_genmueta;
     std::vector<float> vec_genphoeta;
     TLorentzVector genmuP4,  genphoP4;
     std::vector<TLorentzVector> vec_genmuP4;
     std::vector<TLorentzVector> vec_genphoP4;
     //std::cout << "debug3 : " << std::endl;
     
     double wt = 1.0;
     num_pu_vtx = Pileup_nPU;
     //cout << "PU vertices = " << num_pu_vtx << endl;
     //cout << "PU weight = " << dopuweight(num_pu_vtx) << endl;
     //cout << "PU weight from nanoaod branch = " << puWeight << endl;
     //     if (isMC()) wt = wt * Generator_weight * dopuweight(num_pu_vtx);
     if (isMC()) wt = wt * Generator_weight * puWeight;
     else wt = 1.0;
     
     //std::cout << "weight : " << wt << endl
     //	       << "Generator_weight : " << Generator_weight << endl
     //	       << "dopuweight(num_pu_vtx) : " << dopuweight(num_pu_vtx) << endl;

     double evtweight = -99.;

     if (Generator_weight != 1.0) evtweight = Generator_weight;
     else evtweight = 1.0;
     
     evtWeightSum_ += evtweight;
     

     outfile->cd();
     outfile->cd("ZtoJPsiGamma");

     fillHist1D("generatorweight", Generator_weight);
     fillHist1D("puweight", puWeight);
     fillHist2D("pudist", num_pu_vtx, puWeight);
     //fillHist1D("puweight", dopuweight(num_pu_vtx));
     //     fillHist1D("genweight", genWeight);

     fillHist1D("evtCutFlow", 0);
     fillHist1D("evtCutFlowWt", 0, wt);

     // gen part

#if 0          
     std::vector<int> vec_genmupdg;
     for (UInt_t i=0; i< nGenPart ; i++){
       //  std::cout << "Gen particle with id : " << i
       //<< " gen pdg id : " << GenPart_pdgId[i] 
       //<< " gen mass : " << GenPart_mass[i]
       //<< " gen status : " << GenPart_status[i]
       //<< " gen mother id : " << GenPart_genPartIdxMother[i]
       //<< std::endl;
       if (GenPart_pdgId[i] == 22 ){
	 int genphomo1indx = GenPart_genPartIdxMother[i];
	 if (GenPart_pdgId[genphomo1indx] == 25){
	   //genPhoPt->Fill(GenPart_pt[i]);
	   //vec_phopt.push_back(GenPart_pt[i]);
	   //vec_phoeta.push_back(GenPart_eta[i]);
	   genphoP4.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
	   vec_genphoP4.push_back(genphoP4);
	 }
       }
       
       if ((GenPart_pdgId[i] == 13 || GenPart_pdgId[i] == -13))
	 {
	   int genmu1mo1indx = GenPart_genPartIdxMother[i];
	   //std::cout << "genmu1mo1indx : " << genmu1mo1indx << std::endl;
	   if (GenPart_pdgId[genmu1mo1indx] == 443) {
	     int genmu1mo2indx = GenPart_genPartIdxMother[genmu1mo1indx];
	     //std::cout << "genmu1mo2indx : " << genmu1mo2indx << std::endl;
	     if (GenPart_pdgId[genmu1mo2indx] == 25){
	       //vec_mupt.push_back(GenPart_pt[i]);
	       //vec_mueta.push_back(GenPart_eta[i]);
	       genmuP4.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
	       vec_genmupdg.push_back(GenPart_pdgId[i]);
	       vec_genmuP4.push_back(genmuP4);
	       //genMu1Pt->Fill(GenPart_pt[i]);
	       //std::cout << "gen mu pt : " << GenPart_pt[i] << std::endl;
	     }	
	     else if (GenPart_pdgId[genmu1mo2indx] == 443){
	       int genmu1mo3indx = GenPart_genPartIdxMother[genmu1mo2indx];
	       //std::cout << "genmu1mo3indx : " << genmu1mo3indx << std::endl;
	       if (GenPart_pdgId[genmu1mo3indx] == 25){
		 //vec_mupt.push_back(GenPart_pt[i]);
		 //vec_mueta.push_back(GenPart_eta[i]);
		 genmuP4.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
		 vec_genmupdg.push_back(GenPart_pdgId[i]);
		 vec_genmuP4.push_back(genmuP4);
		 //genMu1Pt->Fill(GenPart_pt[i]);                                                                                             
		 //std::cout << "gen mu pt : " << GenPart_pt[i] << std::endl;
	       }
	     }
	   }
	 }
       
     }
     
     if ( vec_genmuP4.size() > 1) std::sort(vec_genmuP4.begin(), vec_genmuP4.end(), PtComparatorTL<TLorentzVector>());
     if ( vec_genphoP4.size() > 1) std::sort(vec_genphoP4.begin(), vec_genphoP4.end(), PtComparatorTL<TLorentzVector>());
     
     //     if (! ((vec_mupt.size() > 1 && vec_mupt[0] >= vec_mupt[1]) || (vec_mupt.size() > 1 && vec_mupt[0] < vec_mupt[1]))) cout << "jentry : " << jentry << endl;
     
     if ( vec_genmuP4.size() > 0) {
       fillHist1D("genMu1Pt", vec_genmuP4[0].Pt());
        fillHist1D("genMu1Eta", vec_genmuP4[0].Eta());
     }
     
     if ( vec_genmuP4.size() > 1) {
       fillHist1D("genMu2Pt", vec_genmuP4[1].Pt());
       fillHist1D("genMu2Eta", vec_genmuP4[1].Eta());
       fillHist1D("gendRMu1Mu2", vec_genmuP4[0].DeltaR(vec_genmuP4[1]));
     }
     if (vec_genphoP4.size() > 0) {
       fillHist1D("genPhoPt", vec_genphoP4[0].Pt());
       fillHist1D("genPhoEta", vec_genphoP4[0].Eta());
     }
     if ( vec_genmuP4.size() > 0 && vec_genphoP4.size() > 0 ){
       fillHist1D("gendRMu1Pho", vec_genphoP4[0].DeltaR(vec_genmuP4[0]));
     }
     if ( vec_genmuP4.size() > 1 && vec_genphoP4.size() > 0 ){
       fillHist1D("gendRdiMuPho", vec_genphoP4[0].DeltaR((vec_genmuP4[0]+vec_genmuP4[1])));
       fillHist1D("gendRMu2Pho", vec_genphoP4[0].DeltaR(vec_genmuP4[1]));
     }
     
     TLorentzVector genz, genjpsi, genmupo;
     float jpsihel;
     if (vec_genmuP4.size() > 1) {
       genjpsi = vec_genmuP4[0] + vec_genmuP4[1];
       if (vec_genmupdg[0] == -13) genmupo = vec_genmuP4[0];
       else if (vec_genmupdg[1] == -13) genmupo = vec_genmuP4[1];
       jpsihel = costheta1(genmupo, genjpsi);
       fillHist1D("costheta1", jpsihel);
       fillHist1D("costheta1_rewt", jpsihel, (3./2)*(1-jpsihel*jpsihel));
       //       fillHist1D("costheta1_rewt", jpsihel, (3./4)*(1+jpsihel*jpsihel)); // for higgs
     }
     float zhel;
     if (vec_genmuP4.size() > 1 && vec_genphoP4.size() > 0) {
       genjpsi = vec_genmuP4[0] + vec_genmuP4[1];
       genz = vec_genmuP4[0] + vec_genmuP4[1] + vec_genphoP4[0];
       zhel = costhetastar(genjpsi, genz);
       fillHist1D("costhetas", zhel);
       fillHist1D("costhetas_rewt", zhel, (3./4)*(1+zhel*zhel));
     }
#endif
     //if (vec_genmuP4.size() < 2 || vec_genphoP4.size() < 1) continue;
     //     fillHist1D("evtCutFlow", 1);
     //fillHist1D("evtCutFlowWt", 1, wt);

     //if (std::abs(vec_genmuP4[0].Eta()) > 2.5 || std::abs(vec_genmuP4[1].Eta()) > 2.5 || std::abs(vec_genphoP4[0].Eta()) > 2.5  ) continue;
     //fillHist1D("evtCutFlow", 2);
     //fillHist1D("evtCutFlowWt", 2, wt);

     fillHist1D("nPV", PV_npvs, wt);
     fillHist1D("nPVgood", PV_npvsGood, wt);
     if (PV_npvsGood < 1) continue;
     fillHist1D("evtCutFlow", 1);
     fillHist1D("evtCutFlowWt", 1, wt);

     bool HLT;
     //     if (hlt() == 0) HLT = HLT_Mu17_Photon30_IsoCaloId;
     if (hlt() == 0) HLT = HLT_Mu17_Photon30_CaloIdL_L1ISO;
     else if (hlt() == 1) HLT = HLT_Dimuon25_Jpsi;
     else if (hlt() == 2) HLT = HLT_DoubleMu20_7_Mass0to30_Photon23;
     else if (hlt() == 3) HLT = HLT_Dimuon20_Jpsi;
     else if (hlt() == 4) HLT = HLT_Mu17_Photon30_IsoCaloId;
     else std::cout << "HLT not found" << std::endl;
     //  if (!(HLT_DoubleMu20_7_Mass0to30_Photon23 || HLT_SingleMu24 || HLT_SingleMu27 || HLT_Dimuon25_Jpsi || HLT_Mu17_Photon30_IsoCaloId || HLT_DoubleMu4_3_Jpsi)) continue;
     if (!HLT) continue;
     //if (!HLT_Mu17_Photon30_IsoCaloId) continue;
     //     if (!HLT_Dimuon25_Jpsi) continue;
     //     if (!HLT_DoubleMu20_7_Mass0to30_Photon23) continue;
     fillHist1D("evtCutFlow", 2);
     fillHist1D("evtCutFlowWt", 2, wt);
     

     //std::vector<float> vec_mupt;
       //std::vector<float> vec_phopt;
     TLorentzVector muLP4, muTP4, muTIsoP4, phoP4;
     std::vector<TLorentzVector> vec_muLP4, vec_muTP4, vec_muTIsoP4;
     std::vector<int> vec_muid;
     std::vector<TLorentzVector> vec_phoP4;
     
     //std::cout << "debug _1 : " << std::endl;
     
     outfile->cd();
     outfile->cd("ObjectSelection");
     float mupt_max = -999.;
     int mu1id = -999;
     for (int i=0; i < nMuon; i++){
       fillHist1D("muCutFlow", 0);
       
       fillHist1D("muPt", Muon_pt[i], wt);
       fillHist1D("muEta", Muon_eta[i], wt);
       fillHist1D("muDxy", Muon_dxy[i], wt);
       fillHist1D("muDz", Muon_dz[i], wt);
       fillHist1D("muSIP3d", Muon_sip3d[i], wt);
       fillHist1D("muSIP3dincm", Muon_ip3d[i], wt);
       //       fillHist1D("muLooseid", Muon_looseId[i], wt);
       if (!Muon_looseId[i]) continue; // Muon POG Loose ID
       fillHist1D("muCutFlow", 1);
       muLP4.SetPtEtaPhiM( Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
       vec_muLP4.push_back(muLP4);

       //std::cout << "Muon_highPtId : " << Muon_highPtId[i] << std::endl; 
       
       //fillHist1D("muIspfcand", Muon_isPFcand[i], wt);
       //fillHist1D("muHighptid", Muon_highPtId[i], wt);
       //fillHist1D("muTightid", Muon_tightId[i], wt);

       //       if (!(Muon_pt[i] < 200 && Muon_isPFcand[i])) continue;
       //if (Muon_pt[i] > 200) continue;
       //       if (Muon_pt[i] > 200 && Muon_highPtId[i] != 1) continue;
       //if (Muon_pt[i] < 200 && !Muon_isPFcand[i]) continue;
       //if (Muon_pt[i] > 200 && !(Muon_highPtId[i] == 1 || Muon_isPFcand[i])) continue;
 
       //if (!(Muon_isPFcand[i] || (Muon_highPtId[i] == 1 && Muon_pt[i] > 200))) continue; // HZZ4L Tight ID
       if (!Muon_tightId[i]) continue;
       fillHist1D("muCutFlow", 2);
       
       muTP4.SetPtEtaPhiM( Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
       vec_muTP4.push_back(muTP4);

       if (Muon_pt[i] > mupt_max){
	 mupt_max = Muon_pt[i];
	 mu1id = i;
       }
       
     }

     int mu2id = -999.;
     int mu2ndmaxpt = -999.;
     if (nMuon > 0) {
       for (int j=0; j < nMuon; j++){
	 if (j==mu1id) continue;
	 if (Muon_pt[j] > mu2ndmaxpt){
	   mu2ndmaxpt = Muon_pt[j];
	   mu2id = j;
	 }
       }
     }
     outfile->cd();
     outfile->cd("ObjectSelection");     
     
     float phopt_max = -999.;
     int pho1id = -999;     

     for (int i = 0 ; i < nPhoton ; i++){
       fillHist1D("photonCutFlow", 0);

       //fillHist1D("phoMvaidwp90", Photon_mvaID_WP90[i], wt);     
       if (Photon_mvaID_WP90[i] == false) continue;
       fillHist1D("photonCutFlow", 1);
       
       fillHist1D("phoPt", Photon_pt[i], wt);
       if (Photon_pt[i] <= 33) continue;
       fillHist1D("photonCutFlow", 2);
       
       if (Photon_isScEtaEB[i]) fillHist1D("nPhoScEta", 0, wt); // 0 = ScEtaEB, 1 = ScEtaEE
       else if (Photon_isScEtaEE[i]) fillHist1D("nPhoScEta", 1, wt);

       //       fillHist1D("Photon_isScEtaEB", Photon_isScEtaEB[i], wt);
       //fillHist1D("Photon_isScEtaEE", Photon_isScEtaEE[i], wt);
       if (! (Photon_isScEtaEB[i] || Photon_isScEtaEE[i])) continue;
       fillHist1D("photonCutFlow", 3);

       fillHist1D("phoR9", Photon_r9[i], wt);
       if (Photon_isScEtaEB[i]) fillHist1D("phoR9_ScEB", Photon_r9[i], wt);
       if (Photon_isScEtaEE[i]) fillHist1D("phoR9_ScEE", Photon_r9[i], wt);

       phoP4.SetPtEtaPhiM(Photon_pt[i], Photon_eta[i], Photon_phi[i], Photon_mass[i]);
       vec_phoP4.push_back(phoP4);
       
       if (Photon_pt[i] > phopt_max){
	 phopt_max = Photon_pt[i];
	 pho1id = i;
       }


     }
     //     std::cout << "# photons  : " << vec_phoP4.size() << std::endl;



     outfile->cd();
     outfile->cd("ZtoJPsiGamma");
     
     //if (vec_muTP4.size()>0) std::cout << "lead mu pt from P4 : " << vec_muTP4[0].Pt() << "lead mu pt from index : " << Muon_pt[mu1id] << std::endl;
     fillHist1D("nLoosemuons", vec_muLP4.size(), wt);
     fillHist1D("nTightmuons", vec_muTP4.size(), wt);
     if ( vec_muTP4.size() < 2 ) continue;
     fillHist1D("evtCutFlow", 3);
     fillHist1D("evtCutFlowWt", 3, wt);
    
     std::sort(vec_muTP4.begin(), vec_muTP4.end(), PtComparatorTL<TLorentzVector>());
    
     // std::cout << "lead mu pt from P4 : " << vec_muTP4[0].Pt() << "lead mu pt from index : " << Muon_pt[mu1id] << std::endl;
     
     //std::cout << "(Muon_pfRelIso03_all[mu1id])  : " << Muon_pfRelIso03_all[mu1id] << std::endl;

     //     std::cout << "mu1 charge " << Muon_charge[mu1id] << "  mu2 charge " << Muon_charge[mu2id] << std::endl;
     //     if (Muon_charge[mu1id] + Muon_charge[mu2id] != 0) continue;
     //std::cout << "mu1 charge " << Muon_charge[mu1id] << "  mu2 charge " << Muon_charge[mu2id] << std::endl;
     //std::cout << "end of this event" << std::endl;

     fillHist1D("muPfRelIso03all", Muon_pfRelIso03_all[mu1id], wt);
     fillHist1D("muPfRelIso03chg", Muon_pfRelIso03_chg[mu1id], wt);
     fillHist1D("muPfRelIso04all", Muon_pfRelIso04_all[mu1id], wt);

     fillHist1D("muPfRelIso03all", Muon_pfRelIso03_all[mu1id], wt);
     fillHist1D("muPfRelIso03allByMuSubleadpt", Muon_pfRelIso03_all[mu1id] / vec_muTP4[1].Pt(), wt);
     if (Muon_pfRelIso03_all[mu1id] >= 0.35) continue;
     fillHist1D("evtCutFlow", 4);
     fillHist1D("evtCutFlowWt", 4, wt);
     
     fillHist1D("muLeadpt", vec_muTP4[0].Pt(), wt);
     //     fillHist1D("muLeadpt", vec_muTP4[0].Pt(), wt); //test
 
     if (vec_muTP4[0].Pt() <= 20) continue;
     fillHist1D("evtCutFlow", 5);
     fillHist1D("evtCutFlowWt", 5, wt);

     fillHist1D("muSubleadpt", vec_muTP4[1].Pt(), wt);     
     if (vec_muTP4[1].Pt() <= 4) continue;
     fillHist1D("evtCutFlow", 6);
     fillHist1D("evtCutFlowWt", 6, wt);


    fillHist1D("nPhotons", vec_phoP4.size(), wt);
     if ( vec_phoP4.size() < 1) continue;
     fillHist1D("evtCutFlow", 7);
     fillHist1D("evtCutFlowWt", 7, wt);
     std::sort(vec_phoP4.begin(), vec_phoP4.end(), PtComparatorTL<TLorentzVector>());

     
     outfile->cd();
     outfile->cd("ZtoJPsiGamma");

     fillHist1D("dRmu1pho", vec_muTP4[0].DeltaR(vec_phoP4[0]), wt);
     if (vec_muTP4[0].DeltaR(vec_phoP4[0]) <= 1 ) continue;
     fillHist1D("evtCutFlow", 8);
     fillHist1D("evtCutFlowWt", 8, wt);

 
     fillHist1D("dRmu2pho", vec_muTP4[1].DeltaR(vec_phoP4[0]), wt);
     if (vec_muTP4[1].DeltaR(vec_phoP4[0]) <= 1 ) continue;
     fillHist1D("evtCutFlow", 9);     
     fillHist1D("evtCutFlowWt", 9, wt);
     
     fillHist1D("dRdimupho", (vec_muTP4[0] + vec_muTP4[1]).DeltaR(vec_phoP4[0]), wt);
      if ((vec_muTP4[0] + vec_muTP4[1]).DeltaR(vec_phoP4[0]) <= 2 ) continue;
     fillHist1D("evtCutFlow", 10);  
     fillHist1D("evtCutFlowWt", 10, wt);
     
     fillHist1D("dPhidimupho", (vec_muTP4[0] + vec_muTP4[1]).DeltaPhi(vec_phoP4[0]), wt);
     if (std::fabs((vec_muTP4[0] + vec_muTP4[1]).DeltaPhi(vec_phoP4[0])) <= 1.5 ) continue;
     fillHist1D("evtCutFlow", 11);
     fillHist1D("evtCutFlowWt", 11, wt);

     //fillHist1D("jpsimass_bmasscut", (vec_muTP4[0] + vec_muTP4[1]).M());
     fillHist1D("jpsimass_bmasscut", (vec_muTP4[0] + vec_muTP4[1]).M(), wt);
     
     if ( (vec_muTP4[0] + vec_muTP4[1]).M() < 3.0 || (vec_muTP4[0] + vec_muTP4[1]).M() > 3.2) continue;
 
     fillHist1D("evtCutFlow", 12);
     fillHist1D("evtCutFlowWt", 12, wt);

     fillHist1D("jpsimass_amasscut", (vec_muTP4[0] + vec_muTP4[1]).M(), wt);     

     fillHist1D("zmass_bmasscut", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), wt );

     if ((vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M() <= ZtoJPsiGamma::cutValue(SRCutMap(), "low") || (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M() >= ZtoJPsiGamma::cutValue(SRCutMap(), "high")) continue;
     fillHist1D("evtCutFlow", 13);
     fillHist1D("evtCutFlowWt", 13, wt);

     fillHist1D("zmass_amasscut", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), wt );


     fillHist1D("ratiocutdimu", (vec_muTP4[0] + vec_muTP4[1]).Pt()/(vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), wt);
     if ((vec_muTP4[0] + vec_muTP4[1]).Pt()/(vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M() <= ZtoJPsiGamma::cutValue(ratioCutMap(), "ptratio")) continue;
     fillHist1D("evtCutFlow", 14);
     fillHist1D("evtCutFlowWt", 14, wt);

     fillHist1D("ratiocutpho", vec_phoP4[0].Pt()/(vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), wt);
     if (vec_phoP4[0].Pt()/(vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M() <= ZtoJPsiGamma::cutValue(ratioCutMap(), "etratio")) continue;
     fillHist1D("evtCutFlow", 15);          
     fillHist1D("evtCutFlowWt", 15, wt);

     finaljpsim = (vec_muTP4[0] + vec_muTP4[1]).M();
     finalm = (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M();
     //

     if (vec_phoP4[0].Pt()/(vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M() > ZtoJPsiGamma::cutValue(ratioCutMap(), "etratio")) { fillHist1D("nPhoCat", 0); fillHist1D("nPhoCatWt", 0, wt);}
     if (Photon_isScEtaEB[pho1id] && Photon_r9[pho1id] >= 0.94 )   { 
       fillHist1D("nPhoCat", 1); 
       fillHist1D("nPhoCatWt", 1, wt);
       fillHist1D("muLeadpt_cat1", vec_muTP4[0].Pt(), wt);
       fillHist1D("muLeadeta_cat1", vec_muTP4[0].Eta(), wt);
       fillHist1D("muSubleadpt_cat1", vec_muTP4[1].Pt(), wt);
       fillHist1D("muSubleadeta_cat1", vec_muTP4[1].Eta(), wt);
       fillHist1D("phoPt_cat1", vec_phoP4[0].Pt(), wt);
       fillHist1D("phoEta_cat1", vec_phoP4[0].Eta(), wt);
       fillHist1D("dRmu1pho_cat1", vec_muTP4[0].DeltaR(vec_phoP4[0]), wt);
       fillHist1D("dRmu2pho_cat1", vec_muTP4[1].DeltaR(vec_phoP4[0]), wt);
       fillHist1D("dRdimupho_cat1", (vec_muTP4[0] + vec_muTP4[1]).DeltaR(vec_phoP4[0]), wt);
       fillHist1D("dPhidimupho_cat1", (vec_muTP4[0] + vec_muTP4[1]).DeltaPhi(vec_phoP4[0]), wt);
       fillHist1D("jpsimass_cat1", (vec_muTP4[0] + vec_muTP4[1]).M(), wt);
       fillHist1D("dRmu1mu2_cat1", vec_muTP4[0].DeltaR(vec_muTP4[1]), wt); 
       fillHist1D("diMuPt_cat1", (vec_muTP4[0] + vec_muTP4[1]).Pt(), wt);
       fillHist1D("finalZPt_cat1", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).Pt(), wt);
       fillHist1D("zmass_cat1", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), wt);
       fillHist1D("vtxdist_x" , PV_x - (vec_muTP4[0] + vec_muTP4[1]).X(), wt);
       fillHist1D("vtxdist_y" , PV_y - (vec_muTP4[0] + vec_muTP4[1]).Y(), wt);
       fillHist1D("vtxdist_z" , PV_z - (vec_muTP4[0] + vec_muTP4[1]).Z(), wt);
       finalm1 = (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M();
       T1->Fill();
     }
     else if (Photon_isScEtaEB[pho1id] && Photon_r9[pho1id] < 0.94 ){
       fillHist1D("nPhoCat", 2); 
       fillHist1D("nPhoCatWt", 2, wt); 
       fillHist1D("muLeadpt_cat2", vec_muTP4[0].Pt(), wt);
       fillHist1D("muLeadeta_cat2", vec_muTP4[0].Eta(), wt);
       fillHist1D("muSubleadpt_cat2", vec_muTP4[1].Pt(), wt);
       fillHist1D("muSubleadeta_cat2", vec_muTP4[1].Eta(), wt);
       fillHist1D("phoPt_cat2", vec_phoP4[0].Pt(), wt);
       fillHist1D("phoEta_cat2", vec_phoP4[0].Eta(), wt);
       fillHist1D("dRmu1pho_cat2", vec_muTP4[0].DeltaR(vec_phoP4[0]), wt);
       fillHist1D("dRmu2pho_cat2", vec_muTP4[1].DeltaR(vec_phoP4[0]), wt);
       fillHist1D("dRdimupho_cat2", (vec_muTP4[0] + vec_muTP4[1]).DeltaR(vec_phoP4[0]), wt);
       fillHist1D("dPhidimupho_cat2", (vec_muTP4[0] + vec_muTP4[1]).DeltaPhi(vec_phoP4[0]), wt);
       fillHist1D("jpsimass_cat2", (vec_muTP4[0] + vec_muTP4[1]).M(), wt);
       fillHist1D("dRmu1mu2_cat2", vec_muTP4[0].DeltaR(vec_muTP4[1]), wt);
       fillHist1D("diMuPt_cat2", (vec_muTP4[0] + vec_muTP4[1]).Pt(), wt);
       fillHist1D("finalZPt_cat2", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).Pt(), wt);
       fillHist1D("zmass_cat2", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), wt);
       finalm2 = (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M();
       T2->Fill();
     } 
     else if (Photon_isScEtaEE[pho1id] ) {
       fillHist1D("nPhoCat", 3); 
       fillHist1D("nPhoCatWt", 3, wt);
       fillHist1D("muLeadpt_cat3", vec_muTP4[0].Pt(), wt);
       fillHist1D("muLeadeta_cat3", vec_muTP4[0].Eta(), wt);
       fillHist1D("muSubleadpt_cat3", vec_muTP4[1].Pt(), wt);
       fillHist1D("muSubleadeta_cat3", vec_muTP4[1].Eta(), wt);
       fillHist1D("phoPt_cat3", vec_phoP4[0].Pt(), wt);
       fillHist1D("phoEta_cat3", vec_phoP4[0].Eta(), wt);
       fillHist1D("dRmu1pho_cat3", vec_muTP4[0].DeltaR(vec_phoP4[0]), wt);
       fillHist1D("dRmu2pho_cat3", vec_muTP4[1].DeltaR(vec_phoP4[0]), wt);
       fillHist1D("dRdimupho_cat3", (vec_muTP4[0] + vec_muTP4[1]).DeltaR(vec_phoP4[0]), wt);
       fillHist1D("dPhidimupho_cat3", (vec_muTP4[0] + vec_muTP4[1]).DeltaPhi(vec_phoP4[0]), wt);
       fillHist1D("jpsimass_cat3", (vec_muTP4[0] + vec_muTP4[1]).M(), wt);
       fillHist1D("dRmu1mu2_cat3", vec_muTP4[0].DeltaR(vec_muTP4[1]), wt);
       fillHist1D("diMuPt_cat3", (vec_muTP4[0] + vec_muTP4[1]).Pt(), wt);
       fillHist1D("finalZPt_cat3", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).Pt(), wt);
       fillHist1D("zmass_cat3", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), wt);
       finalm3 = (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M();
       T3->Fill();
     }
     //final plots
     
     T->Fill();
     
     fillHist1D("jpsimass_afinalcut", (vec_muTP4[0] + vec_muTP4[1]).M(), wt);
     
     //     fillHist1D("zmass_afinalcut", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M());     
     fillHist1D("zmass_afinalcut", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), wt);
     // fillHist1D("zmass_afinalcut", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), wt*((intLumi * xsec )/ 200000));
     
     //To check jpsi is coming from prompt vertex
     
     
     
     /*     fillHist1D("patdRmu1mu2_bHLT", vec_muTP4[0].DeltaR(vec_muTP4[1]));
	    fillHist1D("patMmu1mu2_bHLT", (vec_muTP4[0] + vec_muTP4[1]).M());
	    fillHist1D("patPTmu1mu2_bHLT", (vec_muTP4[0] + vec_muTP4[1]).Pt());
	    
	    //    if (!HLT_Mu17_Photon30_IsoCaloId) continue;
	    //   if (!HLT_Dimuon25_Jpsi) continue;
	    if (!HLT_DoubleMu20_7_Mass0to30_Photon23) continue;
	    
	    fillHist1D("patMmu1mu2_aHLT", (vec_muTP4[0] + vec_muTP4[1]).M());
	    fillHist1D("patPTmu1mu2_aHLT", (vec_muTP4[0] + vec_muTP4[1]).Pt());    
	    
	    if ( (vec_muTP4[0]+vec_muTP4[1]).M() < 3.0 || (vec_muTP4[0]+vec_muTP4[1]).M() > 3.2) continue;
	    //    if ( (vec_muTP4[0]+vec_muTP4[1]).M() > 3.2) continue;
	    
	    fillHist1D("patdRmu1mu2_aHLT", vec_muTP4[0].DeltaR(vec_muTP4[1]));*/
     
     
     
     
   }
   saveHistograms();
   outfile->cd();
   outfile->cd("ObjectSelection");
   objectEfficiency();
   showCategoryYield();
   outfile->cd();
   outfile->cd("ZtoJPsiGamma");
   T->Write();
   eventEfficiency();
   outfile->Close();
   cout << "Now the ouput file is closed" << endl;
   
   
}
