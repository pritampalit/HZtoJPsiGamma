#include <iostream>
#include <string>
#include <memory>
#include <utility>
#include "TROOT.h"
#include "TStopwatch.h"
#include "HtoJPsiGamma.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;

int main(int argc, char* argv[]) {
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " jobFile " << endl;
    exit(1);
    }
  string jobFile(argv[1]);

  //  TChain* chain = new TChain("Events");
  // Create  analysis object                                                                                                                
  HtoJPsiGamma anaH;

  // Read job input                                                 
  int nFiles = 0;
  bool succeed = anaH.readJob(jobFile, nFiles);
  if (!succeed) exit(2);
  //  if (anaH.getEntries() <= 0) {
  //cerr << "No events present in the input chain, exiting ...!" << endl;
  //exit(3);
  //}


  gROOT->SetBatch(kTRUE);

  // Now go                                                                                                                                  
  TStopwatch timer;
  cout << "Start event loop now with " << nFiles << " file(s)" << endl;
  timer.Start();

  // Initialize analysis                                                                                                                   
  //if (!anaH.beginJob()) exit(4);
  
  //for (int i=0; i< fileList_.length(); i++){
  //chain->Add(i);
  //}
  
  // anaH(chain);

  anaH.Loop();

  // Analysis over                                                                                                                           
  //anaH.endJob();
  //anaH.closeFiles();

  timer.Stop();
  cout << "Realtime/CpuTime = " << timer.RealTime() << "/" << timer.CpuTime() << endl;

  exit(0);
}

