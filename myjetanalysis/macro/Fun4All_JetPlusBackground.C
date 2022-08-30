#ifndef MACRO_FUN4ALLJETANA_C
#define MACRO_FUN4ALLJETANA_C

// CreateFileList.pl -n 1000 -type 6  DST_VERTEX DST_CALO_G4HIT DST_CALO_CLUSTER DST_BBC_G4HIT

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

/* #include <GlobalVariables.C> */
/* #include <G4_Global.C> */
#include <g4jets/FastJetAlgo.h>
#include <g4jets/JetReco.h>
#include <g4jets/TowerJetInput.h>
#include <g4jets/TruthJetInput.h>

#include <g4centrality/PHG4CentralityReco.h>

#include <jetbackground/FastJetAlgoSub.h>

// here you need your package name (set in configure.ac)
#include <jetplusbackground/JetPlusBackground.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libjetplusbackground.so)

void Fun4All_JetPlusBackground(
    const int nevnt = 10, 
    const double min_calo_pt=0.02, 
    const int verbosity=1,
    const char* fout_name="out/JetPlusBackground.root")
{
  gSystem->Load("libjetplusbackground");
  gSystem->Load("libg4dst");

  Fun4AllServer *se = Fun4AllServer::instance();

  PHG4CentralityReco *cent = new PHG4CentralityReco();
  cent->Verbosity(0);
  cent->GetCalibrationParameters().ReadFromFile("centrality", "xml", 0, 0, string(getenv("CALIBRATIONROOT")) + string("/Centrality/"));
  se->registerSubsystem( cent );

  // change lower pt and eta cut to make them visible using the example
  //  pythia8 file
  int print_stats_freq = 20;
  JetPlusBackground *myJetAnalysis = new JetPlusBackground(min_calo_pt, nevnt, print_stats_freq, "AntiKt_Tower_r04", fout_name);
  myJetAnalysis->setPtRange(5, 100);
  myJetAnalysis->setEtaRange(-1.1, 1.1);
  myJetAnalysis->add_input(new TowerJetInput(Jet::CEMC_TOWER));
  myJetAnalysis->add_input(new TowerJetInput(Jet::HCALIN_TOWER));
  myJetAnalysis->add_input(new TowerJetInput(Jet::HCALOUT_TOWER));
  se->registerSubsystem(myJetAnalysis);

  // need truth jets
  // need calo  jets
  // need event info
  // need primary vertex
// $ CreateFileList.pl -run 4 -n 1000 -type 11 -embed DST_VERTEX DST_CALO_G4HIT DST_CALO_CLUSTER DST_TRUTH_JET DST_BBC_G4HIT
//    PHG4CentralityReco::InitRun : cannot find G4HIT_BBC, will not use MBD centrality
//    PHG4CentralityReco::InitRun : cannot find G4HIT_EPD, will not use sEPD centrality
    
  /* Fun4AllInputManager *intrue = new Fun4AllDstInputManager("DSTtruth"); */
  /* intrue->AddListFile("dst_truth_jet.list",1); // adding the option "1" confirms to use, even if file is large */
  /* se->registerInputManager(intrue); */

  Fun4AllInputManager *incalo = new Fun4AllDstInputManager("DSTcalo");
  incalo->AddListFile("lists/dst_calo_g4hit.list",1);
  se->registerInputManager(incalo);

  Fun4AllInputManager *incalocluster = new Fun4AllDstInputManager("DSTcalocluster");
  incalocluster->AddListFile("lists/dst_calo_cluster.list",1);
  se->registerInputManager(incalocluster);

  Fun4AllInputManager *invertex = new Fun4AllDstInputManager("DSTvertex");
  invertex->AddListFile("lists/dst_vertex.list",1);
  se->registerInputManager(invertex);
    
  Fun4AllInputManager *inbbc = new Fun4AllDstInputManager("DSTbbc");
  inbbc->AddListFile("lists/dst_bbc_g4hit.list",1);
  se->registerInputManager(inbbc);

  myJetAnalysis->Verbosity(verbosity);
  se->run(nevnt);
  se->End();
  delete se;
  gSystem->Exit(0);
}

#endif
