#include "TracksAndCaloJets.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>

// fastjet includes
#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/BackgroundEstimatorBase.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <g4eval/JetEvalStack.h>

#include <trackbase_historic/SvtxTrackMap.h>

#include <centrality/CentralityInfo.h>

//#include <g4jets/JetMap.h>
#include <g4jets/JetMapv1.h>
#include <g4jets/JetInput.h>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>

#include <TRandom3.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace std;
using namespace fastjet;

TracksAndCaloJets::TracksAndCaloJets(
      const double _min_calo_pt,
      const int total_jobs,
      const int n_print_freq,
      const std::string& recojetname,
      const std::string& truthjetname,
      const std::string& outputfilename
    )
  : SubsysReco("TracksAndCaloJets_" + recojetname + "_" + truthjetname)
  , min_calo_pt {_min_calo_pt}
  , m_recoJetName(recojetname)
  , m_truthJetName(truthjetname)
  , m_outputFileName(outputfilename)
  , m_etaRange (-1, 1)
  , m_ptRange  (5,  100)
  , m_T  (nullptr)
  , m_id (-1)
  , m_CaloJetEta   {}
  , m_CaloJetPhi   {}
  , m_CaloJetE     {}
  , m_CaloJetPt    {}
  , m_CaloJetArea  {}
  , m_TruthJetEta  {}
  , m_TruthJetPhi  {}
  , m_TruthJetE    {}
  , m_TruthJetPt   {}
  , m_TruthJetArea {}
  , _inputs        {}
  , print_stats{n_print_freq, total_jobs}
{ 
}

TracksAndCaloJets::~TracksAndCaloJets()
{ 
  for (unsigned int i = 0; i < _inputs.size(); ++i) delete _inputs[i];
  _inputs.clear();
    print_stats.set_get_stats();
   cout << " Max memory used: " << print_stats.max_mem/1000. << " MB " << endl;

}

int TracksAndCaloJets::Init(PHCompositeNode* topNode)
{
  if (Verbosity() >= TracksAndCaloJets::VERBOSITY_SOME)
    cout << "TracksAndCaloJets::Init - Outoput to " << m_outputFileName << endl;

  PHTFileServer::get().open(m_outputFileName, "RECREATE");

  //Tree
  m_T = new TTree("T", "TracksAndCaloJets Tree");

  //      int m_event;
  m_T->Branch("id",          &m_id);
  m_T->Branch("rho",         &m_rho);
  m_T->Branch("rho_sigma",   &m_rho_sigma);
  m_T->Branch("centrality",  &m_centrality);
  m_T->Branch("impactparam", &m_impactparam);

  m_T->Branch("CaloJetEta",    &m_CaloJetEta);
  m_T->Branch("CaloJetPhi",    &m_CaloJetPhi);
  m_T->Branch("CaloJetE",      &m_CaloJetE);
  m_T->Branch("CaloJetPt",     &m_CaloJetPt);
  m_T->Branch("CaloJetArea",   &m_CaloJetArea);

  //Truth Jets
  m_T->Branch("TruthJetEta",  &m_TruthJetEta);
  m_T->Branch("TruthJetPhi",  &m_TruthJetPhi);
  m_T->Branch("TruthJetE",    &m_TruthJetE);
  m_T->Branch("TruthJetPt",   &m_TruthJetPt);
  m_T->Branch("TruthJetArea", &m_TruthJetArea);

  return Fun4AllReturnCodes::EVENT_OK;
}

int TracksAndCaloJets::End(PHCompositeNode* topNode)
{
  cout << "TracksAndCaloJets::End - Output to " << m_outputFileName << endl;
  PHTFileServer::get().cd(m_outputFileName);

  /* m_hInclusiveE->Write(); */
  /* m_hInclusiveEta->Write(); */
  /* m_hInclusivePhi->Write(); */
  m_T->Write();

  return Fun4AllReturnCodes::EVENT_OK;
}

int TracksAndCaloJets::InitRun(PHCompositeNode* topNode)
{
  topNode->print();
  cout << " Input Selections:" << endl;
  for (unsigned int i = 0; i < _inputs.size(); ++i) _inputs[i]->identify();
  return Fun4AllReturnCodes::EVENT_OK;
}

int TracksAndCaloJets::process_event(PHCompositeNode* topNode)
{
  /* cout << "Verbosity: " << Verbosity() << endl; */
  // statistics on how the program is doing
  print_stats.call();
  //interface to truth jets
  JetMap* jetsMC = findNode::getClass<JetMap>(topNode, m_truthJetName);
  if (!jetsMC )
  {
    std::cout
      << "MyJetAnalysis::process_event - Error can not find DST Truth JetMap node "
      << m_truthJetName << std::endl;
    exit(-1);
  }

  // get the inputs for reconstructed jets (from /direct/sphenix+u/dstewart/vv/coresoftware/simulation/g4simulation/g4jets/JetReco.cc
  std::vector<Jet *> inputs;  // owns memory
  for (unsigned int iselect = 0; iselect < _inputs.size(); ++iselect)
  {
    std::vector<Jet *> parts = _inputs[iselect]->get_input(topNode);
    for (unsigned int ipart = 0; ipart < parts.size(); ++ipart)
    {
      inputs.push_back(parts[ipart]);
      inputs.back()->set_id(inputs.size() - 1);  // unique ids ensured
    }
  }

  auto& particles=inputs;

  // /direct/sphenix+u/dstewart/vv/coresoftware/offline/packages/jetbackground/FastJetAlgoSub.cc ::58
  std::vector<fastjet::PseudoJet> particles_pseudojets;
  int   smallptcutcnt =0;//FIXME
  for (unsigned int ipart = 0; ipart < particles.size(); ++ipart)
  {
    float this_e = particles[ipart]->get_e();

    if (this_e == 0.) continue;

    float this_px = particles[ipart]->get_px();
    float this_py = particles[ipart]->get_py();
    float this_pz = particles[ipart]->get_pz();

    if (this_e < 0)
    {
      // make energy = +1 MeV for purposes of clustering
      float e_ratio = 0.001 / this_e;

      this_e  = this_e * e_ratio;
      this_px = this_px * e_ratio;
      this_py = this_py * e_ratio;
      this_pz = this_pz * e_ratio;
    }

    fastjet::PseudoJet pseudojet(this_px, this_py, this_pz, this_e);

    pseudojet.set_user_index(ipart);

    float _pt = pseudojet.perp();
    if (_pt < min_calo_pt) {
      /* cout << " CUT SMALL: " << _pt << " < " << min_calo_pt << endl; */
      ++smallptcutcnt;
    } else {
      particles_pseudojets.push_back(pseudojet);
    }
  }
  for (auto &p : particles) delete p;

  //centrality
   CentralityInfo* cent_node = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
   if (!cent_node)
   {
     std::cout
       << "MyJetAnalysis::process_event - Error can not find centrality node "
       << std::endl;
     exit(-1);
   }

  //get the event centrality/impact parameter from HIJING
   m_centrality  =  cent_node->get_centile(CentralityInfo::PROP::bimp);
   m_impactparam =  cent_node->get_quantity(CentralityInfo::PROP::bimp);

  vector<Jet*> truth_jets;
  for (auto& jet : jetsMC->vec(Jet::SORT::PT)) {  // will return jets in order of descending pT
    float pt  = jet->get_pt();
    float eta = jet->get_eta();
    if  (pt < m_ptRange.first
        || pt  > m_ptRange.second
        || eta < m_etaRange.first
        || eta > m_etaRange.second) continue;
    truth_jets.push_back(jet);
  }

  Jet* leadJet    = (truth_jets.size()>0 ? truth_jets[0] : nullptr);
  Jet* subLeadJet = (truth_jets.size()>1 ? truth_jets[1] : nullptr);
  if (false) cout << leadJet->get_pt() << " " << subLeadJet->get_pt() << endl;

  for (auto jet : truth_jets) {
    m_TruthJetPt .push_back(jet->get_pt());
    m_TruthJetEta.push_back(jet->get_eta());
    /* cout << " olives: A2 MC " << truthjet->get_eta() << endl; */
    m_TruthJetPhi.push_back(jet->get_phi());
    m_TruthJetE  .push_back(jet->get_e());
  }

  if (Verbosity()>5) cout << "Starting background density calc" << endl;
  JetDefinition jet_def(kt_algorithm, 0.4);     //  JET DEFINITION

  const double ghost_max_rap { 2.0 };
  const double ghost_R = 0.01;
  const double jet_R = 0.4;
  AreaDefinition area_def_bkgd( active_area_explicit_ghosts, GhostedAreaSpec(ghost_max_rap, 1, ghost_R));
  JetDefinition jet_def_bkgd(kt_algorithm, jet_R); // <--
  Selector selector_rm2 = SelectorAbsEtaMax(0.6) * (!SelectorNHardest(2)); // <--
  fastjet::JetMedianBackgroundEstimator bge_rm2 {selector_rm2, jet_def_bkgd, area_def_bkgd};
  bge_rm2.set_particles(particles_pseudojets);

  m_rho = bge_rm2.rho();
  m_rho_sigma = bge_rm2.sigma();

  if (Verbosity()>5) cout << "Starting clustered jets" << endl;
  // cluster the measured jets:
  double max_rap = 2.0;
  fastjet::Selector jetrap         = fastjet::SelectorAbsEtaMax(0.6);
  fastjet::Selector not_pure_ghost = !SelectorIsPureGhost();
  fastjet::Selector selection      = jetrap && not_pure_ghost;
  AreaDefinition area_def( active_area_explicit_ghosts, GhostedAreaSpec(max_rap, 1, ghost_R));
  JetDefinition jet_def_antikt(antikt_algorithm, jet_R);
  fastjet::ClusterSequenceArea clustSeq(particles_pseudojets, jet_def_antikt, area_def);
  vector<PseudoJet> jets = sorted_by_pt( selection( clustSeq.inclusive_jets(m_ptRange.first) ));
  for (auto jet : jets) {
    m_CaloJetEta  .push_back( jet.eta());
    m_CaloJetPhi  .push_back( jet.phi_std());
    m_CaloJetE    .push_back( jet.E());
    m_CaloJetPt   .push_back( jet.pt());
    m_CaloJetArea .push_back( jet.area());
  }

  m_T->Fill();
  clear_vectors();

  return Fun4AllReturnCodes::EVENT_OK;
}

void TracksAndCaloJets::clear_vectors() {
  m_CaloJetEta.clear();
  m_CaloJetPhi.clear();
  m_CaloJetE.clear();
  m_CaloJetPt.clear();
  m_CaloJetArea.clear();
  
  m_TruthJetEta.clear();
  m_TruthJetPhi.clear();
  m_TruthJetE.clear();
  m_TruthJetPt.clear();
  m_TruthJetArea.clear();
}
