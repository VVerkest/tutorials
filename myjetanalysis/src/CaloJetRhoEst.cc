#include "CaloJetRhoEst.h"

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

#include <g4jets/JetMap.h>
#include <g4jets/JetInput.h>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace std;
using namespace fastjet;

CaloJetRhoEst::CaloJetRhoEst(const std::string& recojetname, const std::string& truthjetname, const std::string& outputfilename)
  : SubsysReco("CaloJetRhoEst_" + recojetname + "_" + truthjetname)
  , m_recoJetName(recojetname)
  , m_truthJetName(truthjetname)
  , m_outputFileName(outputfilename)
  , m_etaRange (-1, 1)
  , m_ptRange  (5,  100)
  , m_T  (nullptr)
  , m_id (-1)
  , m_eta       {}
  , m_phi       {}
  , m_e         {}
  , m_pt        {}
  , m_area      {}
  , m_truthEta  {}
  , m_truthPhi  {}
  , m_truthE    {}
  , m_truthPt   {}
  , m_truthArea {}
  , _inputs {}
{ }

CaloJetRhoEst::~CaloJetRhoEst()
{ 
  for (unsigned int i = 0; i < _inputs.size(); ++i) delete _inputs[i];
  _inputs.clear();
}

int CaloJetRhoEst::Init(PHCompositeNode* topNode)
{
  if (Verbosity() >= CaloJetRhoEst::VERBOSITY_SOME)
    cout << "CaloJetRhoEst::Init - Outoput to " << m_outputFileName << endl;

  PHTFileServer::get().open(m_outputFileName, "RECREATE");

  //Tree
  m_T = new TTree("T", "CaloJetRhoEst Tree");

  //      int m_event;
  m_T->Branch("m_id",          &m_id);
  m_T->Branch("m_rho",         &m_rho);
  m_T->Branch("m_rho_sigma",   &m_rho_sigma);
  m_T->Branch("m_centrality",  &m_centrality);
  m_T->Branch("m_impactparam", &m_impactparam);

  m_T->Branch("m_rawEta",    &m_eta);
  m_T->Branch("m_rawPhi",    &m_phi);
  m_T->Branch("m_rawE",      &m_e);
  m_T->Branch("m_rawPt",     &m_pt);
  m_T->Branch("m_rawArea",   &m_area);

  //Truth Jets
  m_T->Branch("m_truthEta",  &m_truthEta);
  m_T->Branch("m_truthPhi",  &m_truthPhi);
  m_T->Branch("m_truthE",    &m_truthE);
  m_T->Branch("m_truthPt",   &m_truthPt);
  m_T->Branch("m_truthArea", &m_truthArea);

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloJetRhoEst::End(PHCompositeNode* topNode)
{
  cout << "CaloJetRhoEst::End - Output to " << m_outputFileName << endl;
  PHTFileServer::get().cd(m_outputFileName);

  /* m_hInclusiveE->Write(); */
  /* m_hInclusiveEta->Write(); */
  /* m_hInclusivePhi->Write(); */
  m_T->Write();

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloJetRhoEst::InitRun(PHCompositeNode* topNode)
{
  /* m_jetEvalStack = shared_ptr<JetEvalStack>(new JetEvalStack(topNode, m_recoJetName, m_truthJetName)); */
  /* m_jetEvalStack->get_stvx_eval_stack()->set_use_initial_vertex(initial_vertex); */
  topNode->print();
  cout << " Input Selections:" << endl;
  for (unsigned int i = 0; i < _inputs.size(); ++i) _inputs[i]->identify();
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloJetRhoEst::process_event(PHCompositeNode* topNode)
{
  /* return Fun4AllReturnCodes::EVENT_OK; // FIXME :: just printing the nodes for now in order to find them */
  /* ++m_id; */
//  JetMap* jets = findNode::getClass<JetMap>(topNode, m_recoJetName);
//  if (!jets)
//  {
//    std::cout
//      << "MyJetAnalysis::process_event - Error can not find DST Reco JetMap node "
//      << m_recoJetName << std::endl;
//    exit(-1);
//  }

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

  // now make pseudojet particles 
  //      (as in from /direct/sphenix+u/dstewart/vv/coresoftware/simulation/g4simulation/g4jets/JetReco.cc ::94 ->
  //                  /direct/sphenix+u/dstewart/vv/coresoftware/offline/packages/jetbackground/FastJetAlgoSub.h :: get_jets
  auto& particles=inputs;
  /* cout << " particles: " << particles.size() << endl; */

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

      /* if (_verbosity > 5) */
      /* { */
      /*   std::cout << " FastJetAlgoSub input particle with negative-E, original kinematics px / py / pz / E = "; */
      /*   std::cout << particles[ipart]->get_px() << " / " << particles[ipart]->get_py() << " / " << particles[ipart]->get_pz() << " / " << particles[ipart]->get_e() << std::endl; */
      /*   std::cout << " --> entering with modified kinematics px / py / pz / E = " << this_px << " / " << this_py << " / " << this_pz << " / " << this_e << std::endl; */
      /* } */
    }

    fastjet::PseudoJet pseudojet(this_px, this_py, this_pz, this_e);

    pseudojet.set_user_index(ipart);

    float _pt = pseudojet.perp();
    if (_pt < 0.002) {
      /* cout << " CUT SMALL: " << _pt << endl; */
      ++smallptcutcnt;
    } else {
      particles_pseudojets.push_back(pseudojet);
    }
  }

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

  //get reco jets
  // cout << " olives A0 " << endl;
//  for (JetMap::Iter iter = jets->begin(); iter != jets->end(); ++iter)
//  {
//    Jet* jet = iter->second;
//    float pt = jet->get_pt();
//    float eta = jet->get_eta();
//    if  (pt < m_ptRange.first  || pt  > m_ptRange.second
//        || eta < m_etaRange.first || eta > m_etaRange.second) continue;
//    m_pt.push_back(pt);
//    m_eta.push_back(eta);
//    m_phi.push_back(jet->get_phi());
//    m_e.push_back(jet->get_e());
//  }
    
  // for now, it appears that the pT is stored in ascending pT order
  // to avoid the guarantee, just check for lead and sublead as goes along
  // Make a vector of Jet*, sorted in descending pt
  vector<Jet*> truth_jets;
  for (JetMap::Iter iter = jetsMC->begin(); iter != jetsMC->end(); ++iter)
  {
    Jet* truthjet = iter->second;
    float pt = truthjet->get_pt();
    float eta = truthjet->get_eta();
    if  (pt < m_ptRange.first  
        || pt  > m_ptRange.second 
        || eta < m_etaRange.first 
        || eta > m_etaRange.second) continue;
    truth_jets.push_back(truthjet);
  }
  std::sort(truth_jets.begin(), truth_jets.end(), [](Jet* a, Jet* b) { return a->get_pt() > b->get_pt(); });


  Jet* leadJet    = (truth_jets.size()>0 ? truth_jets[0] : nullptr);
  Jet* subLeadJet = (truth_jets.size()>1 ? truth_jets[1] : nullptr);
  for (auto jet : truth_jets) {
    m_truthPt .push_back(jet->get_pt());
    m_truthEta.push_back(jet->get_eta());
    /* cout << " olives: A2 MC " << truthjet->get_eta() << endl; */
    m_truthPhi.push_back(jet->get_phi());
    m_truthE  .push_back(jet->get_e());
  }

  if (false) cout << leadJet->get_pt() << " " << subLeadJet->get_pt() << endl;
    
  JetDefinition jet_def(cambridge_algorithm, 0.4);     //  JET DEFINITION
<<<<<<< HEAD
    
  Selector leadCircle = SelectorCircle(0.4);
  if(have_lead) { leadCircle.set_reference(leadJet); }
  Selector subCircle = SelectorCircle(0.4);
  if(have_sub) { subCircle.set_reference(subLeadJet); }
  Selector bgRapRange = SelectorRapRange( -0.6, 0.6 );
  Selector bgSelector = bgRapRange && !leadCircle && !subCircle;
  double ghost_maxrap = 1.0;
  AreaDefinition area_def(active_area, GhostedAreaSpec(ghost_maxrap));
//  JetMedianBackgroundEstimator UE( bgSelector, jet_def, area_def);
//  UE.set_jets(pseudojets);
//  cout<<UE.rho()<<endl;
    
  clear_vectors();
=======
>>>>>>> 91920cd47ce641f43bc91865642ef89a177063d2

  /* Selector leadCircle = SelectorCircle(0.4); */
  /* if(have_lead) { leadCircle.set_reference(leadJet); } */
  /* Selector subCircle = SelectorCircle(0.4); */
  /* if(have_sub) { subCircle.set_reference(subLeadJet); } */
  /* Selector bgRapRange = SelectorRapRange( -0.6, 0.6 ); */
  /* Selector bgSelector = bgRapRange && !leadCircle && !subCircle; */
  /* double ghost_maxrap = 4.0; */

  /* Selector bgRapRange = SelectorRapRange( -0.6, 0.6 ); */
  /* Selector bgSelector = bgRapRange && !leadCircle && !subCircle; */
  const double ghost_max_rap { 4.0 };
  const double ghost_R = 0.01;
  const double jet_R = 0.4;
  AreaDefinition area_def_bkgd( active_area_explicit_ghosts, GhostedAreaSpec(ghost_max_rap, 1, ghost_R));
  JetDefinition jet_def_bkgd(cambridge_algorithm, jet_R); // <--
  Selector selector_rm2 = SelectorAbsEtaMax(0.6) * (!SelectorNHardest(2)); // <--
  fastjet::JetMedianBackgroundEstimator bge_rm2 {selector_rm2, jet_def_bkgd, area_def_bkgd};
  bge_rm2.set_particles(particles_pseudojets);

  m_rho = bge_rm2.rho();
  m_rho_sigma = bge_rm2.sigma();
  /* cout << " got: " << rho << " " << rho_sigma << endl; */


  // cluster the measured jets:
  double max_rap = 1.6;
  fastjet::Selector jetrap         = fastjet::SelectorAbsEtaMax(0.6);
  fastjet::Selector not_pure_ghost = !SelectorIsPureGhost();
  fastjet::Selector selection      = jetrap && not_pure_ghost;
  AreaDefinition area_def( active_area_explicit_ghosts, GhostedAreaSpec(max_rap, 1, ghost_R));
  JetDefinition jet_def_antikt(antikt_algorithm, jet_R);
  fastjet::ClusterSequenceArea clustSeq(particles_pseudojets, jet_def_antikt, area_def);
  vector<PseudoJet> jets = sorted_by_pt( selection( clustSeq.inclusive_jets(m_ptRange.first) ));
  for (auto jet : jets) {
    m_eta  .push_back( jet.eta());
    m_phi  .push_back( jet.phi_std());
    m_e    .push_back( jet.E());
    m_pt   .push_back( jet.pt());
    m_area .push_back( jet.area());
  }

  m_T->Fill();
  clear_vectors();
  return Fun4AllReturnCodes::EVENT_OK;
}


void CaloJetRhoEst::clear_vectors() {
  m_eta.clear();
  m_phi.clear();
  m_e.clear();
  m_pt.clear();
  m_area.clear();

  
  m_truthEta.clear();
  m_truthPhi.clear();
  m_truthE.clear();
  m_truthPt.clear();
  m_truthArea.clear();
}
