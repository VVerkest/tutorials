// Veronica Verkest
// September 21, 2022

// Takes in PYTHIA and truth jets (without embedding) and smears the PYTHIA jets
// by sampling from the HIJING deltaPt distribution. The delta pT distribution
// comes from the HIJING BG event with a manually-embedded high-pT jet. Using BG
// estimation, delta pT = pT,calo - rho*A - pT,truth, where pT,calo is of a jet
// geometrically matched to the embedded jet. All is done in centrality bins.

#include <vector>
#include <iostream>
#include <TH1D.h>
#include <TH2D.h>

// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ CONSTANTS ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
const int nbins_centrality = 10;
const double bins_centrality[nbins_centrality+1] = { 0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100. };
const string name_centrality[nbins_centrality] = { "_0_10", "_10_20", "_20_30", "_30_40", "_40_50", "_50_60", "_60_70", "_70_80", "_80_90", "_90_100" };
const TString title_centrality[nbins_centrality] = { "0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100%" };
const int color_centrality[nbins_centrality] = { 51, 56, 61, 66, 71, 76, 81, 86, 91, 96};

const int nbins_rhopt = 80;
const double bins_rhopt[nbins_rhopt+1] = { -40., -39., -38., -37., -36., -35., -34., -33., -32., -31., -30., -29., -28., -27., -26., -25., -24., -23., -22., -21., -20., -19., -18., -17., -16., -15., -14., -13., -12., -11., -10., -9., -8., -7., -6., -5., -4., -3., -2., -1., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40. };

const double R=0.4;

// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ FUNCTIONS ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
//TString gausEq( double mean, double sigma ) {
//
//    TString mu = Form("%f",mean);
//    TString sig = Form("%f",sigma);
//
//    TString exp = "-0.5*pow( (x - " + mu + ")/" + sig + ",2)";
////    TString exp = "-0.5*pow( (x - "; exp += mu; exp += ")/"; exp += sig; exp += ",2)";
//    TString denom = "(" + sig + "*sqrt(2.*" + Form("%f",M_PI) + "))";
//    TString eq = "exp(" + exp + ")/" + denom;
//    return eq;
//};

void FillSmearedResponse(double jetPt, double truthPt, TH2D *hSmeared, TH1D *hDeltaPt) {
    double mean_dPt = 0.;
    for (int i=1; i<=hDeltaPt->GetNbinsX(); ++i) {
        if (hDeltaPt->GetBinContent(i)==0) { continue; }
        double delta_pt = hDeltaPt->GetBinCenter(i);
        double weight = hDeltaPt->GetBinContent(i);
        mean_dPt += (jetPt+delta_pt)*weight;
    }
    hSmeared->Fill(mean_dPt,truthPt);
};

int get_centrality_bin( float centrality ) {
    if ( (centrality<bins_centrality[0]) || (centrality>bins_centrality[nbins_centrality]) ) {
        std::cout<<"ERROR: CENTRALITY MUST BE GIVEN AS A FLOAT OR DOUBLE BETWEEN 0 AND 100"<<std::endl;
        return -99;
    }
    for (int i=0; i<nbins_centrality; ++i) {
        if ( (centrality>=bins_centrality[i]) && (centrality<=bins_centrality[i+1]) ) { return i; }
    }
    return -99; // you shouldn't get here, but necessary
};

double delta_eta( double eta_1, double eta_2 ) { return fabs(eta_1-eta_2); };

double delta_phi( double phi_1, double phi_2 ) {
    double dphi = fabs( phi_1 - phi_2 );
    while (dphi>2.*M_PI) { dphi-= 2.*M_PI; }
    return dphi;
};

double delta_R ( double eta_1, double phi_1, double eta_2, double phi_2 ) {
    return sqrt( delta_eta(eta_1,eta_2)*delta_eta(eta_1,eta_2) + delta_phi(phi_1,phi_2)*delta_phi(phi_1,phi_2) );
};

bool match_jet ( double eta_1, double phi_1, double eta_2, double phi_2, double R=0.4 ) {
    return delta_R( eta_1, phi_1, eta_2, phi_2 )<R;
};


// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ MACRO ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
void SmearJetResponse(){
 
    TH1::SetDefaultSumw2(); TH2::SetDefaultSumw2(); TH3::SetDefaultSumw2();

    TString name, title; // free temporary variables

    //             HIJING (JET PLUS BACKGROUND)
    TFile *hFile = new TFile("out/jet_bg/JetPlusBackground_136.root","READ"); // "h" corresponds to HIJING (or jetPlusBackground)
    TTree *t_jetBG = (TTree*)hFile->Get("T");

    int h_id;
    float h_rho, h_rho_sigma, h_centrality, h_impactparam, rhoBias_lead, rhoBias_sub;
    vector<float> *h_CaloJetEta = NULL;
    vector<float> *h_CaloJetPhi = NULL;
    vector<float> *h_CaloJetE = NULL;
    vector<float> *h_CaloJetPt = NULL;
    vector<float> *h_CaloJetArea = NULL;
    float embEta_A, embPhi_A, embPt_A, embEta_B, embPhi_B, embPt_B;
    
    t_jetBG->SetBranchAddress("id",&h_id);
    t_jetBG->SetBranchAddress("rho",&h_rho);
    t_jetBG->SetBranchAddress("rho_sigma",&h_rho_sigma);
    t_jetBG->SetBranchAddress("centrality",&h_centrality);
    t_jetBG->SetBranchAddress("impactparam",&h_impactparam);
    t_jetBG->SetBranchAddress("rhoBias_lead",&rhoBias_lead);
    t_jetBG->SetBranchAddress("rhoBias_sub",&rhoBias_sub);
    t_jetBG->SetBranchAddress("CaloJetEta",&h_CaloJetEta);
    t_jetBG->SetBranchAddress("CaloJetPhi",&h_CaloJetPhi);
    t_jetBG->SetBranchAddress("CaloJetE",&h_CaloJetE);
    t_jetBG->SetBranchAddress("CaloJetPt",&h_CaloJetPt);
    t_jetBG->SetBranchAddress("CaloJetArea",&h_CaloJetArea);
    t_jetBG->SetBranchAddress("embEta_A",&embEta_A);
    t_jetBG->SetBranchAddress("embPhi_A",&embPhi_A);
    t_jetBG->SetBranchAddress("embPt_A",&embPt_A);
    t_jetBG->SetBranchAddress("embEta_B",&embEta_B);
    t_jetBG->SetBranchAddress("embPhi_B",&embPhi_B);
    t_jetBG->SetBranchAddress("embPt_B",&embPt_B);
    
    TH1D *h_delta_pt[nbins_centrality];
    
    for (int i=0; i<nbins_centrality; ++i) {
        name = "h_delta_pt" + name_centrality[i];
        h_delta_pt[i] = new TH1D(name,";p_{T}^{calo,sub} - #rho * A_{calo,sub} - p_{T}^{truth}",101,-50.5,50.5);
        h_delta_pt[i]->SetLineColor(color_centrality[i]);
        h_delta_pt[i]->SetMarkerColor(color_centrality[i]);
        h_delta_pt[i]->SetMarkerStyle(4);
        h_delta_pt[i]->SetMarkerSize(0.5);
    }

    
    for (int ientry=0; ientry<t_jetBG->GetEntries(); ++ientry) {
        
        t_jetBG->GetEntry(ientry);
        
//        vector<int> break_entry = { 373514, 1306154, 2386357, 4276830, 5905514, 6937999, 7225771, 7646415, 12170773, 15010003, 15763070 };
//        if ( find(break_entry.begin(), break_entry.end(), ientry) != break_entry.end() ) { continue; }

        bool have_sublead = h_CaloJetPt->size()-1; // check for a second hard jet
        int cent = get_centrality_bin( h_centrality );

        double leadPt = h_CaloJetPt->at(0);
        double leadEta = h_CaloJetEta->at(0);
        double leadPhi = h_CaloJetPhi->at(0); //        double leadE = h_CaloJetE->at(0);
        double leadArea = h_CaloJetArea->at(0);
        
        double subPt = 0.;
        double subEta = 0.;
        double subPhi = 0.; //        double subE = NULL;
        double subArea = 0.;
        if (have_sublead) {  // we need >1 jet to fill these
            subPt = h_CaloJetPt->at(1);
            subEta = h_CaloJetEta->at(1);
            subPhi = h_CaloJetPhi->at(1); //        double subE = CaloJetE->at(1);
            subArea = h_CaloJetArea->at(1);
        }
        
        bool A_in_lead = match_jet( embEta_A, embPhi_A, leadEta, leadPhi, R );
        bool B_in_lead = match_jet( embEta_B, embPhi_B, leadEta, leadPhi, R );

        bool A_in_sub = false;
        bool B_in_sub = false;
        if (have_sublead) {
            A_in_sub = match_jet( embEta_A, embPhi_A, subEta, subPhi, R );
            B_in_sub = match_jet( embEta_B, embPhi_B, subEta, subPhi, R );
        }

        if ( (A_in_lead && B_in_sub) || (B_in_lead && A_in_sub) ) { // if A and B are the 2 highest pT jets, fill
            if ( A_in_lead && B_in_sub ) {
                h_delta_pt[cent]->Fill( (leadPt-h_rho*leadArea) - embPt_A );
                h_delta_pt[cent]->Fill( (subPt-h_rho*subArea) - embPt_B );
            }
            else if ( B_in_lead && A_in_sub ) {
                h_delta_pt[cent]->Fill( (leadPt-h_rho*leadArea) - embPt_B );
                h_delta_pt[cent]->Fill( (subPt-h_rho*subArea) - embPt_A );
            }
        }
        
        
    }

    for (int i=0; i<nbins_centrality; ++i) {
        h_delta_pt[i]->Scale(1./h_delta_pt[i]->Integral());
        h_delta_pt[i]->SetAxisRange(0.00001,2.,"Y");
        h_delta_pt[i]->Draw("SAME");
    }
    

    //             PYTHIA JET NO EMBED (NO EMBED)
    TFile *pFile = new TFile("out/CaloJetRho_noEmbed_Sept20.root","READ"); // "p" corresponds to PYTHIA (or noEmbed)
    TTree *t_noEmb = (TTree*)pFile->Get("T");

    int p_id;
    float p_rho, p_rho_sigma, p_centrality, p_impactparam;
    vector<float> *p_CaloJetEta = NULL;
    vector<float> *p_CaloJetPhi = NULL;
    vector<float> *p_CaloJetE = NULL;
    vector<float> *p_CaloJetPt = NULL;
    vector<float> *p_CaloJetArea = NULL;
    vector<float> *TruthJetEta = NULL;
    vector<float> *TruthJetPhi = NULL;
    vector<float> *TruthJetE = NULL;
    vector<float> *TruthJetPt = NULL;
    vector<float> *TruthJetArea = NULL;

    t_noEmb->SetBranchAddress("id",&p_id);
    t_noEmb->SetBranchAddress("rho",&p_rho);
    t_noEmb->SetBranchAddress("rho_sigma",&p_rho_sigma);
    t_noEmb->SetBranchAddress("centrality",&p_centrality);
    t_noEmb->SetBranchAddress("impactparam",&p_impactparam);
    t_noEmb->SetBranchAddress("CaloJetEta",&p_CaloJetEta);
    t_noEmb->SetBranchAddress("CaloJetPhi",&p_CaloJetPhi);
    t_noEmb->SetBranchAddress("CaloJetE",&p_CaloJetE);
    t_noEmb->SetBranchAddress("CaloJetPt",&p_CaloJetPt);
    t_noEmb->SetBranchAddress("CaloJetArea",&p_CaloJetArea);
    t_noEmb->SetBranchAddress("TruthJetEta",&TruthJetEta);
    t_noEmb->SetBranchAddress("TruthJetPhi",&TruthJetPhi);
    t_noEmb->SetBranchAddress("TruthJetE",&TruthJetE);
    t_noEmb->SetBranchAddress("TruthJetPt",&TruthJetPt);
    t_noEmb->SetBranchAddress("TruthJetArea",&TruthJetArea);

    TH2D *hResp_noEmbed = new TH2D("hResp_noEmbed",";p_{T}^{truth};p_{T}^{calo}",100,0.,100.,100,0.,100.);
    TH2D *hResp_noEmbed_smear[nbins_centrality];
    
    for (int i=0; i<nbins_centrality; ++i) {
        name = "hResp_noEmbed_smear" + name_centrality[i];
        hResp_noEmbed_smear[i] = new TH2D(name,";p_{T}^{truth};HIJING BG-smeared p_{T}^{calo}",100,0.,100.,100,0.,100.);
    }

    for (int ientry=0; ientry<t_noEmb->GetEntries(); ++ientry) {
        
        t_noEmb->GetEntry(ientry);

        if (TruthJetPt->size()<=0 || p_CaloJetPt->size()<=0 || TruthJetEta->size()<=0 || p_CaloJetEta->size()<=0 || TruthJetPhi->size()<=0 || p_CaloJetPhi->size()<=0 || TruthJetPt->at(0)<=0 ) { continue; }

        double truthPt = TruthJetPt->at(0);
        double truthEta = TruthJetEta->at(0);
        double truthPhi = TruthJetPhi->at(0); //        double truthE = TruthJetE->at(0);
                
        bool match_truth = false;
        
        for (int i=0; i<p_CaloJetPt->size(); ++i) {

            if (match_truth) { continue; }

            double jetPt = p_CaloJetPt->at(0);
            double jetEta = p_CaloJetEta->at(0);
            double jetPhi = p_CaloJetPhi->at(0); //        double jetE = p_CaloJetE->at(0);
            double jetArea = p_CaloJetArea->at(0);

            if (match_jet(truthEta,truthPhi,jetEta,jetPhi,0.4)) {
                hResp_noEmbed->Fill(jetPt,truthPt);
                for (int i=0; i<nbins_centrality; ++i) { FillSmearedResponse(jetPt,truthPt,hResp_noEmbed_smear[i],h_delta_pt[i]); }
                match_truth = true;
            }
        }

        
    } // end PYTHIA loop

    TH1D *yproj = (TH1D*) hResp_noEmbed->ProjectionY("yproj",1,hResp_noEmbed->GetNbinsY(),"E");
    TH1D *xproj = (TH1D*) hResp_noEmbed->ProjectionX("xproj",1,hResp_noEmbed->GetNbinsX(),"E");
    for (int ix=1; ix<=hResp_noEmbed->GetNbinsX(); ++ix) {
        for (int iy=1; iy<=hResp_noEmbed->GetNbinsY(); ++iy) {
            if (yproj->GetBinContent(iy)==0) { continue; }
            hResp_noEmbed->SetBinContent(ix, iy, hResp_noEmbed->GetBinContent(ix,iy)/yproj->GetBinContent(iy) );
        }
    }
    
    TCanvas *c[nbins_centrality];

    for (int i=0; i<nbins_centrality; ++i) {
        
        TH1D *yproj = (TH1D*) hResp_noEmbed_smear[i]->ProjectionY("yproj",1,hResp_noEmbed_smear[i]->GetNbinsY(),"E");
        for (int ix=1; ix<=hResp_noEmbed_smear[i]->GetNbinsX(); ++ix) {
            for (int iy=1; iy<=hResp_noEmbed_smear[i]->GetNbinsY(); ++iy) {
                if (yproj->GetBinContent(iy)==0) { continue; }
                hResp_noEmbed_smear[i]->SetBinContent(ix, iy, hResp_noEmbed_smear[i]->GetBinContent(ix,iy)/yproj->GetBinContent(iy) );
            }
        }
        
        name = "c" + name_centrality[i];
        c[i] = new TCanvas(name);
        c[i]->SetLogz();
        hResp_noEmbed_smear[i]->Draw("COLZ");
    }
    
    
    TFile *outFile = new TFile("analyze/SmearJetResponse.root","RECREATE");

    hResp_noEmbed->Write();
    for (int i=0; i<nbins_centrality; ++i) {
        hResp_noEmbed_smear[i]->Write();
    }
    
}

