// Veronica Verkest
// Augusgt 31, 2022

#include <vector>
#include <iostream>
#include <TH1D.h>

const int nbins_centrality = 10;
const double bins_centrality[nbins_centrality+1] = { 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0 };
const string name_centrality[nbins_centrality] = { "_0_10", "_10_20", "_20_30", "_30_40", "_40_50", "_50_60", "_60_70", "_70_80", "_80_90", "_90_100" };
const string title_centrality[nbins_centrality] = { "0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100%" };
const int color_centrality[nbins_centrality] = { 51, 56, 61, 66, 71, 76, 81, 86, 91, 96};

const int nbins_rhopt = 80;
const double bins_rhopt[nbins_rhopt+1] = { -40.0, -39.0, -38.0, -37.0, -36.0, -35.0, -34.0, -33.0, -32.0, -31.0, -30.0, -29.0, -28.0, -27.0, -26.0, -25.0, -24.0, -23.0, -22.0, -21.0, -20.0, -19.0, -18.0, -17.0, -16.0, -15.0, -14.0, -13.0, -12.0, -11.0, -10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0 };
const int nbins_rhopt_40 = 40;
const double bins_dpt_40[nbins_rhopt_40+1] = { -40.0, -38.0, -36.0, -34.0, -32.0, -30.0, -28.0, -26.0, -24.0, -22.0, -20.0, -18.0, -16.0, -14.0, -12.0, -10.0, -8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0 };

const double R=0.4;

int get_centrality_bin( float centrality ) {
    if ( (centrality<bins_centrality[0]) || (centrality>bins_centrality[nbins_centrality]) ) {
        std::cout<<"ERROR: CENTRALITY MUST BE GIVEN AS A FLOAT OR DOUBLE BETWEEN 0 AND 100"<<std::endl;
        return -99;
    }
    for (int i=0; i<nbins_centrality; ++i) {
        if ( (centrality>bins_centrality[i]) && (centrality<bins_centrality[i+1]) ) { return i; }
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

void ExploreTree_JetPlusBackground() {
    
    TH1::SetDefaultSumw2(); TH2::SetDefaultSumw2(); TH3::SetDefaultSumw2();
    
    TString name, title; // free temporary variables

    TFile *inFile = new TFile("out/jet_bg/JetPlusBackground_136.root","READ");
    TTree *T = (TTree*)inFile->Get("T");

    int id;
    float rho, rho_sigma, centrality, impactparam, rhoBias_lead, rhoBias_sub;
    vector<float> *CaloJetEta = NULL;
    vector<float> *CaloJetPhi = NULL;
    vector<float> *CaloJetE = NULL;
    vector<float> *CaloJetPt = NULL;
    vector<float> *CaloJetArea = NULL;
    float embEta_A, embPhi_A, embPt_A, embEta_B, embPhi_B, embPt_B;
    
    T->SetBranchAddress("id",&id);
    T->SetBranchAddress("rho",&rho);
    T->SetBranchAddress("rho_sigma",&rho_sigma);
    T->SetBranchAddress("centrality",&centrality);
    T->SetBranchAddress("impactparam",&impactparam);
    T->SetBranchAddress("rhoBias_lead",&rhoBias_lead);
    T->SetBranchAddress("rhoBias_sub",&rhoBias_sub);
    T->SetBranchAddress("CaloJetEta",&CaloJetEta);
    T->SetBranchAddress("CaloJetPhi",&CaloJetPhi);
    T->SetBranchAddress("CaloJetE",&CaloJetE);
    T->SetBranchAddress("CaloJetPt",&CaloJetPt);
    T->SetBranchAddress("CaloJetArea",&CaloJetArea);
    T->SetBranchAddress("embEta_A",&embEta_A);
    T->SetBranchAddress("embPhi_A",&embPhi_A);
    T->SetBranchAddress("embPt_A",&embPt_A);
    T->SetBranchAddress("embEta_B",&embEta_B);
    T->SetBranchAddress("embPhi_B",&embPhi_B);
    T->SetBranchAddress("embPt_B",&embPt_B);
        
    TH1D *h_rhoBias_lead[nbins_centrality];
    TH1D *h_ptBias_lead[nbins_centrality];
    TH1D *h_rhoBias_sub[nbins_centrality];
    TH1D *h_ptBias_sub[nbins_centrality];
    TH2D *h_where_are_jets = new TH2D("h_where_are_jets","",3,0.,3.,3,0.,3.);
    
    h_where_are_jets->GetXaxis()->SetBinLabel(1,"A in lead");
    h_where_are_jets->GetXaxis()->SetBinLabel(2,"A in sub");
    h_where_are_jets->GetXaxis()->SetBinLabel(3,"A not found");
    h_where_are_jets->GetYaxis()->SetBinLabel(1,"B in lead");
    h_where_are_jets->GetYaxis()->SetBinLabel(2,"B in sub");
    h_where_are_jets->GetYaxis()->SetBinLabel(3,"B not found");
    
    for (int i=0; i<nbins_centrality; ++i) {
        name = "h_rhoBias_lead" + name_centrality[i];
        h_rhoBias_lead[i] = new TH1D(name,";(1/A_{calo,lead}) [p_{T}^{probe} - (p_{T}^{calo,lead} - #rho * A_{calo,lead})]",nbins_rhopt,bins_rhopt);
        name = "h_ptBias_lead" + name_centrality[i];
        h_ptBias_lead[i] = new TH1D(name,";p_{T}^{probe} - (p_{T}^{calo,lead} - #rho * A_{calo,lead})",nbins_rhopt,bins_rhopt);
        name = "h_rhoBias_sub" + name_centrality[i];
        h_rhoBias_sub[i] = new TH1D(name,";(1/A_{calo,sub}) [p_{T}^{probe} - (p_{T}^{calo,sub} - #rho * A_{calo,sub})]",nbins_rhopt,bins_rhopt);
        name = "h_ptBias_sub" + name_centrality[i];
        h_ptBias_sub[i] = new TH1D(name,";p_{T}^{probe} - (p_{T}^{calo,sub} - #rho * A_{calo,sub})",nbins_rhopt,bins_rhopt);
        
        h_rhoBias_lead[i]->SetLineColor(color_centrality[i]);
        h_rhoBias_lead[i]->SetMarkerColor(color_centrality[i]);
        h_rhoBias_lead[i]->SetMarkerStyle(20);
        
        h_ptBias_lead[i]->SetLineColor(color_centrality[i]);
        h_ptBias_lead[i]->SetMarkerColor(color_centrality[i]);
        h_ptBias_lead[i]->SetMarkerStyle(20);
        
        h_rhoBias_sub[i]->SetLineColor(color_centrality[i]);
        h_rhoBias_sub[i]->SetMarkerColor(color_centrality[i]);
        h_rhoBias_sub[i]->SetMarkerStyle(21);
        
        h_ptBias_sub[i]->SetLineColor(color_centrality[i]);
        h_ptBias_sub[i]->SetMarkerColor(color_centrality[i]);
        h_ptBias_sub[i]->SetMarkerStyle(21);
//        [i]->Write();
//        h_rhoBias_sub[i]->Write();
//        h_ptBias_sub[i]->Write();
    }

    for (int ientry=0; ientry<T->GetEntries(); ++ientry) {
        
        T->GetEntry(ientry);
        
        bool have_sublead = CaloJetPt->size()-1; // check for a second hard jet
        int cent = get_centrality_bin( centrality );

        double leadPt = CaloJetPt->at(0);
        double leadEta = CaloJetEta->at(0);
        double leadPhi = CaloJetPhi->at(0); //        double leadE = CaloJetE->at(0);
        double leadArea = CaloJetArea->at(0);
        
        double subPt = 0.0;
        double subEta = 0.0;
        double subPhi = 0.0; //        double subE = NULL;
        double subArea = 0.0;
        if (have_sublead) {  // we need >1 jet to fill these
            subPt = CaloJetPt->at(1);
            subEta = CaloJetEta->at(1);
            subPhi = CaloJetPhi->at(1); //        double subE = CaloJetE->at(1);
            subArea = CaloJetArea->at(1);
        }
        
        bool A_in_lead = match_jet( embEta_A, embPhi_A, leadEta, leadPhi, R );
        bool B_in_lead = match_jet( embEta_B, embPhi_B, leadEta, leadPhi, R );

        bool A_in_sub = false;
        bool B_in_sub = false;
        if (have_sublead) {
            A_in_sub = match_jet( embEta_A, embPhi_A, subEta, subPhi, R );
            B_in_sub = match_jet( embEta_B, embPhi_B, subEta, subPhi, R );
        }

        double fill_A_bin, fill_B_bin;
        if (A_in_lead) { fill_A_bin = 0.5; }
        else if (A_in_sub) { fill_A_bin = 1.5; }
        else { fill_A_bin = 2.5; }
        
        if (B_in_lead) { fill_B_bin = 0.5; }
        else if (B_in_sub) { fill_B_bin = 1.5; }
        else { fill_B_bin = 2.5; }
        
        h_where_are_jets->Fill(fill_A_bin,fill_B_bin); // cout<<A_in_lead<<" \t "<<A_in_sub<<" \t "<<B_in_lead<<" \t "<<B_in_sub<<endl;
        
        if ( (A_in_lead && B_in_sub) || (B_in_lead && A_in_sub) ) { // if A and B are the 2 highest pT jets, fill
            h_rhoBias_lead[cent]->Fill(rhoBias_lead);
            h_rhoBias_sub[cent]->Fill(rhoBias_sub);
            if ( A_in_lead && B_in_sub ) {
                h_ptBias_lead[cent]->Fill( embPt_A - leadPt );
                h_ptBias_sub[cent]->Fill( embPt_B - subPt );
            }
            else if ( B_in_lead && A_in_sub ) {
                h_ptBias_lead[cent]->Fill( embPt_B - leadPt );
                h_ptBias_sub[cent]->Fill( embPt_A - subPt );
            }
        }


    }
    
    TFile *outFile = new TFile("analyze/exploreTree.root","RECREATE");

    h_where_are_jets->Write();
    for (int i=0; i<nbins_centrality; ++i) {
        h_rhoBias_lead[i]->Write();
        h_ptBias_lead[i]->Write();
        h_rhoBias_sub[i]->Write();
        h_ptBias_sub[i]->Write();
    }

    h_rhoBias_lead[0]->Scale(1./h_rhoBias_lead[0]->Integral());
    h_rhoBias_lead[0]->GetYaxis()->SetRangeUser( 0., 1.);
    h_rhoBias_lead[0]->SetAxisRange( 0., 1., "Y");
    h_rhoBias_lead[0]->Scale(1./h_rhoBias_lead[0]->Integral());
    h_rhoBias_lead[0]->GetYaxis()->SetRangeUser( 0., 1.);
    h_rhoBias_lead[0]->SetAxisRange( 0., 1., "Y");
    for (int i=0; i<nbins_centrality; ++i) {
        h_rhoBias_lead[i]->Scale(1./h_rhoBias_lead[i]->Integral());
        h_rhoBias_lead[i]->SetAxisRange( 0., 1., "Y");
        h_rhoBias_lead[i]->Draw("PSAME");

    }
    new TCanvas();
    for (int i=0; i<nbins_centrality; ++i) {
        h_rhoBias_sub[i]->Scale(1./h_rhoBias_sub[i]->Integral());
        h_rhoBias_sub[i]->SetAxisRange( 0., 1., "Y");
        h_rhoBias_sub[i]->Draw("PSAME");
//        h_ptBias_lead[i]->;
//        h_ptBias_sub[i]->;
    }
}
