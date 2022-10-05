// Veronica Verkest
// September 21, 2022

// Takes in PYTHIA and truth jets (without embedding) and smears the PYTHIA jets
// by sampling from the HIJING deltaPt distribution. The delta pT distribution
// comes from the HIJING BG event with a manually-embedded high-pT jet. Using BG
// estimation, delta pT = pT,calo - rho*A - pT,truth, where pT,calo is of a jet
// geometrically matched to the embedded jet. All is done in centrality bins.

#include <vector>
#include <iostream>
#include <cstring>
#include <string>
#include <TH1D.h>
#include <TH2D.h>

// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ CONSTANTS ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
const int nbins_centrality = 10;
const double bins_centrality[nbins_centrality+1] = { 0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100. };
const string name_centrality[nbins_centrality] = { "_0_10", "_10_20", "_20_30", "_30_40", "_40_50", "_50_60", "_60_70", "_70_80", "_80_90", "_90_100" };
const string title_centrality[nbins_centrality] = { "0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100%" };
const int color_centrality[nbins_centrality] = { 51, 54, 59, 62, 65, 70, 83, 91, 95, 99};

const int nbins_rhopt = 80;
const double bins_rhopt[nbins_rhopt+1] = { -40., -39., -38., -37., -36., -35., -34., -33., -32., -31., -30., -29., -28., -27., -26., -25., -24., -23., -22., -21., -20., -19., -18., -17., -16., -15., -14., -13., -12., -11., -10., -9., -8., -7., -6., -5., -4., -3., -2., -1., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40. };

const double R=0.4;

const int nbins_pt = 2;
const double bins_pt[nbins_pt+1] = { 30.,40.,60.};
const string name_pt[nbins_pt] = {"_30_40GeV","_40_60GeV"};
const string title_pt[nbins_pt] = {"30-40 GeV p_{T}^{reco}","40-60 GeV p_{T}^{reco}"};
const int marker_pt[nbins_pt] = { 20, 24 };

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

void drawText(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  //tex->SetTextFont(42);
  tex->SetNDC();
  tex->Draw();
};

void FillSmearedResponse(double jetPt, double truthPt, TH2D *hSmeared, TH1D *hDeltaPt) {
    double mean_dPt = 0.;
    for (int i=1; i<=hDeltaPt->GetNbinsX(); ++i) {
        if (hDeltaPt->GetBinContent(i)==0) { continue; }
        double delta_pt = hDeltaPt->GetBinCenter(i);
        double weight = hDeltaPt->GetBinContent(i);
        mean_dPt += (delta_pt)*weight;
    }
    hSmeared->Fill(truthPt,jetPt+mean_dPt);
};

TH2D* smear_TH2_by_TH1(const TH2D* hg_in, TH1D* pdf, const string title="smeared") {
/* #include <cmath> */
    // Assumptions:
    //      x-axis bin_widths in hg_in and pdf are constant and equal to each other
    //      x-axis 0 occurse in the middle of a pdf bin, and is contained in the range of pdf
    pdf->Scale(1./pdf->Integral());
    auto hg_out = (TH2D*) hg_in->Clone(title.c_str());
    hg_out->Reset();

    const double eps = 0.25;
    const TAxis* pdf_axis = pdf->GetXaxis();
    double loc0 = pdf_axis->GetBinCenter(1) / pdf_axis->GetBinWidth(1); // must be close to whole number
    loc0 += (loc0<0) ? -eps : eps;
    const int bin_offset = (int) loc0;
    const int n_pdf = pdf->GetNbinsX();

    const int ncols = hg_in->GetNbinsX();
    for (int row = 1; row <= hg_in->GetNbinsY(); ++row) {
        for (int col = 1; col <= ncols; ++col) {
            double sum_col = 0;
            double sumsq_err = 0;
            for (int i_pdf = 1; i_pdf <= n_pdf; ++i_pdf) {
                const int i_from = col - (i_pdf + bin_offset - 1);
                if (i_from < 1 || i_from > ncols) continue;
                sum_col += hg_in->GetBinContent(i_from,row) * pdf->GetBinContent(i_pdf);
                sumsq_err += TMath::Sq(hg_in->GetBinError(i_from,row) * pdf->GetBinContent(i_pdf));
            }
            if (sum_col > 0) {
                hg_out->SetBinContent(col,row,sum_col);
                hg_out->SetBinError  (col,row, TMath::Sqrt(sumsq_err));
            }
        }
    }
    hg_out->SetTitle(";HIJING BG-smeared p_{T}^{calo};p_{T}^{truth}");
    return hg_out;
};

void normalize_truth_axis( TH2D *hResponse, TH1D *hPyTruth ){
    TH1D *hSmeared_truthPt = (TH1D*) hResponse->ProjectionY("hSmeared_truthPt",1,-1,"E");
    
    TH1D *h_ratio = (TH1D*) hPyTruth->Clone("h_ratio");
    h_ratio->Divide(hSmeared_truthPt);
    for (int ix=1; ix<hResponse->GetNbinsX(); ++ix){
        for (int iy=1; iy<hResponse->GetNbinsY(); ++iy){
            if (h_ratio->GetBinContent(iy)>0. && h_ratio->GetBinContent(iy)<10000) {
                double newContent = hResponse->GetBinContent(ix,iy) * h_ratio->GetBinContent(iy);
                double newError = hResponse->GetBinError(ix,iy) * h_ratio->GetBinError(iy);
                hResponse->SetBinContent(ix,iy,newContent);
                hResponse->SetBinError(ix,iy,newError);
            }
        }
    }
    hSmeared_truthPt->Delete();
    h_ratio->Delete();
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
    
    auto *c1 = new TCanvas("c1","",700,500);

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

        if (h_CaloJetPt->at(0)<30.) { continue; }

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

    TH2D *hResp_noEmbed = new TH2D("hResp_noEmbed",";p_{T}^{calo};p_{T}^{truth}",100,0.,100.,100,0.,100.);
    TH2D *hResp_noEmbed_smear[nbins_centrality];
    
    for (int i=0; i<nbins_centrality; ++i) {
        name = "hResp_noEmbed_smear" + name_centrality[i];
        hResp_noEmbed_smear[i] = new TH2D(name,";HIJING BG-smeared p_{T}^{calo};p_{T}^{truth}",100,0.,100.,100,0.,100.);
    }

    for (int ientry=0; ientry<t_noEmb->GetEntries(); ++ientry) {
        
        t_noEmb->GetEntry(ientry);

        if (TruthJetPt->size()<=0 || p_CaloJetPt->size()<=0 || TruthJetEta->size()<=0 || p_CaloJetEta->size()<=0 || TruthJetPhi->size()<=0 || p_CaloJetPhi->size()<=0 || TruthJetPt->at(0)<=0 ) { continue; }

        if (TruthJetPt->at(0)<30.) { continue; }
        
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
//                for (int i=0; i<nbins_centrality; ++i) { FillSmearedResponse(jetPt,truthPt,hResp_noEmbed_smear[i],h_delta_pt[i]); }
                match_truth = true;
            }
        }

        
    } // end PYTHIA loop
    
    for (int i=0; i<nbins_centrality; ++i) {
        hResp_noEmbed_smear[i] = smear_TH2_by_TH1(hResp_noEmbed, h_delta_pt[i], Form("%s",hResp_noEmbed_smear[i]->GetName()));
    }
    
//    TH1D *yproj = (TH1D*) hResp_noEmbed->ProjectionY("yproj",1,hResp_noEmbed->GetNbinsY(),"E");
//    TH1D *xproj = (TH1D*) hResp_noEmbed->ProjectionX("xproj",1,hResp_noEmbed->GetNbinsX(),"E");
//    for (int ix=1; ix<=hResp_noEmbed->GetNbinsX(); ++ix) {
//        for (int iy=1; iy<=hResp_noEmbed->GetNbinsY(); ++iy) {
//            if (yproj->GetBinContent(iy)==0) { continue; }
//            hResp_noEmbed->SetBinContent(ix, iy, hResp_noEmbed->GetBinContent(ix,iy)/yproj->GetBinContent(iy) );
//        }
//    }
    
    TCanvas *c[nbins_centrality];

    for (int i=0; i<nbins_centrality; ++i) {
        
//        TH1D *yproj = (TH1D*) hResp_noEmbed_smear[i]->ProjectionY("yproj",1,hResp_noEmbed_smear[i]->GetNbinsY(),"E");
//        for (int ix=1; ix<=hResp_noEmbed_smear[i]->GetNbinsX(); ++ix) {
//            for (int iy=1; iy<=hResp_noEmbed_smear[i]->GetNbinsY(); ++iy) {
//                if (yproj->GetBinContent(iy)==0) { continue; }
//                hResp_noEmbed_smear[i]->SetBinContent(ix, iy, hResp_noEmbed_smear[i]->GetBinContent(ix,iy)/yproj->GetBinContent(iy) );
//            }
//        }
        
        name = "c" + name_centrality[i];
        c[i] = new TCanvas(name);
        c[i]->SetLogz();
        hResp_noEmbed_smear[i]->Draw("COLZ");
    }
    
    
    
    
    
    
    //             PYTHIA JET WITH EMBED (CALO RHO)
    TFile *embFile = new TFile("out/CaloJetRho_Aug16.root","READ"); // "p" corresponds to PYTHIA (or noEmbed)
    TTree *t_caloRho = (TTree*)embFile->Get("T");

    int c_id;
    float c_rho, c_rho_sigma, c_centrality, c_impactparam;
    vector<float> *c_CaloJetEta = NULL;
    vector<float> *c_CaloJetPhi = NULL;
    vector<float> *c_CaloJetE = NULL;
    vector<float> *c_CaloJetPt = NULL;
    vector<float> *c_CaloJetArea = NULL;
    TruthJetEta = NULL;
    TruthJetPhi = NULL;
    TruthJetE = NULL;
    TruthJetPt = NULL;
    TruthJetArea = NULL;

    t_caloRho->SetBranchAddress("id",&c_id);
    t_caloRho->SetBranchAddress("rho",&c_rho);
    t_caloRho->SetBranchAddress("rho_sigma",&c_rho_sigma);
    t_caloRho->SetBranchAddress("centrality",&c_centrality);
    t_caloRho->SetBranchAddress("impactparam",&c_impactparam);
    t_caloRho->SetBranchAddress("CaloJetEta",&c_CaloJetEta);
    t_caloRho->SetBranchAddress("CaloJetPhi",&c_CaloJetPhi);
    t_caloRho->SetBranchAddress("CaloJetE",&c_CaloJetE);
    t_caloRho->SetBranchAddress("CaloJetPt",&c_CaloJetPt);
    t_caloRho->SetBranchAddress("CaloJetArea",&c_CaloJetArea);
    t_caloRho->SetBranchAddress("TruthJetEta",&TruthJetEta);
    t_caloRho->SetBranchAddress("TruthJetPhi",&TruthJetPhi);
    t_caloRho->SetBranchAddress("TruthJetE",&TruthJetE);
    t_caloRho->SetBranchAddress("TruthJetPt",&TruthJetPt);
    t_caloRho->SetBranchAddress("TruthJetArea",&TruthJetArea);

    TH2D *hResp_caloRho[nbins_centrality];
    TH2D *hResp_caloJet[nbins_centrality];

    for (int i=0; i<nbins_centrality; ++i) {
        name = "hResp_caloRho" + name_centrality[i];
        hResp_caloRho[i] = new TH2D(name,";p_{T}^{reco} from embedding [GeV];p_{T}^{truth} [GeV]",100,0.,100.,100,0.,100.);
        name = "hResp_caloJet" + name_centrality[i];
        hResp_caloJet[i] = new TH2D(name,";p_{T}^{reco} - #rho #cdot A [GeV];p_{T}^{truth} [GeV]",100,0.,100.,100,0.,100.);
    }

    for (int ientry=0; ientry<t_caloRho->GetEntries(); ++ientry) {
        
        t_caloRho->GetEntry(ientry);

        if (TruthJetPt->size()<=0 || c_CaloJetPt->size()<=0 || TruthJetEta->size()<=0 || c_CaloJetEta->size()<=0 || TruthJetPhi->size()<=0 || c_CaloJetPhi->size()<=0 || TruthJetPt->at(0)<=0 ) { continue; }

        if (TruthJetPt->at(0)<30.) { continue; }

        int cent = get_centrality_bin( c_centrality );

        double truthPt = TruthJetPt->at(0);
        double truthEta = TruthJetEta->at(0);
        double truthPhi = TruthJetPhi->at(0); //        double truthE = TruthJetE->at(0);
                
        bool match_truth = false;
        
        for (int i=0; i<c_CaloJetPt->size(); ++i) {

            if (match_truth) { continue; }

            double jetPt = c_CaloJetPt->at(0);
            double jetEta = c_CaloJetEta->at(0);
            double jetPhi = c_CaloJetPhi->at(0); //        double jetE = c_CaloJetE->at(0);
            double jetArea = c_CaloJetArea->at(0);

            if (match_jet(truthEta,truthPhi,jetEta,jetPhi,0.4)) {
                hResp_caloRho[cent]->Fill(jetPt,truthPt);
                double correctedPt = jetPt - c_rho*jetArea;
                hResp_caloJet[cent]->Fill(correctedPt, truthPt);
                match_truth = true;
            }
        }

        
    } // end PYTHIA embed loop
    

// now, hResp_caloJet[i] and hResp_noEmbed_smear[i] should be *roughly* similar response matrices
// normalize to PYTHIA truth pT spectrum
    
    TH1D *hPythiaTruth = (TH1D*) hResp_noEmbed->ProjectionY("hPythiaTruth",1,-1,"E");
//    hPythiaTruth->GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN_{jets}}{dp_{T}^{truth}} [GeV^{-1}]");
    hPythiaTruth->SetTitle(";p_{T}^{truth} [GeV];1/N_{jets} dN_{jets}/dp_{T}^{truth} [GeV^{-1}]");

    for (int i=0; i<nbins_centrality; ++i) {
        normalize_truth_axis( hResp_noEmbed_smear[i], hPythiaTruth );
        normalize_truth_axis( hResp_caloJet[i], hPythiaTruth );
    }
    

    // PROJECT HISTOS HERE
        
    TH1D *hTruthPt_smear[nbins_centrality][nbins_pt];
    TH1D *hTruthPt_embed[nbins_centrality][nbins_pt];
    TH1D *hTruthPt[nbins_pt];

    for (int p=0; p<nbins_pt; ++p) {

        int binlo = hResp_noEmbed->GetXaxis()->FindBin(bins_pt[p] + 0.5);
        int binhi = hResp_noEmbed->GetXaxis()->FindBin(bins_pt[p+1] + 0.5) - 1;

        name = "hTruthPt" + name_pt[p];
        hTruthPt[p] = (TH1D*) hResp_noEmbed->ProjectionY(name,binlo,binhi,"E");
        hTruthPt[p]->Scale(1./hTruthPt[p]->Integral());
        hTruthPt[p]->SetMarkerStyle(marker_pt[p]);
        hTruthPt[p]->SetMarkerColor(kBlack);
        hTruthPt[p]->SetLineColor(kBlack);

        for (int i=0; i<nbins_centrality; ++i) {

            name = "hTruthPt_smear" + name_pt[p] + "_" + name_centrality[i] + "percent";
            hTruthPt_smear[i][p] = (TH1D*) hResp_noEmbed_smear[i]->ProjectionY(name,binlo,binhi,"E");
            name = "hTruthPt_embed" + name_pt[p] + "_" + name_centrality[i] + "percent";
            hTruthPt_embed[i][p] = (TH1D*) hResp_caloJet[i]->ProjectionY(name,binlo,binhi,"E");
        }
    }


    TH1D *hTruthPt_ratio[nbins_centrality][nbins_pt];
    for (int i=0; i<nbins_centrality; ++i) {
        for (int p=0; p<nbins_pt; ++p) {
            hTruthPt_smear[i][p]->Scale(1./hTruthPt_smear[i][p]->Integral());
            hTruthPt_embed[i][p]->Scale(1./hTruthPt_embed[i][p]->Integral());
            hTruthPt_ratio[i][p] = (TH1D*)hTruthPt_smear[i][p]->Clone();
            name = "hTruthPt_ratio" + name_pt[p] + "_" + name_centrality[i] + "percent";
            hTruthPt_ratio[i][p]->SetName(name);
            hTruthPt_ratio[i][p]->Divide(hTruthPt_embed[i][p]);
            hTruthPt_ratio[i][p]->SetTitle(";p_{T}^{truth};smeared / embedded");
            hTruthPt_ratio[i][p]->SetLineColor(color_centrality[i]);
            hTruthPt_ratio[i][p]->SetMarkerColor(color_centrality[i]);
            hTruthPt_ratio[i][p]->SetMarkerStyle(marker_pt[p]);
            hTruthPt_smear[i][p]->SetLineColor(color_centrality[i]);
            hTruthPt_smear[i][p]->SetMarkerColor(color_centrality[i]);
            hTruthPt_smear[i][p]->SetMarkerStyle(marker_pt[1]);
            hTruthPt_embed[i][p]->SetLineColor(color_centrality[i]);
            hTruthPt_embed[i][p]->SetMarkerColor(color_centrality[i]);
            hTruthPt_embed[i][p]->SetMarkerStyle(marker_pt[0]);
        }
    }

    
    auto *c0 = new TCanvas("c0","",700,500);
    for (int i=0; i<nbins_centrality; ++i) {
        hTruthPt_ratio[i][0]->Draw("pSAME");
    }
    c0->BuildLegend();

    auto *c2 = new TCanvas("c2","",700,500);
    for (int i=0; i<nbins_centrality; ++i) {
        hTruthPt_ratio[i][1]->Draw("pSAME");
    }
    c2->BuildLegend();
    

    
    auto *c3 = new TCanvas("c3","",700,500);

    for (int p=0; p<nbins_pt; ++p) {
        hTruthPt[p]->Draw("P");
        for (int i=0; i<nbins_centrality; ++i) {
            hTruthPt[p]->SetAxisRange(0.,.1,"Y");
//            hTruthPt[p]->Draw("PSAME");
//            if (i==0) { hTruthPt[p]->Draw("P"); }
            hTruthPt_smear[i][p]->Draw("PSAME");
            hTruthPt_embed[i][p]->Draw("PSAME");
        }
        c3->BuildLegend();
        name = "plots/TruthPt" + name_pt[p] + ".pdf";
        c3->SaveAs(name,"PDF");
    }
    
    
    for (int p=0; p<nbins_pt; ++p) {
        hTruthPt[p]->SetLineColor(kBlack);
        hTruthPt[p]->SetMarkerColor(kBlack);
        hTruthPt[p]->SetMarkerStyle(24);
        for (int i=0; i<nbins_centrality; ++i) {
            hTruthPt_smear[i][p]->SetLineColor(9);
            hTruthPt_smear[i][p]->SetMarkerColor(9);
            hTruthPt_smear[i][p]->SetMarkerStyle(29);
            hTruthPt_embed[i][p]->SetLineColor(8);
            hTruthPt_embed[i][p]->SetMarkerColor(8);
            hTruthPt_embed[i][p]->SetMarkerStyle(33);
        }
    }
    for (int p=0; p<nbins_pt; ++p) {
        for (int i=0; i<nbins_centrality; ++i) {
            hTruthPt[p]->Draw("P");
            hTruthPt_smear[i][p]->Draw("PSAME");
            hTruthPt_embed[i][p]->Draw("PSAME");
            c3->BuildLegend(.6,.7,.98,.98);
            name = "plots/TruthPt" + name_pt[p] + name_centrality[i] + ".pdf";
            c3->SaveAs(name,"PDF");
        }
    }
    
    
    
    TFile *outFile = new TFile("analyze/SmearJetResponse.root","RECREATE");
    for (int i=0; i<nbins_centrality; ++i) { h_delta_pt[i]->Write(); }
    hResp_noEmbed->Write();
    for (int i=0; i<nbins_centrality; ++i) {
        hResp_noEmbed_smear[i]->Write();
    }
    for (int i=0; i<nbins_centrality; ++i) {
        hResp_caloJet[i]->Write();
    }
    for (int p=0; p<nbins_pt; ++p) {
        for (int i=0; i<nbins_centrality; ++i) {
            hTruthPt_embed[i][p]->Write();
            hTruthPt_smear[i][p]->Write();
            hTruthPt_ratio[i][p]->Write();
        }
    }
    for (int p=0; p<nbins_pt; ++p) { hTruthPt[p]->Write(); }
    
    
    
    
    
    // ratio plots
    
    auto *c4 = new TCanvas("c4","",500,500);
    c4->SetBottomMargin(0.2);
    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetTopMargin(0.025);
    pad1->SetBottomMargin(0);
    pad1->SetRightMargin(0.025);
    pad1->SetLeftMargin(0.12);
    pad1->Draw();
    
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0);
    pad2->SetRightMargin(0.025);
    pad2->SetBottomMargin(0.3);
    pad2->SetLeftMargin(0.12);
    pad2->Draw();


    for (int p=0; p<1; ++p) {
        hTruthPt[p]->SetTitle(";p_{T}^{truth} [GeV];1/N_{jets} dN_{jets}/dp_{T}^{truth} [GeV^{-1}]");
//        for (int i : { 0, 4, 9 }) {
        for (int i=0; i<nbins_centrality; ++i) {

            pad1->cd();
            hTruthPt[p]->SetAxisRange(0.000001,1.,"Y");
            hTruthPt[p]->SetStats(0);
            hTruthPt[p]->GetYaxis()->SetLabelSize(.03);
            hTruthPt[p]->GetYaxis()->SetTitleSize(.05);
            hTruthPt[p]->GetYaxis()->SetTitleOffset(1.05);
            hTruthPt[p]->Draw("P");
            hTruthPt_smear[i][p]->Draw("PSAME");
            hTruthPt_embed[i][p]->Draw("PSAME");
            
            cout<<hTruthPt[p]->Integral()<<" \t"<<hTruthPt_smear[i][p]->Integral()<<" \t"<<hTruthPt_embed[i][p]->Integral()<<endl;

            title = hTruthPt[p]->GetTitle();
            hTruthPt[p]->SetTitle("pythia truth");
            hTruthPt_smear[i][p]->SetTitle("smeared");
            hTruthPt_embed[i][p]->SetTitle("embedded");
            auto leg = (TLegend*)pad1->BuildLegend(.58,.34,.86,.63);
            leg->SetFillColorAlpha(1,0.);
            leg->SetLineColorAlpha(1,0.);
            hTruthPt[p]->SetTitle(title);

            drawText(title_pt[p].c_str(), .6, 0.8, 20);
            drawText(title_centrality[i].c_str(), .65, 0.72, 20);

            c4->cd();

            pad2->cd();
            hTruthPt_ratio[i][p]->SetAxisRange(0.8,1.2,"Y");
            hTruthPt_ratio[i][p]->GetXaxis()->SetTitleSize(.1);
            hTruthPt_ratio[i][p]->GetXaxis()->SetTitleOffset(1.);
            hTruthPt_ratio[i][p]->GetXaxis()->SetLabelSize(.08);
            hTruthPt_ratio[i][p]->GetYaxis()->SetTitleSize(.1);
            hTruthPt_ratio[i][p]->GetYaxis()->SetTitleOffset(.4);
            hTruthPt_ratio[i][p]->GetYaxis()->SetLabelSize(.08);
            hTruthPt_ratio[i][p]->SetLineColor(kBlack);
            hTruthPt_ratio[i][p]->SetMarkerColor(kBlack);
            hTruthPt_ratio[i][p]->SetMarkerSize(0.7);
            hTruthPt_ratio[i][p]->SetStats(0);
            hTruthPt_ratio[i][p]->Draw("P");

            TF1 *unity_line = new TF1("unity_line","1.",0,100);
            unity_line->SetLineStyle(3);
            unity_line->SetLineColor(kBlack);
            unity_line->Draw("SAME");

            c4->cd();
      
            name = "plots/TruthPtSpectraWithRatio" + name_pt[p] + name_centrality[i] + ".pdf";
            c4->SaveAs(name,"PDF");
            
            pad1->SetLogy();
            name = "plots/TruthPtSpectraWithRatio" + name_pt[p] + name_centrality[i] + "_logy.pdf";
            leg->SetX1NDC(.22);
            leg->SetY1NDC(.15);
            leg->SetX2NDC(.50);
            leg->SetY2NDC(.41);
            c4->SaveAs(name,"PDF");
        }
    }
    
    
    
}
