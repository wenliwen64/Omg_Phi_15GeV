#include "StrAnalyMaker.hh"
#include <TMath.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLine.h>
#include <TProfile.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <iostream>
#include <string>
#include <fstream>

ClassImp(StrAnalyMaker)
StrAnalyMaker::StrAnalyMaker(std::string par_type):mParticleType(par_type), pdgmass_xi(1.67245), mKCentBin(2), mKPtBin(6){
    mCentString[0] = "10%-60%";
    mCentString[1] = "0%-10%";

    mPtBd[0] = 0.7;
    mPtBd[1] = 1.2;
    mPtBd[2] = 1.6;
    mPtBd[3] = 2.0;
    mPtBd[4] = 2.4;
    mPtBd[5] = 2.8;
    mPtBd[6] = 3.6;

    mXRawSpectra[0] = 0.95;
    mXRawSpectra[1] = 1.40;
    mXRawSpectra[2] = 1.80;
    mXRawSpectra[3] = 2.20;
    mXRawSpectra[4] = 2.60;
    mXRawSpectra[5] = 3.20;

    mXRawSpectraError[0] = 0;
    mXRawSpectraError[1] = 0;
    mXRawSpectraError[2] = 0;
    mXRawSpectraError[3] = 0;
    mXRawSpectraError[4] = 0;
    mXRawSpectraError[5] = 0;

    mDptSpectra[0] = 0.5;
    mDptSpectra[1] = 0.4;
    mDptSpectra[2] = 0.4;
    mDptSpectra[3] = 0.4;
    mDptSpectra[4] = 0.4;
    mDptSpectra[5] = 0.8;

    mXCorrSpectra[0][0] = 0.95;
    mXCorrSpectra[0][1] = 1.40;
    mXCorrSpectra[0][2] = 1.80;
    mXCorrSpectra[0][3] = 2.20;
    mXCorrSpectra[0][4] = 2.60;
    mXCorrSpectra[0][5] = 3.20;

    mXCorrSpectra[1][0] = 0.95;
    mXCorrSpectra[1][1] = 1.40;
    mXCorrSpectra[1][2] = 1.80;
    mXCorrSpectra[1][3] = 2.20;
    mXCorrSpectra[1][4] = 2.60;
    mXCorrSpectra[1][5] = 3.20;

// For 11GeV data comparison
    mDptSpectra11GeV[0] = 0.4;
    mDptSpectra11GeV[1] = 0.4;
    mDptSpectra11GeV[2] = 0.4;
    mDptSpectra11GeV[3] = 0.4;
    mDptSpectra11GeV[4] = 0.4;
    mDptSpectra11GeV[5] = 0.8;

    mXCorrSpectra11GeV[0][0] = 1.0113;
    mXCorrSpectra11GeV[0][1] = 1.39489;
    mXCorrSpectra11GeV[0][2] = 1.78865;
    mXCorrSpectra11GeV[0][3] = 2.18517;
    mXCorrSpectra11GeV[0][4] = 2.58297;
    mXCorrSpectra11GeV[0][5] = 3.12641;

    mXCorrSpectra11GeV[1][0] = 1.01367;
    mXCorrSpectra11GeV[1][1] = 1.3961;
    mXCorrSpectra11GeV[1][2] = 1.7898;
    mXCorrSpectra11GeV[1][3] = 2.18636;
    mXCorrSpectra11GeV[1][4] = 2.58419;
    mXCorrSpectra11GeV[1][5] = 3.13102;

    mYCorrSpectra11GeV[0][0] = 0.000970335;
    mYCorrSpectra11GeV[0][1] = 0.000395753;
    mYCorrSpectra11GeV[0][2] = 0.000234294;
    mYCorrSpectra11GeV[0][3] = 6.37403e-05;
    mYCorrSpectra11GeV[0][4] = 1.01077e-05;
    mYCorrSpectra11GeV[0][5] = 2.61219e-06;

    mYCorrSpectra11GeV[1][0] = 0.00367224;
    mYCorrSpectra11GeV[1][1] = 0.00140577;
    mYCorrSpectra11GeV[1][2] = 0.000653242;
    mYCorrSpectra11GeV[1][3] = 0.000262759;
    mYCorrSpectra11GeV[1][4] = 7.855e-05;
    mYCorrSpectra11GeV[1][5] = 8.46299e-06;

    //Statistical Error Only
    mYCorrSpectra11GeVError[0][0] = 0.000211822;
    mYCorrSpectra11GeVError[0][1] = 5.50952e-05;
    mYCorrSpectra11GeVError[0][2] = 2.50402e-05;
    mYCorrSpectra11GeVError[0][3] = 9.66436e-06;
    mYCorrSpectra11GeVError[0][4] = 3.27398e-06;
    mYCorrSpectra11GeVError[0][5] = 9.57011e-07;

    mYCorrSpectra11GeVError[1][0] = 0.00148974;
    mYCorrSpectra11GeVError[1][1] = 0.00031946;
    mYCorrSpectra11GeVError[1][2] = 0.0001118;
    mYCorrSpectra11GeVError[1][3] = 5.20221e-05;
    mYCorrSpectra11GeVError[1][4] = 2.23603e-05;
    mYCorrSpectra11GeVError[1][5] = 4.5805e-06;

    std::cout << "StrAnalyMaker Constructor v0.01 2015-08-30 " << std::endl;
}

StrAnalyMaker::~StrAnalyMaker(){}

void StrAnalyMaker::Init(std::string overview_filename, std::string dat_filename, std::string rotbg_filename, std::string fpeff_filename, std::string expeff_filename){
    std::cout << "!!! InitializationII for AuAu14.5GeV " << mParticleType << " Analysis" << std::endl;
    // Initialize TFile pointers 
    mOverviewFile = new TFile(overview_filename.c_str(), "read");
    if(mOverviewFile->IsZombie()){
        std::cout << "ERROR!!! No overview file!" << std::endl;
	exit(1);
    }
    mDatFile = new TFile(dat_filename.c_str(), "read");
    if(mDatFile->IsZombie()){
    std::cout << "ERROR!!! No data file!" << std::endl;
	exit(1);
    }
    mRotBgFile = new TFile(rotbg_filename.c_str(), "read");
    if(mRotBgFile->IsZombie()){
	std::cout << "ERROR!!! No rotational file!" << std::endl;
	exit(1);
    }
    mFpEffFile = new TFile(fpeff_filename.c_str(), "read");
    if(mFpEffFile->IsZombie()){
	std::cout << "ERROR!!! No fpeff file!" << std::endl;
	exit(1);
    }
    mExpEffFile = new TFile(expeff_filename.c_str(), "read");
    if(mExpEffFile->IsZombie()){
	std::cout << "ERROR!!! No expeff file!" << std::endl;
	exit(1);
    }


    // Get initial efficiency
    effInit();

    // Initialize Levy Function
    levyInit(); 

    // Get the event number of weighted and unweighted 
    nEventsInit();

    std::cout << "happy!" << std::endl;
    // Initialize branching ratio, Lambda->p+pi * Omega->Lambda+k
    mBr = 0.678*0.639;

    // Initialize signal counting range
    mSigRangeLeft = 1.66;
    mSigRangeRight = 1.685;

    //Initialize rotational background parameters
    rotBgAnalysisInit();
}

void StrAnalyMaker::effInit(){
    for(int i = 0; i < mKCentBin; i++){
	char originEffName[50];
	sprintf(originEffName, "effxiptcent%d", i);
	mHFpEffFine[i] = (TH1F*)mFpEffFile->Get(originEffName); 
	mHExpEffFine[i] = (TH1F*)mExpEffFile->Get(originEffName);
	for(int j = 0; j < mKPtBin; j++){
	    mFpEff[i][j] = 1.0;
	    mFpEffError[i][j] = 0.0;
	    mExpEff[i][j] = 1.0;
	    mExpEffError[i][j] = 0.0;
	    mEff[i][j] = 1.0;
	    mEffError[i][j] = 0.0;
	}
    }
}

void StrAnalyMaker::levyInit(){
    mLevy = new TF1("levy", "[0]*pow(1+(sqrt(x*x+1.67245*1.67245)-1.67245)/([1]*[2]),-[1])*([1]-1)*([1]-2)/(2*3.14159265*[1]*[2]*([1]*[2]+1.67245*([1]-2)))", 0., 8.);
    mLevyPt = new TF1("levyPt", "x*[0]*pow(1+(sqrt(x*x+1.67245*1.67245)-1.67245)/([1]*[2]),-[1])*([1]-1)*([1]-2)/(2*3.14159265*[1]*[2]*([1]*[2]+1.67245*([1]-2)))", 0., 8.);
    mLevyPt2 = new TF1("levyPt2", "x*x*[0]*pow(1+(sqrt(x*x+1.67245*1.67245)-1.67245)/([1]*[2]),-[1])*([1]-1)*([1]-2)/(2*3.14159265*[1]*[2]*([1]*[2]+1.67245*([1]-2)))", 0., 8.);

    mLevy->SetParName(0, "dN/dy");
    mLevy->SetParName(1, "a");
    mLevy->SetParName(2, "T");
    mLevy->SetLineColor(4);
    mLevy->SetLineStyle(2);
    mLevyPar[0][0] = 0.01;   
    mLevyPar[1][0] = 0.04;

    mLevyPar[0][1] = 6.0e+06;
    mLevyPar[1][1] = 6.1e+06;

    mLevyPar[0][2] = 0.22;
    mLevyPar[1][2] = 0.29;
}

//void StrAnalyMaker::expInit(){
    //mExp = new TF1("exp", );
//}

void StrAnalyMaker::nEventsInit(){
/*
    TH1F* h_centbin9_unweighted = (TH1F*)mOverviewFile->Get("h_centbin9_after0");
    TH1F* h_centbin9_weighted = (TH1F*)mOverviewFile->Get("h_centbin9_after1");

    mNEventsUnweighted[0] = h_centbin9_unweighted->GetBinContent(8) + h_centbin9_unweighted->GetBinContent(7) + h_centbin9_unweighted->GetBinContent(6) + h_centbin9_unweighted->GetBinContent(5) + h_centbin9_unweighted->GetBinContent(4);
    mNEventsWeighted[0] = h_centbin9_weighted->GetBinContent(8) + h_centbin9_weighted->GetBinContent(7) + h_centbin9_weighted->GetBinContent(6) + h_centbin9_weighted->GetBinContent(5) + h_centbin9_weighted->GetBinContent(4);
    mNEventsUnweighted[1] = h_centbin9_unweighted->GetBinContent(10) + h_centbin9_unweighted->GetBinContent(9);
    mNEventsWeighted[1] = h_centbin9_weighted->GetBinContent(10) + h_centbin9_weighted->GetBinContent(9);

    for(int i = 0; i < mKCentBin; i++){
	std::cout << mNEventsUnweighted[i] << "cent" << i << " nevents unweighted!" << std::endl;
	std::cout << mNEventsWeighted[i] << "cent" << i << " nevents weighted!" << std::endl;
    }
*///TODO: should get the overview soon
   TH1F* h_centbin9_unweighted = (TH1F*)mOverviewFile->Get("h_centbin9_after1");
   TH1F* h_centbin9_weighted = (TH1F*)mOverviewFile->Get("h_centbin9_after1");

    mNEventsUnweighted[0] = h_centbin9_unweighted->GetBinContent(8) + h_centbin9_unweighted->GetBinContent(7) + h_centbin9_unweighted->GetBinContent(6) + h_centbin9_unweighted->GetBinContent(5) + h_centbin9_unweighted->GetBinContent(4);
    mNEventsWeighted[0] = h_centbin9_weighted->GetBinContent(8) + h_centbin9_weighted->GetBinContent(7) + h_centbin9_weighted->GetBinContent(6) + h_centbin9_weighted->GetBinContent(5) + h_centbin9_weighted->GetBinContent(4);
    mNEventsUnweighted[1] = h_centbin9_unweighted->GetBinContent(10) + h_centbin9_unweighted->GetBinContent(9);
    mNEventsWeighted[1] = h_centbin9_weighted->GetBinContent(10) + h_centbin9_weighted->GetBinContent(9);

    for(int i = 0; i < mKCentBin; i++){
	std::cout << mNEventsUnweighted[i] << "cent" << i << " nevents unweighted!" << std::endl;
	std::cout << mNEventsWeighted[i] << "cent" << i << " nevents weighted!" << std::endl;
    }

}

void StrAnalyMaker::rotBgAnalysisInit(){
    std::cout << "!!! Analyze Rotational Background Initialization" << std::endl;
    // Initialize normalization range
    mRotNormLeftLowB = 1.625;//pdgmass_xi - 0.05;
    mRotNormLeftHighB = 1.655;//pdgmass_xi - 0.015;

    mRotNormRightLowB = 1.69;//pdgmass_xi + 0.015;
    mRotNormRightHighB = 1.72;//pdgmass_xi + 0.05;
}

Double_t StrAnalyMaker::compRotNormFactor(Int_t centbin, Int_t ptbin,  TH1F* hdat, TH1F* hrot){
    std::cout << "!!! Compute Rot Norm Factor!" << std::endl;

    Int_t ratio_l1 = hrot->FindBin(mRotNormLeftLowB);
    Int_t ratio_l2 = hrot->FindBin(mRotNormRightLowB);
    Int_t ratio_u1 = hrot->FindBin(mRotNormLeftHighB);
    Int_t ratio_u2 = hrot->FindBin(mRotNormRightHighB);
 
    mRotScale_ratio[centbin][ptbin] = (hrot->Integral(ratio_l1, ratio_u1) + hrot->Integral(ratio_l2, ratio_u2)) / (hdat->Integral(ratio_l1, ratio_u1) + hdat->Integral(ratio_l2, ratio_u2));    
    mRotScale_ratio[centbin][ptbin] = 1.; // In low energies, use 1 as the normalization factor
    std::cout << "------Norm Factor For cent" << centbin << "pt" << ptbin << "is " << mRotScale_ratio[centbin][ptbin] << std::endl;
    return mRotScale_ratio[centbin][ptbin];
}

void StrAnalyMaker::plotRotInvMassWithData(Int_t centbin, Int_t ptbin, TH1F* hdat, TH1F* hrot, Double_t scale){
    
    std::cout << "!!! Plot Inv Mass" << std::endl;
    hdat->SetMarkerStyle(8);
    hdat->Draw("Hist");
    hdat->GetXaxis()->SetTitle("InvMass(GeV)");
    hdat->GetYaxis()->SetTitle("Counts");

    TH1F* hrot_copy = (TH1F*)hrot->Clone();
    hrot_copy->Sumw2();
    hrot_copy->Scale(1./scale);
    hrot_copy->SetLineColor(2);
    hrot_copy->SetFillColor(2);
    hrot_copy->SetFillStyle(3354);
    gPad->SetTicks(1, 1);
    hrot_copy->Draw("Hist sames");

    TLine* lline = new TLine(mSigRangeLeft, 0, mSigRangeLeft, hdat->GetMaximum());
    TLine* uline  = new TLine(mSigRangeRight, 0, mSigRangeRight, hdat->GetMaximum());
    lline->SetLineColor(4);
    lline->SetLineWidth(2);
    lline->SetLineStyle(10);
    uline->SetLineColor(4);
    uline->SetLineWidth(2);
    uline->SetLineStyle(10);
    lline->Draw("sames");
    uline->Draw("sames");

    TLine* l1line = new TLine(mRotNormLeftLowB, 0, mRotNormLeftLowB, hdat->GetMaximum());
    TLine* u1line = new TLine(mRotNormLeftHighB, 0, mRotNormLeftHighB, hdat->GetMaximum());
    l1line->SetLineColor(6);
    l1line->SetLineWidth(2);
    l1line->SetLineStyle(10);
    u1line->SetLineColor(6);
    u1line->SetLineWidth(2);
    u1line->SetLineStyle(10);
    l1line->Draw("sames");
    u1line->Draw("sames");

    TLine* l2line = new TLine(mRotNormRightLowB, 0, mRotNormRightLowB, hdat->GetMaximum());
    TLine* u2line = new TLine(mRotNormRightHighB, 0, mRotNormRightHighB, hdat->GetMaximum());
    l2line->SetLineColor(6);
    l2line->SetLineWidth(2);
    l2line->SetLineStyle(10);
    u2line->SetLineColor(6);
    u2line->SetLineWidth(2);
    u2line->SetLineStyle(10);
    l2line->Draw("sames");
    u2line->Draw("sames"); 
    
    char plotname[50];
    sprintf(plotname, "../%s_plots/rot_%spt%dcent%d.pdf", mParticleType.c_str(), mParticleType.c_str(), ptbin+1, centbin);
    gPad->SaveAs(plotname);
}

void StrAnalyMaker::compRawSigCounts(Int_t centbin, Int_t ptbin, TH1F* hdat, TH1F* hrot, Double_t scale){
    //TH1F* hrot_copy = (TH1F*)hrot->Clone();
    Int_t sigRangeLeftBin = hdat->FindBin(mSigRangeLeft);
    Int_t sigRangeRightBin = hdat->FindBin(mSigRangeRight);
    Int_t datcounts = hdat->Integral(sigRangeLeftBin, sigRangeRightBin);
    Int_t rotcounts = hrot->Integral(sigRangeLeftBin, sigRangeRightBin);
    mRawSigCounts[centbin][ptbin] = datcounts - rotcounts/scale;
    mRawSigCountsError[centbin][ptbin] = sqrt(datcounts + (rotcounts/scale));
    std::cout << "!!! mRawSigCounts for cent " << centbin << "pt" << ptbin << " is " << mRawSigCounts[centbin][ptbin] << " error = " << mRawSigCountsError[centbin][ptbin] << std::endl;
}

void StrAnalyMaker::compRawSpectra(){
    double PI = 3.1415926; 
    for(int i = 0; i < mKCentBin; i++){
        for(int j = 0; j < mKPtBin; j++){
            mYRawSpectra[i][j] = 1/(2*PI) * mRawSigCounts[i][j] / mXRawSpectra[j] / mDptSpectra[j] / mNEventsWeighted[i] / mBr;
            mYRawSpectraError[i][j] = 1/(2*PI) * mRawSigCountsError[i][j] / mXRawSpectra[j] / mDptSpectra[j] / mNEventsWeighted[i] / mBr;
            std::cout << "mYRawSpectra for cent" << i << "pt" << j <<" is " << mYRawSpectra[i][j] << "mYRawSpectraError is " << mYRawSpectraError[i][j] << std::endl;
	}
    }
}

void StrAnalyMaker::plotRawSpectra(){ 
    TCanvas* rawspectra_can = new TCanvas("rawspectra_can", "rawspectra_can");
    rawspectra_can->SetLogy();
    rawspectra_can->SetTicks(1, 1);

    TGraphErrors* GRawSpectra_1060 = new TGraphErrors(6, mXRawSpectra, mYRawSpectra[0], mXRawSpectraError, mYRawSpectraError[0]); 
    GRawSpectra_1060->SetMarkerSize(1.0);
    GRawSpectra_1060->SetMarkerStyle(34);
    GRawSpectra_1060->SetMarkerColor(2);
    GRawSpectra_1060->SetMaximum(10E-1);
    GRawSpectra_1060->SetMinimum(10E-9);
    GRawSpectra_1060->GetXaxis()->SetLimits(0.00, 3.60);
    GRawSpectra_1060->SetTitle("#Omega^{+} Spectra, Au+Au 14.5GeV");
    GRawSpectra_1060->GetYaxis()->SetTitle("#frac{d^{2}N}{2#piNP_{T}dP_{T}dy}(GeV/c)^{-2}");
    GRawSpectra_1060->GetXaxis()->SetTitle("P_{T}(GeV/c)");
    GRawSpectra_1060->GetYaxis()->SetTitleOffset(1.3);
    GRawSpectra_1060->Draw("AP");


    TGraphErrors* GRawSpectra_010 = new TGraphErrors(6, mXRawSpectra, mYRawSpectra[1], mXRawSpectraError, mYRawSpectraError[1]); 
    GRawSpectra_010->SetMarkerSize(1.0);
    GRawSpectra_010->SetMarkerStyle(34);
    GRawSpectra_010->SetMarkerColor(1);
    GRawSpectra_010->Draw("P same");

    TLegend* leg = new TLegend(0.65, 0.65, 0.85, 0.85);
    leg->SetBorderSize(0);
    leg->AddEntry(GRawSpectra_010, "0-10%", "p");
    leg->AddEntry(GRawSpectra_1060, "10-60%", "p");
    leg->Draw("sames");

    std::string plotname = "../" + mParticleType + "_plots/" + mParticleType + "_rawspectra.pdf";
    rawspectra_can->SaveAs(plotname.c_str());  
}


void StrAnalyMaker::analyzeEff(){
    TProfile* pFpEff[2];
    TProfile* pExpEff[2];
    for(int i = 0; i < mKCentBin; i++){
        char coarseFpProfileName[50];
        char coarseExpProfileName[50];
        sprintf(coarseFpProfileName, "fpeffprofile_cent%d", i);
        sprintf(coarseExpProfileName, "expeffprofile_cent%d", i);
        pFpEff[i] = new TProfile(coarseFpProfileName, coarseFpProfileName, mKPtBin, mPtBd);
        pExpEff[i] = new TProfile(coarseExpProfileName, coarseExpProfileName, mKPtBin, mPtBd);
        pFpEff[i]->BuildOptions(0, 0, "s");
        pExpEff[i]->BuildOptions(0, 0, "s");

	Int_t noBins = mHFpEffFine[i]->GetSize() - 2;
        for(int j = 0; j < noBins; j++){
            Double_t effPt = mHFpEffFine[i]->GetBinCenter(j+1);
            Double_t fpEff = mHFpEffFine[i]->GetBinContent(j+1);
            Double_t expEff = mHExpEffFine[i]->GetBinContent(j+1);
            Double_t weight = getSpectraWeight(i, effPt);

            //std::cout << "original eff data = " << effPt << " " << fpEff << " " << expEff << " " << weight << std::endl; 
            pFpEff[i]->Fill(effPt, fpEff, weight);
            pExpEff[i]->Fill(effPt, expEff, weight);
	}
        
        for(int k = 0; k < mKPtBin; k++){
	    mExpEff[i][k] = pExpEff[i]->GetBinContent(k+1);
	    mExpEffError[i][k] = pExpEff[i]->GetBinError(k+1);
	    mFpEff[i][k] = pFpEff[i]->GetBinContent(k+1);
	    mFpEffError[i][k] = pFpEff[i]->GetBinError(k+1);
            if(k < 2){
                mEff[i][k] = pExpEff[i]->GetBinContent(k+1);
                mEffError[i][k] = pExpEff[i]->GetBinError(k+1);
            }
            else{
                mEff[i][k] = pFpEff[i]->GetBinContent(k+1);
                mEffError[i][k] = pFpEff[i]->GetBinError(k+1);
	    }
            std::cout << "!!! mEff for cent" << i << "pt" << k << " is " << mEff[i][k] << std::endl;
	}
    }

    // Plot eff error flat v.s. exp 
    Double_t relative_efferr_fp[mKCentBin][mKPtBin];
    Double_t relative_efferr_exp[mKCentBin][mKPtBin];
    for(int i = 0; i < mKCentBin; i++){
        for(int j = 0; j < mKPtBin; j++){
	    relative_efferr_fp[i][j] = mFpEffError[i][j] / mFpEff[i][j];
	    relative_efferr_exp[i][j] = mExpEffError[i][j] / mExpEff[i][j];
	    std::cout << "=========> " << relative_efferr_exp[i][j] << " -> exp; " << relative_efferr_fp[i][j] << " -> fp;"  << "<==========" << std::endl;
	}
    }

    TGraph* g_fp;
    TGraph* g_exp;
    TCanvas* c = new TCanvas();
    c->Divide(2, 1);
    for(int i = 0; i < mKCentBin; i++){
	c->cd(i+1);
	g_fp = new TGraph(mKPtBin, mXRawSpectra, relative_efferr_fp[i]);
	char title[50];
	sprintf(title, "#Omega^{-} efferr_comparison(r-exp, k-fp)cent%d", i);
	g_fp->SetTitle(title);
	g_fp->GetXaxis()->SetTitle("pT(GeV)");
	g_fp->GetYaxis()->SetTitle("eff_err");
	g_fp->SetMarkerStyle(20);
	g_fp->SetMarkerColor(1);
	g_fp->Draw("AP");

	g_exp = new TGraph(mKPtBin, mXRawSpectra, relative_efferr_exp[i]);
	g_exp->SetMarkerStyle(20);
	g_exp->SetMarkerColor(2);
	g_exp->Draw("P same");

	gPad->SetTicks(1, 1);
    } 

    char plotname[50];
    sprintf(plotname, "../%s_plots/%s_efferr_comparison.pdf", mParticleType.c_str(), mParticleType.c_str());
    c->SaveAs(plotname);

}

std::string StrAnalyMaker::getCentString(Int_t i){
    return mCentString[i];
}

void StrAnalyMaker::plotEff(){
    TCanvas* canEff = new TCanvas("canEff", "canEff");
    canEff->SetLogy();
    canEff->SetTicks(1, 1);
    TLegend* leg = new TLegend(0.55, 0.45, 0.75, 0.65);
    leg->SetBorderSize(0);
    for(int j = 0; j < mKCentBin; j++){
        Int_t i =  mKCentBin - 1 - j;
	TGraphErrors* geff = new TGraphErrors(mKPtBin, mXRawSpectra, mEff[i], 0, mEffError[i]);
        geff->SetMarkerSize(1.0);
        geff->SetMarkerStyle(20);
        geff->SetMarkerColor(i+1);
        geff->SetMaximum(0.1);
        geff->SetMinimum(0.0);
        geff->GetXaxis()->SetLimits(0.0, 3.6);
        if(i == 1){
            geff->SetTitle("#Omega^{+} Efficiency, Au+Au 14.5GeV"); 
            geff->GetYaxis()->SetTitle("efficiency");
            geff->GetXaxis()->SetTitle("P_{T}(GeV/c)");
            geff->Draw("AP");
	}
        else{
            geff->Draw("P same");
	}
        std::string legend = getCentString(i);
        leg->AddEntry(geff, legend.c_str(), "p");
    }
    leg->Draw("sames");

    char plotname[50]; 
    sprintf(plotname, "../%s_plots/final_eff_combined.pdf", mParticleType.c_str());
    canEff->SaveAs(plotname);
}

Double_t StrAnalyMaker::getSpectraWeight(Int_t centbin, Double_t pt){
    Double_t wgt;
    //std::cout << "!!! mLevyPar = " << mLevyPar[centbin][0] << " " << mLevyPar[centbin][1] << " " << mLevyPar[centbin][2] << std::endl;
    mLevyPt->SetParameters(mLevyPar[centbin]);
    wgt = mLevyPt->Eval(pt);
    return wgt;
}

void StrAnalyMaker::compCorrSpectra(){
    //Compute Data points
    for(int i = 0;  i < mKCentBin; i++){
        for(int j = 0; j < mKPtBin; j++){
            mYCorrSpectra[i][j] = mYRawSpectra[i][j] / mEff[i][j]; 
            Double_t realitive_rawyerror = mYRawSpectraError[i][j]/mYRawSpectra[i][j];
            Double_t realitive_efferror = mEffError[i][j]/mEff[i][j];
            mYCorrSpectraError[i][j] = mYCorrSpectra[i][j] * sqrt(realitive_rawyerror*realitive_rawyerror); 
            //mYCorrSpectraError[i][j] = mYCorrSpectra[i][j] * sqrt(realitive_rawyerror*realitive_rawyerror + realitive_efferror*realitive_efferror); 
            std::cout << "Calculation========YCorrSpectraCent" << i << "Pt" << j<< "= " << mYCorrSpectra[i][j] << " with error = " << mYCorrSpectraError[i][j] << " efferror =  " << mEffError[i][j] << "============" << std::endl;
	}
    } 

    //Fitting Levy Function, obtain the X position and output the fitting par
    for(int i = 0; i < mKCentBin; i++){
        TGraphErrors* g = new TGraphErrors(mKPtBin, mXCorrSpectra[i], mYCorrSpectra[i], 0, mYCorrSpectraError[i]);
        mLevy->SetParameters(mLevyPar[i]);
        if(i == 0)
	    mLevy->SetParameter(1, -6e+06);
	g->Fit(mLevy, "REM0");

	mLevyPar[i][0] = mLevy->GetParameter(0);
	mLevyParError[i][0] = mLevy->GetParError(0);
	mLevyPar[i][1] = mLevy->GetParameter(1);
	mLevyParError[i][1] = mLevy->GetParError(1);
	mLevyPar[i][2] = mLevy->GetParameter(2);
	mLevyParError[i][2] = mLevy->GetParError(2);
        mLevyPt->SetParameters(mLevyPar[i]);
        mLevyPt2->SetParameters(mLevyPar[i]);

	for(int j = 0; j < mKPtBin; j++){
	    mXCorrSpectra[i][j] = mLevyPt2->Integral(mPtBd[j], mPtBd[j+1]) / mLevyPt->Integral(mPtBd[j], mPtBd[j+1]);
	}
    } 
}

void StrAnalyMaker::compDndy(){
    // Initialize measured and unmeasured pt range
    Double_t leftlow = 0.;
    Double_t lefthigh = 0.7;
    Double_t rightlow = 3.6;
    Double_t righthigh = 15.0;
    for(int i = 0; i < mKCentBin; i++){
	mLevy->SetParameters(mLevyPar[i]); 
	mLevyPt->SetParameters(mLevyPar[i]); 
        mDndy[i] = 0;
        mDndyError[i] = 0;
        Double_t tmp_dndy_err02 = 0;
        Double_t tmp_dndy_err1 = 0;
        Double_t tmp_dndy_err2 = 0;
        for(int j = 0; j < mKPtBin; j++){
	    mDndy[i] += mYCorrSpectra[i][j]*2*3.1415926*mDptSpectra[j]*mXCorrSpectra[i][j];
            tmp_dndy_err02 += (mYCorrSpectraError[i][j]*2*3.1415926*mDptSpectra[j]*mXCorrSpectra[i][j])*(mYCorrSpectraError[i][j]*2*3.1415926*mDptSpectra[j]*mXCorrSpectra[i][j]);
	}
	mDndy[i] += 2*3.1415926*(mLevyPt->Integral(leftlow, lefthigh) + mLevyPt->Integral(rightlow, righthigh));
        tmp_dndy_err1 = mLevyParError[i][0]/mLevyPar[i][0]*2*3.1415926*mLevyPt->Integral(leftlow, lefthigh);
        tmp_dndy_err2 = mLevyParError[i][0]/mLevyPar[i][0]*2*3.1415926*mLevyPt->Integral(rightlow, righthigh);
        mDndyError[i] += sqrt(tmp_dndy_err02 + tmp_dndy_err1*tmp_dndy_err1 + tmp_dndy_err2*tmp_dndy_err2);//TODO:
    }
}

Double_t StrAnalyMaker::getDndy(Int_t centbin){
    return mDndy[centbin];
}
 
Double_t StrAnalyMaker::getDndyError(Int_t centbin){
    return mDndyError[centbin];
}

void StrAnalyMaker::plotCorrSpectra(){
    TCanvas* canCorrSpectra = new TCanvas("canCorrSpectra", "canCorrSpectra");
    canCorrSpectra->SetLogy();
    canCorrSpectra->SetTicks(1, 1);
    TLegend* leg = new TLegend(0.65, 0.65, 0.85, 0.85);
    leg->SetBorderSize(0);
    for(int j = 0; j < mKCentBin; j++){
        Int_t i = mKCentBin - j - 1;
	TGraphErrors* gerr = new TGraphErrors(mKPtBin, mXCorrSpectra[i], mYCorrSpectra[i], 0, mYCorrSpectraError[i]); 
        gerr->SetMarkerSize(1.0);
        gerr->SetMarkerStyle(20);
        gerr->SetMarkerColor(i+1);
	std::string centString = getCentString(i);
        if(i == 1){
	    gerr->SetMaximum(10E-2);
	    gerr->SetMinimum(10E-8);
	    gerr->GetXaxis()->SetLimits(0.0, 3.60);
	    gerr->SetTitle("#Omega^{+} Spectra, Au+Au 14.5GeV");
	    gerr->GetYaxis()->SetTitle("#frac{d^{2}N}{2#piNP_{T}dP_{T}dy}(GeV/c)^{-2}");
	    gerr->GetXaxis()->SetTitle("P_{T}(GeV/c)");
	    gerr->GetYaxis()->SetTitleOffset(1.3);
	    gerr->Draw("AP");
	}
        else
            gerr->Draw("P same");
        
        mLevy->SetParameters(mLevyPar[i]);
        TF1* levy_copy = (TF1*)mLevy->Clone();
        levy_copy->SetLineStyle(2);
        levy_copy->SetLineColor(4);
        levy_copy->Draw("sames");
        leg->AddEntry(gerr, centString.c_str(), "p");
    }
    leg->AddEntry(mLevy, "Levy Function", "l");
    leg->Draw("sames");

    char plotname[50]; 
    sprintf(plotname, "../%s_plots/finalCorrSpectra_%s.pdf", mParticleType.c_str(), mParticleType.c_str());
    canCorrSpectra->SaveAs(plotname);
    //canCorrSpectra->SaveAs("../omg_plots/finalCorrSpectra.gif");
    //canCorrSpectra->SaveAs("../omg_plots/finalCorrSpectra.eps");
    //canCorrSpectra->SaveAs("../omg_plots/finalCorrSpectra.jpg");
}

double IntegralError(int npar, double * c, double * errPar, 
   double * covMatrix = 0) {   
// calculate the error on the integral given the parameter error and 
// the integrals of the gradient functions c[] 

   double err2 = 0; 
   for (int i = 0; i < npar; ++i) { 
      if (covMatrix == 0) { // assume error are uncorrelated
         err2 += c[i] * c[i] * errPar[i] * errPar[i]; 
      } else {
         double s = 0; 
         for (int j = 0; j < npar; ++j) {
            s += covMatrix[i*npar + j] * c[j]; 
         }
         err2 += c[i] * s; 
      }
   }

   return TMath::Sqrt(err2);
}

void StrAnalyMaker::compYields(){
    for(int i = 0; i < mKCentBin; i++){
        for(int j = 0; j < mKPtBin; j++){
	    mYields[i][j] = mYCorrSpectra[i][j]*2*3.1415926*mDptSpectra[j]*mXCorrSpectra[i][j];
            mYieldsError[i][j] = mYCorrSpectraError[i][j]*2*3.1415926*mDptSpectra[j]*mXCorrSpectra[i][j];

	    mYields11GeV[i][j] = mYCorrSpectra11GeV[i][j]*2*3.1415926*mDptSpectra11GeV[j]*mXCorrSpectra11GeV[i][j];
            mYields11GeVError[i][j] = mYCorrSpectra11GeVError[i][j]*2*3.1415926*mDptSpectra11GeV[j]*mXCorrSpectra11GeV[i][j];
            std::cout << "yields = " << mYields[i][j] << std::endl;
	}
    }

}

void StrAnalyMaker::compare11GeV(){
    TCanvas* canCompare11GeV = new TCanvas("canCompare11GeV", "canCompare11GeV");
    canCompare11GeV->SetLogy();
    canCompare11GeV->SetTicks(1, 1);
    TLegend* leg = new TLegend(0.65, 0.65, 0.85, 0.85);
    leg->SetBorderSize(0);
    for(int j = 0; j < mKCentBin; j++){
        Int_t i = mKCentBin - j - 1;
        TGraphErrors* gerr = new TGraphErrors(mKPtBin, mXCorrSpectra[i], mYields[i], 0, mYieldsError[i]);
        gerr->SetMarkerSize(1.2); 
        gerr->SetMarkerStyle(20);
        gerr->SetMarkerColor(i+1);
        if(i == 1){
            gerr->SetMinimum(10e-8);
            gerr->SetMaximum(10e-1);
            gerr->GetXaxis()->SetLimits(0.0, 3.60);
            gerr->GetXaxis()->SetTitle("pT(GeV/c)");
            gerr->GetYaxis()->SetTitle("Yields");
            std::string title = "#Omega^{+} Yileds, Au+Au 14.5GeV";
            gerr->SetTitle(title.c_str());
            gerr->Draw("AP same");
	}
        else
            gerr->Draw("P same");

        std::string centString = getCentString(i);
        leg->AddEntry(gerr, centString.c_str(), "p");
    }    

    for(int j = 0; j < mKCentBin; j++){
        Int_t i = mKCentBin - j - 1;
        TGraphErrors* gerr = new TGraphErrors(mKPtBin, mXCorrSpectra11GeV[i], mYields11GeV[i], 0, mYields11GeVError[i]);
        gerr->SetMarkerSize(1.2); 
        gerr->SetMarkerStyle(25);
        gerr->SetMarkerColor(i+1);
	gerr->Draw("P same");

        std::string centString = getCentString(i) + "11GeV";
        leg->AddEntry(gerr, centString.c_str(), "p");
    }
    leg->Draw("same");
    
    char plotName[50];
    sprintf(plotName, "../%s_plots/compare11GeV_%s.pdf", mParticleType.c_str(), mParticleType.c_str());
    gPad->SaveAs(plotName);
}
void StrAnalyMaker::Analyze(){
    std::cout << "Load infile_dat/rot successfully!" << std::endl;
    for(int i = 0; i < mKPtBin; i++){
        char hist_name_rot_010[200]; 
        char hist_name_rot_1060[200];
        char hist_name_dat_010[200];
        char hist_name_dat_1060[200];
	sprintf(hist_name_rot_010, "sig_xipt%dcent_010", i+1);
	sprintf(hist_name_rot_1060, "sig_xipt%dcent_1060", i+1); 
	sprintf(hist_name_dat_010, "sig_xipt%dcent_010", i+1); 
	sprintf(hist_name_dat_1060, "sig_xipt%dcent_1060", i+1); 

        TH1F* hrot_010 = (TH1F*)mRotBgFile->Get(hist_name_rot_010);
        TH1F* hdat_010 = (TH1F*)mDatFile->Get(hist_name_dat_010);
        TH1F* hrot_1060 = (TH1F*)mRotBgFile->Get(hist_name_rot_1060);
        TH1F* hdat_1060 = (TH1F*)mDatFile->Get(hist_name_dat_1060);
       
	hrot_010->Rebin(4);
        hdat_010->Rebin(4);
        hrot_1060->Rebin(4);
        hdat_1060->Rebin(4);

        hrot_010->Sumw2();
        hdat_010->Sumw2();
        hrot_1060->Sumw2();
        hdat_1060->Sumw2();

        double rot_scale_1060 = compRotNormFactor(0, i, hdat_1060, hrot_1060);
        double rot_scale_010 = compRotNormFactor(1, i, hdat_010, hrot_010);
        //rot_scale_1060 = 1.;
        //rot_scale_010 = 1.;

        plotRotInvMassWithData(0, i, hdat_1060, hrot_1060, rot_scale_1060);
        plotRotInvMassWithData(1, i, hdat_010, hrot_010, rot_scale_010);
       
        compRawSigCounts(0, i, hdat_1060, hrot_1060, rot_scale_1060); 
        compRawSigCounts(1, i, hdat_010, hrot_010, rot_scale_010); 

    }
    compRawSpectra(); 
    plotRawSpectra();

    Double_t dndy = 0;
    Double_t deltaPar = 100000.;
    for(int i = 0; i < mKCentBin; i++){
        deltaPar = 100000.;
	while(deltaPar > 0.001){
	    Double_t original_dndy = dndy;      
	    analyzeEff(); // Update the mEff and mEffError
	    compCorrSpectra(); // Compute and Fit and get the dndy
            compDndy();
	    dndy = mLevyPar[i][0]; // use the new efficiency data to update the fitting results
	    deltaPar = dndy - original_dndy; 
	}
    }

    plotEff();
    plotCorrSpectra();
    Double_t realdndy0 = getDndy(0); 
    Double_t realdndy1 = getDndy(1); 
    Double_t realdndy0err = getDndyError(0); 
    Double_t realdndy1err = getDndyError(1); 

    std::cout << "real dndy is " << realdndy0 << " and " << realdndy1 << std::endl;
    std::cout << "real dndyerr is " << realdndy0err << " and " << realdndy1err << std::endl;

    compYields();
    compare11GeV();
}
