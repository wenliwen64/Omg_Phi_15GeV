#include "StPhiDownMaker.hh"
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

ClassImp(StPhiDownMaker)
StPhiDownMaker::StPhiDownMaker(std::string par_type):mParticleType(par_type), pdgmass_phi(1.01945), mKCentBin(9), mKPtBin(11){
    mCentString[0] = "70%-80%";
    mCentString[1] = "60%-70%";
    mCentString[2] = "50%-60%";
    mCentString[3] = "40%-50%";
    mCentString[4] = "30%-40%";
    mCentString[5] = "20%-30%";
    mCentString[6] = "10%-20%";
    mCentString[7] = "5%-10%";
    mCentString[8] = "0%-5%";

    mPtBd[0] = 0.4;
    mPtBd[1] = 0.5;
    mPtBd[2] = 0.6;
    mPtBd[3] = 0.7;
    mPtBd[4] = 0.8;
    mPtBd[5] = 1.0;
    mPtBd[6] = 1.3;
    mPtBd[7] = 1.7;
    mPtBd[8] = 2.0;
    mPtBd[9] = 2.5;
    mPtBd[10] = 3.5;
    mPtBd[11] = 5.0;

    mXRawSpectra[0] = 0.45;
    mXRawSpectra[1] = 0.55;
    mXRawSpectra[2] = 0.65;
    mXRawSpectra[3] = 0.75;
    mXRawSpectra[4] = 0.90;
    mXRawSpectra[5] = 1.15;
    mXRawSpectra[6] = 1.5;
    mXRawSpectra[7] = 1.85;
    mXRawSpectra[8] = 2.25;
    mXRawSpectra[9] = 3.0;
    mXRawSpectra[10] = 4.25;

    //mXRawSpectraError[0] = 0;
    //mXRawSpectraError[1] = 0;
    //mXRawSpectraError[2] = 0;
    //mXRawSpectraError[3] = 0;
    //mXRawSpectraError[4] = 0;
    //mXRawSpectraError[5] = 0;

    mDptSpectra[0] = 0.1;
    mDptSpectra[1] = 0.1;
    mDptSpectra[2] = 0.1;
    mDptSpectra[3] = 0.1;
    mDptSpectra[4] = 0.2;
    mDptSpectra[5] = 0.3;
    mDptSpectra[6] = 0.4;
    mDptSpectra[7] = 0.3;
    mDptSpectra[8] = 0.5;
    mDptSpectra[9] = 1.0;
    mDptSpectra[10] = 1.5;

    for(int i = 0; i < mKCentBin; i++){
        for(int j = 0; j < mKPtBin; j++)
	    mXCorrSpectra[i][j] = mXRawSpectra[j];
    }

    std::cout << "StPhiDownMaker Constructor v0.01 2015-09-09 " << std::endl;
}

StPhiDownMaker::~StPhiDownMaker(){}

void StPhiDownMaker::Init(std::string dat_filename){
    std::cout << "!!! InitializationII for AuAu14.5GeV " << mParticleType << " Analysis" << std::endl;
    // Initialize TFile pointers 
    //mOverviewFile = new TFile(overview_filename.c_str(), "read");
    mDatFile = new TFile(dat_filename.c_str(), "read");
    //mRotBgFile = new TFile(rotbg_filename.c_str(), "read");
    //mFpEffFile = new TFile(fpeff_filename.c_str(), "read");
    //mExpEffFile = new TFile(expeff_filename.c_str(), "read");

    // Get initial efficiency
    effInit();

    // Initialize Levy Function
    levyInit(); 

    // Get the event number of weighted and unweighted 
    nEventsInit();

    // Initialize branching ratio, phi->K-K+ 
    mBr = 0.489;

    // Initialize signal counting range
    mSigRangeLeft = 1.66;
    mSigRangeRight = 1.685;

    //Initialize rotational background parameters
    rotBgAnalysisInit();
}

void StPhiDownMaker::effInit(){
    std::ifstream infile_flat("weight_phi_fp_15GeV.txt");
    std::ifstream infile_exp("weight_phi_exp_15GeV.txt");
    Int_t centbin = 0;
    Int_t ptbin = 0;
    Float_t dummy = 0;
    Float_t eff = 0;
    Float_t efferr = 0;
    while(infile_flat >> centbin){
	infile_flat >> ptbin >> dummy >> dummy >> dummy >> dummy >> eff >> efferr;         
	mFpEff[centbin][ptbin] = eff;
	mFpEffErr[centbin][ptbin] = efferr;
    }
    while(infile_exp >> centbin){
	infile_exp >> ptbin >> dummy >> dummy >> dummy >> dummy >> eff >> efferr;
	mExpEff[centbin][ptbin] = eff;
	mExpEffErr[centbin][ptbin] = efferr;
    }
}

void StPhiDownMaker::bwFuncInit(){
       total_func = new TF1("total_func", "[0] + [1] * x + 1/(2 * 3.1415926) * [2] * [3] / ((x - [4]) * (x - [4]) + [3] * [3] / 4)", 0.99, 1.05);
       sig_func = new TF1("BW_func", "1/(2 * 3.1415926) * [0] * [1] / ((x - [2]) * (x - [2]) + [1] * [1] / 4)",  0.99, 1.05);
       bg_func = new TF1("bg_func", "[0]+[1]*x", 0.99, 1.05);

       total_func -> SetParName(0, "p0");
       total_func -> SetParName(1, "p1");
       total_func -> SetParName(2, "BW Area");
       total_func -> SetParName(3, "#Gamma");
       total_func -> SetParName(4, "M_{0}");

}

void StPhiDownMaker::levyInit(){
    mLevy = new TF1("levy", "[0]*pow(1+(sqrt(x*x+1.67245*1.67245)-1.67245)/([1]*[2]),-[1])*([1]-1)*([1]-2)/(2*3.14159265*[1]*[2]*([1]*[2]+1.67245*([1]-2)))", 0., 5.);
    mLevyPt = new TF1("levyPt", "x*[0]*pow(1+(sqrt(x*x+1.67245*1.67245)-1.67245)/([1]*[2]),-[1])*([1]-1)*([1]-2)/(2*3.14159265*[1]*[2]*([1]*[2]+1.67245*([1]-2)))", 0., 8.);
    mLevyPt2 = new TF1("levyPt2", "x*x*[0]*pow(1+(sqrt(x*x+1.67245*1.67245)-1.67245)/([1]*[2]),-[1])*([1]-1)*([1]-2)/(2*3.14159265*[1]*[2]*([1]*[2]+1.67245*([1]-2)))", 0., 8.);

    mLevy->SetParName(0, "dN/dy");
    mLevy->SetParName(1, "a");
    mLevy->SetParName(2, "T");
    mLevy->SetLineColor(4);
    mLevy->SetLineStyle(2);
    mLevyPar[0][0] = 0.03;   
    mLevyPar[1][0] = 0.08;

    mLevyPar[0][1] = 2.6e+07;
    mLevyPar[1][1] = -2.6e+06;

    mLevyPar[0][2] = 0.22;
    mLevyPar[1][2] = 0.28;
}

//void StPhiDownMaker::expInit(){
    //mExp = new TF1("exp", );
//}

void StPhiDownMaker::nEventsInit(){
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
   TH1F* h_centbin9_unweighted = (TH1F*)mDatFile->Get("h_centbin9_after1");
   TH1F* h_centbin9_weighted = (TH1F*)mDatFile->Get("h_centbin9_after1");

    for(int i = 0; i < mKCentBin; i++){
        mNEventsWeighted[i] = h_centbin9_weighted->GetBinContent(i+2);  //FIXME
    }
    //mNEventsUnweighted[0] = h_centbin9_unweighted->GetBinContent(8) + h_centbin9_unweighted->GetBinContent(7) + h_centbin9_unweighted->GetBinContent(6) + h_centbin9_unweighted->GetBinContent(5) + h_centbin9_unweighted->GetBinContent(4);
    //mNEventsWeighted[0] = h_centbin9_weighted->GetBinContent(8) + h_centbin9_weighted->GetBinContent(7) + h_centbin9_weighted->GetBinContent(6) + h_centbin9_weighted->GetBinContent(5) + h_centbin9_weighted->GetBinContent(4);
    //mNEventsUnweighted[1] = h_centbin9_unweighted->GetBinContent(10) + h_centbin9_unweighted->GetBinContent(9);
    //mNEventsWeighted[1] = h_centbin9_weighted->GetBinContent(10) + h_centbin9_weighted->GetBinContent(9);

    for(int i = 0; i < mKCentBin; i++){
	std::cout << mNEventsUnweighted[i] << "cent" << i << " nevents unweighted!" << std::endl;
	std::cout << mNEventsWeighted[i] << "cent" << i << " nevents weighted!" << std::endl;
    }

}

void StPhiDownMaker::mixBgAnalysisInit(){
    std::cout << "!!! Analyze Rotational Background Initialization" << std::endl;
    // Initialize normalization range
    mRotNormLeftLowB = //1.625;//pdgmass_xi - 0.05;
    mRotNormLeftHighB = //1.655;//pdgmass_xi - 0.015;

    mRotNormRightLowB = 1.69;//pdgmass_xi + 0.015;
    mRotNormRightHighB = 1.72;//pdgmass_xi + 0.05;
}

Double_t StPhiDownMaker::compMixNormFactor(Int_t centbin, Int_t ptbin,  TH1F* hdat, TH1F* hmix){
    std::cout << "!!! Compute Mixing-Event Norm Factor!" << std::endl;

    Int_t ratio_l1 = hrot->FindBin(mMixNormLeftLowB);
    Int_t ratio_l2 = hrot->FindBin(mMixNormRightLowB);
    Int_t ratio_u1 = hrot->FindBin(mMixNormLeftHighB);
    Int_t ratio_u2 = hrot->FindBin(mMixNormRightHighB);
 
    mMixScale_ratio[centbin][ptbin] = (hrot->Integral(ratio_l1, ratio_u1) + hrot->Integral(ratio_l2, ratio_u2)) / (hdat->Integral(ratio_l1, ratio_u1) + hdat->Integral(ratio_l2, ratio_u2));    
    mMixScale_ratio[centbin][ptbin] = 1.; // In low energies, use 1 as the normalization factor
    std::cout << "------Norm Factor For cent" << centbin << "pt" << ptbin << "is " << mMixScale_ratio[centbin][ptbin] << std::endl;
    return mMixScale_ratio[centbin][ptbin];
}

void StPhiDownMaker::plotMixInvMassWithData(Int_t centbin, Int_t ptbin, TH1F* hdat, TH1F* hmix, Double_t scale){
    
    std::cout << "!!! Plot Inv Mass" << std::endl;
    //hdat->SetMarkerStyle(8);
    hdat->Draw("PE");
    //hdat->SetLineColor(4);
    //hdat->SetFillColorAlpha(4, 0.35);
    hdat->GetXaxis()->SetTitle("InvMass(GeV)");
    hdat->GetYaxis()->SetTitle("Counts");

    TH1F* hmix_copy = (TH1F*)hmix->Clone();
    hmix_copy->Sumw2();
    hmix_copy->Scale(1./scale);
    hmix_copy->SetLineColor(2);
    //hrot_copy->SetFillColorAlpha(2, 0.35);
    //hrot_copy->SetLineWidth(0.5);
    hmix_copy->SetFillColor(2);
    hmix_copy->SetFillStyle(3354);
    gPad->SetTicks(1, 1);
    hrot_copy->Draw("PEsames");

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
    sprintf(plotname, "../%s_plots/mix_%spt%dcent%d.pdf", mParticleType.c_str(), mParticleType.c_str(), ptbin+1, centbin);
    gPad->SaveAs(plotname);
}

void StPhiDownMaker::compRawSigCounts(Int_t centbin, Int_t ptbin, TH1F* hdat, TH1F* hmix, Double_t scale){
    //TH1F* hrot_copy = (TH1F*)hrot->Clone();
    Int_t sigRangeLeftBin = hdat->FindBin(mSigRangeLeft);
    Int_t sigRangeRightBin = hdat->FindBin(mSigRangeRight);
    Int_t datcounts = hdat->Integral(sigRangeLeftBin, sigRangeRightBin);
    Int_t mixcounts = hmix->Integral(sigRangeLeftBin, sigRangeRightBin);
    mRawSigCounts[centbin][ptbin] = datcounts - mixcounts/scale;
    mRawSigCountsError[centbin][ptbin] = sqrt(datcounts + (rotcounts/scale));
    std::cout << "!!! mRawSigCounts for cent " << centbin << "pt" << ptbin << " is " << mRawSigCounts[centbin][ptbin] << " error = " << mRawSigCountsError[centbin][ptbin] << std::endl;
}

void StPhiDownMaker::compRawSpectra(){
    double PI = 3.1415926; 
    for(int i = 0; i < mKCentBin; i++){
        for(int j = 0; j < mKPtBin; j++){
            mYRawSpectra[i][j] = 1/(2*PI) * mRawSigCounts[i][j] / mXRawSpectra[j] / mDptSpectra[j] / mNEventsWeighted[i] / mBr;
            mYRawSpectraError[i][j] = 1/(2*PI) * mRawSigCountsError[i][j] / mXRawSpectra[j] / mDptSpectra[j] / mNEventsWeighted[i] / mBr;
            std::cout << "mYRawSpectra for cent" << i << "pt" << j <<" is " << mYRawSpectra[i][j] << "mYRawSpectraError is " << mYRawSpectraError[i][j] << std::endl;
	}
    }
}

void StPhiDownMaker::plotRawSpectra(){ 
    TCanvas* rawspectra_can = new TCanvas("rawspectra_can", "rawspectra_can");
    rawspectra_can->SetLogy();
    rawspectra_can->SetTicks(1, 1);

    TGraphErrors* GRawSpectra_1060 = new TGraphErrors(6, mXRawSpectra, mYRawSpectra[0], mXRawSpectraError, mYRawSpectraError[0]); 
    GRawSpectra_1060->SetMarkerSize(1.0);
    GRawSpectra_1060->SetMarkerStyle(34);
    GRawSpectra_1060->SetMarkerColor(2);
    GRawSpectra_1060->SetMaximum(10E-1);
    GRawSpectra_1060->SetMinimum(10E-9);
    GRawSpectra_1060->GetXaxis()->SetLimits(0.0, 3.60);
    GRawSpectra_1060->SetTitle("#Omega^{-} Spectra, Au+Au 14.5GeV");
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


void StPhiDownMaker::analyzeEff(){
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
            if(k < 2){
                mEff[i][k] = pExpEff[i]->GetBinContent(k+1);
                mEffError[i][k] = pExpEff[i]->GetBinError(k+1);
            }
            else{
                mEff[i][k] = pFpEff[i]->GetBinContent(k+1);
                mEffError[i][k] = pFpEff[i]->GetBinError(k+1);
	    }
            //std::cout << "!!! mEff for cent" << i << "pt" << k << " is " << mEff[i][k] << std::endl;
	}
    }
}

std::string StPhiDownMaker::getCentString(Int_t i){
    return mCentString[i];
}

void StPhiDownMaker::plotEff(){
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
            geff->SetTitle("#Omega^{-} Efficiency, Au+Au 14.5GeV"); 
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

Double_t StPhiDownMaker::getSpectraWeight(Int_t centbin, Double_t pt){
    Double_t wgt;
    //std::cout << "!!! mLevyPar = " << mLevyPar[centbin][0] << " " << mLevyPar[centbin][1] << " " << mLevyPar[centbin][2] << std::endl;
    mLevyPt->SetParameters(mLevyPar[centbin]);
    wgt = mLevyPt->Eval(pt);
    return wgt;
}

void StPhiDownMaker::compCorrSpectra(){
    //Compute Data points
    for(int i = 0;  i < mKCentBin; i++){
        for(int j = 0; j < mKPtBin; j++){
            mYCorrSpectra[i][j] = mYRawSpectra[i][j] / mEff[i][j]; 
            double relative_rawyerror = mYRawSpectraError[i][j] / mYRawSpectra[i][j];
            double relative_efferror = mEffError[i][j] / mEff[i][j];
            mYCorrSpectraError[i][j] = mYCorrSpectra[i][j] * sqrt(relative_rawyerror*relative_rawyerror); 
            //mYCorrSpectraError[i][j] = mYCorrSpectra[i][j] * sqrt(relative_rawyerror*relative_rawyerror+relative_efferror*relative_efferror); 
            std::cout << "Calculation========YCorrSpectraCent" << i << "Pt" << j<< "= " << mYCorrSpectra[i][j] << "with error = " << mYCorrSpectraError[i][j] << " with efficiency being " << "============" << std::endl; 
	}
    } 

    //Fitting Levy Function, obtain the X position and output the fitting par
    for(int i = 0; i < mKCentBin; i++){
        TGraphErrors* g = new TGraphErrors(mKPtBin, mXCorrSpectra[i], mYCorrSpectra[i], 0, mYCorrSpectraError[i]);
        mLevy->SetParameters(mLevyPar[i]);
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

void StPhiDownMaker::compDndy(){
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
        mDndyFit[i] = 2*3.1415926*mLevyPt->Integral(leftlow, righthigh);
	mDndy[i] += 2*3.1415926*(mLevyPt->Integral(leftlow, lefthigh) + mLevyPt->Integral(rightlow, righthigh));
        tmp_dndy_err1 = mLevyParError[i][0]/mLevyPar[i][0]*2*3.1415926*(mLevyPt->Integral(leftlow, lefthigh));
        tmp_dndy_err2 = mLevyParError[i][0]/mLevyPar[i][0]*2*3.1415926*(mLevyPt->Integral(rightlow, righthigh));
        mDndyError[i] = sqrt(tmp_dndy_err02 + tmp_dndy_err1*tmp_dndy_err1 + tmp_dndy_err2*tmp_dndy_err2);//TODO:
        std::cout << "err = " << tmp_dndy_err02 << " " << tmp_dndy_err1*tmp_dndy_err1 << " " << tmp_dndy_err2*tmp_dndy_err2 << std::endl;
    }
}

Double_t StPhiDownMaker::getDndy(Int_t centbin){
    return mDndy[centbin];
}
 
Double_t StPhiDownMaker::getDndyError(Int_t centbin){
    return mDndyError[centbin];
}

void StPhiDownMaker::plotCorrSpectra(){
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
	    std::string title = "#" + mParticleType + "^{-} Spectra, Au+Au 14.5GeV"; 
	    gerr->SetTitle(title.c_str());
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
    sprintf(plotname, "../%s_plots/finalCorrSpectra.pdf", mParticleType.c_str());
    canCorrSpectra->SaveAs(plotname);
    //canCorrSpectra->SaveAs("../omg_plots/finalCorrSpectra.gif");
    //canCorrSpectra->SaveAs("../omg_plots/finalCorrSpectra.eps");
    //canCorrSpectra->SaveAs("../omg_plots/finalCorrSpectra.jpg");
}

void StPhiDownMaker::Analyze(){
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
    //std::cout << "fit dndy is " << mDndyFit[0] << " and " << mDndyFit[1] << std::endl;
}
