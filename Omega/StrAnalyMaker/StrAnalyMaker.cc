#include "StrAnalyMaker.hh"
#include <TLegend.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLine.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <iostream>
#include <string>
#include <fstream>

ClassImp(StrAnalyMaker)
StrAnalyMaker::StrAnalyMaker(std::string par_type):mParticleType(par_type), pdgmass_xi(1.67245), mKCentBin(2), mKPtBin(6){
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
    mDptSpectra[5] = 0.6;

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

    std::cout << "StrAnalyMaker Constructor v0.01 2015-08-30 " << std::endl;
}

StrAnalyMaker::~StrAnalyMaker(){}

void StrAnalyMaker::Init(std::string overview_filename, std::string dat_filename, std::string rotbg_filename, std::string fpeff_filename, std::string expeff_filename){
    std::cout << "!!! InitializationII for AuAu14.5GeV " << mParticleType << " Analysis" << std::endl;
    // Initialize TFile pointers 
    mOverviewFile = new TFile(overview_filename.c_str(), "read");
    mDatFile = new TFile(dat_filename.c_str(), "read");
    mRotBgFile = new TFile(rotbg_filename.c_str(), "read");
    mFpEffFile = new TFile(fpeff_filename.c_str(), "read");
    mExpEffFile = new TFile(expeff_filename.c_str(), "read");

    // Get initial efficiency
    effInit();

    // Initialize Levy Function
    levyInit(); 

    // Get the event number of weighted and unweighted 
    nEventsInit();

    // Initialize branching ratio
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
	    mExpEff[i][j] = 1.0;
	    mEff[i][j] = 1.0;
	}
    }
}

void StrAnalyMaker::levyInit(){
    mLevy = new TF1("levy", "[0]*pow(1+(sqrt(x*x+1.67245*1.67245)-1.67245)/([1]*[2]),-[1])*([1]-1)*([1]-2)/(2*3.14159265*[1]*[2]*([1]*[2]+1.67245*([1]-2)))", 0., 8.);
    mLevyPt = new TF1("levyPt", "x*[0]*pow(1+(sqrt(x*x+1.67245*1.67245)-1.67245)/([1]*[2]),-[1])*([1]-1)*([1]-2)/(2*3.14159265*[1]*[2]*([1]*[2]+1.67245*([1]-2)))", 0., 8.);
    mLevyPt2 = new TF1("levyPt2", "x*x*[0]*pow(1+(sqrt(x*x+1.67245*1.67245)-1.67245)/([1]*[2]),-[1])*([1]-1)*([1]-2)/(2*3.14159265*[1]*[2]*([1]*[2]+1.67245*([1]-2)))", 0., 8.);

    mLevyPar[0][0] = 0.1;   
    mLevyPar[1][0] = 0.1;

    mLevyPar[0][1] = 2000;
    mLevyPar[1][1] = 2000;

    mLevyPar[0][2] = 0.3;
    mLevyPar[1][2] = 0.3;
}

void StrAnalyMaker::nEventsInit(){
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
}

void StrAnalyMaker::rotBgAnalysisInit(){
    std::cout << "!!! Analyze Rotational Background Initialization" << std::endl;
    // Initialize normalization range
    mRotNormLeftLowB = pdgmass_xi - 0.05;
    mRotNormLeftHighB = pdgmass_xi - 0.015;

    mRotNormRightLowB = pdgmass_xi + 0.015;
    mRotNormRightHighB = pdgmass_xi + 0.05;
}

Double_t StrAnalyMaker::compRotNormFactor(Int_t centbin, Int_t ptbin,  TH1F* hdat, TH1F* hrot){
    std::cout << "!!! Compute Rot Norm Factor!" << std::endl;

    Int_t ratio_l1 = hrot->FindBin(mRotNormLeftLowB);
    Int_t ratio_l2 = hrot->FindBin(mRotNormRightLowB);
    Int_t ratio_u1 = hrot->FindBin(mRotNormLeftHighB);
    Int_t ratio_u2 = hrot->FindBin(mRotNormRightHighB);
 
    //mRotScale_ratio[centbin][ptbin] = (hrot->Integral(ratio_l1, ratio_u1) + hrot->Integral(ratio_l2, ratio_u2)) / (hdat->Integral(ratio_l1, ratio_u1) + hdat->Integral(ratio_l2, ratio_u2));    
    mRotScale_ratio[centbin][ptbin] = 1.;
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
    sprintf(plotname, "../omg_plots/rot_%spt%dcent%d.pdf", mParticleType.c_str(), ptbin+1, centbin);
    gPad->SaveAs(plotname);
}

void StrAnalyMaker::compRawSigCounts(Int_t centbin, Int_t ptbin, TH1F* hdat, TH1F* hrot, Double_t scale){
    //TH1F* hrot_copy = (TH1F*)hrot->Clone();
    Int_t sigRangeLeftBin = hdat->FindBin(mSigRangeLeft);
    Int_t sigRangeRightBin = hdat->FindBin(mSigRangeRight);
    Int_t datcounts = hdat->Integral(sigRangeLeftBin, sigRangeRightBin);
    Int_t rotcounts = hrot->Integral(sigRangeLeftBin, sigRangeRightBin);
    mRawSigCounts[centbin][ptbin] = datcounts - rotcounts/scale;
    mRawSigCountsError[centbin][ptbin] = sqrt(datcounts*datcounts + (rotcounts/scale)*(rotcounts/scale));
    std::cout << "mRawSigCounts for cent " << centbin << "pt" << ptbin << "is " << mRawSigCounts[centbin][ptbin] << std::endl;
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

	Int_t noBins = heff_fp->GetSize() - 2; 
        for(int j = 0; j < noBins; j++){
            Double_t effPt = mHFpEffFine[i]->GetBinCenter(j+1);
            Double_t fpEff = mHFpEffFine[i]->GetBinContent(j+1);
            Double_t expEff = mHExpEffFine[i]->GetBinContent(j+1);
            Double_t weight = getSpectraWeight(i, effPt, fpEff);

            pFpEff[i]->Fill(effPt, fpEff, weight);
            pExpEff[i]->Fill(effPt, expEff, weight);
	}
    }  
}

Double_t StrAnalyMaker::getSpectraWeight(Int_t centbin, Double_t pt, Double_t eff){
    Double_t wgt;
    mLevyPt->SetParameters(mLevyPar[centbin]);  
    wgt = mLevyPt->Eval(pt);
    return wgt; 
}

void StrAnalyMaker::compCorrSpectra(){
    //Compute Data points
    for(int i = 0;  i < mKCentBin; i++){
        for(int j = 0; j < mKPtBin; j++){
            mYCorrSpectra[i][j] = mYRawSpectra[i][j] / mEff[i][j]; 
            mYCorrSpectraError[i][j] = mYCorrSpectra[i][j] * sqrt(mYRawSpectraError[i][j]*mYRawSpectraError[i][j] + mEffError[i][j]*mEffError[i][j]); 
            std::cout << "Calculation========YCorrSpectraCent" << i << "Pt" << j<< "= " << mYCorrSpectra[i][j] << "with error = " << mYCorrSpectraError[i][j] << "============" << std::endl; 
	}
    } 

    //Fitting Levy Function, obtain the X position and output the fitting par
    for(int i = 0; i < mKCentBin; i++){
        TGraphErrors* g = new TGraphErrors(mKPtBin, mXCorrSpectra[i], mYCorrSpectra[i], 0, mYCorrSpectraError[i]);
        g->Fit(mLevy, "REM0");
        mLevyPar[i] = mLevy->GetParameters();
        compDndy();
        mLevyPt->SetParameters(mLevyPar[i]);
        mLevyPt2->SetParameters(mLevyPar[i]);

	for(int j = 0; j < mKPtBin; j++){
	    mXCorrSpectra[i][j] = mLevyPt2->Integral(mPtBd[i], mPtBd[i+1]) / mLevyPt->Integral(mPtBd[i], mPtBd[i+1]);
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
    GRawSpectra_1060->GetXaxis()->SetLimits(0.50, 3.60);
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
    leg->AddEntry(GRawSpectra_1060, "10-60%", "p");
    leg->AddEntry(GRawSpectra_010, "0-10%", "p");
    leg->Draw("sames");

    char plotname[50];
    sprintf(plotname, "../%s_plots/%s_rawspectra.pdf", mParticleType.c_str(), mParticleType);
    rawspectra_can->SaveAs(plotname);  
}

void StrAnalyMaker::compDndy(){
    // Initialize measured and unmeasured pt range
    Double_t leftlow = 0.;
    Double_t lefthigh = 0.7;
    Double_t rightlow = 3.6;
    Double_t righthigh = 8.0;
    for(int i = 0; i < mKCentBin; i++){
	mLevy->SetParameters(mLevyPar[i]); 
        for(int j = 0; j < mKPtBin; j++){
	    mDndy[i] += mYCorrSpectra[i][j]*2*3.1415926*mDptSpectra[j]*mNEventsWeighted[i]*mXCorrSpectra[i][j];
	}
	mDndy[i] += mLevy->Integral(leftlow, lefthigh) + mLevy->Integral(rightlow, righthigh);
        //mDndyError[i] =//TODO:
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

    for(int i = 0; i < mKCentBin; i++){
	TGraphErrors* gerr = new TGraphErrors(mKPtBin, mXCorrSpectra[i], mYCorrSpectra[i], 0, mYCorrSPectraError[i]); 
        gerr->SetMarkerSize(1.0);
        gerr->SetMaerkerStyle(20);
        gerr->SetMarkerColor(i+1);
        gerr->SetMaximum(10E-2);
        gerr->SetMinimum(10E-8);
        gerr->GetXaxis()->SetLimits(0.5, 3.60);

        if(i == 0){
             char title[50];
             
             gerr->SetTitle();
	}
    } 
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
       
	hrot_010->Rebin(2);
        hdat_010->Rebin(2);
        hrot_1060->Rebin(2);
        hdat_1060->Rebin(2);

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

    analyzeEff(); // Set the mEff and mEffError

    Double_t dndy = compCorrSpectra();
    Double_t deltaPar = 100000.;
    for(int i = 0; i < mKCentBin; i++){
        deltaPar = 100000.;
	while(deltaPar < 0.001){
	    Double_t original_dndy = dndy;      
	    analyzeEff(); // Update the efficiency data so
	    compCorrSpectra(); // Compute and Fit and get the dndy
	    dndy = mLevyPar[i][0]; // use the new efficiency data to update the fitting results
	    deltaPar = dndy - original_dndy; 
	}
    }

    plotCorrSpectra();
}
