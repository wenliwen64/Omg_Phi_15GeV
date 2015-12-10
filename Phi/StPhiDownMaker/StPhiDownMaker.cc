#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */
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
    std::cout << GREEN << "<----- StPhiDownMaker Constructor v0.01 2015-12-09 ----->" << RESET << std::endl;
    mCentString[0] = "70%-80%";
    mCentString[1] = "60%-70%";
    mCentString[2] = "50%-60%";
    mCentString[3] = "40%-50%";
    mCentString[4] = "30%-40%";
    mCentString[5] = "20%-30%";
    mCentString[6] = "10%-20%";
    mCentString[7] = "5%-10%";
    mCentString[8] = "0%-5%";

    mCentStringBES[0] = "60%-80%";
    mCentStringBES[1] = "40%-60%";
    mCentStringBES[2] = "30%-40%";
    mCentStringBES[3] = "20%-30%";
    mCentStringBES[4] = "10%-20%";
    mCentStringBES[5] = "0%-10%";

    mNColl[0] = 12.1316;
    mNColl[1] = 25.2342;
    mNColl[2] = 49.8735;
    mNColl[3] = 94.3089;
    mNColl[4] = 167.971;
    mNColl[5] = 282.678;
    mNColl[6] = 454.215;
    mNColl[7] = 633.514;
    mNColl[8] = 787.915;

    mNCollError[0] = 5.08722;
    mNCollError[1] = 8.78222;
    mNCollError[2] = 12.3042;
    mNCollError[3] = 17.6944;
    mNCollError[4] = 22.1882;
    mNCollError[5] = 24.2263;
    mNCollError[6] = 23.5591;
    mNCollError[7] = 20.4984;
    mNCollError[8] = 29.8901;

    mNCollBES[0] = 18.3058;
    mNCollBES[1] = 71.739;
    mNCollBES[2] = 167.971;
    mNCollBES[3] = 282.678;
    mNCollBES[4] = 454.215;
    mNCollBES[5] = 711.486;

    mNCollErrorBES[0] = 6.34844;
    mNCollErrorBES[1] = 15.1698;
    mNCollErrorBES[2] = 22.1882;
    mNCollErrorBES[3] = 24.2263;
    mNCollErrorBES[4] = 23.5591;
    mNCollErrorBES[5] = 27.3585;


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
    mXRawSpectra[6] = 1.50;
    mXRawSpectra[7] = 1.85;
    mXRawSpectra[8] = 2.25;
    mXRawSpectra[9] = 3.0;
    mXRawSpectra[10] = 4.25;

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

    mXRawSpectraOmg[0] = 0.95;
    mXRawSpectraOmg[1] = 1.40;
    mXRawSpectraOmg[2] = 1.80;
    mXRawSpectraOmg[3] = 2.20;
    mXRawSpectraOmg[4] = 2.60;
    mXRawSpectraOmg[5] = 3.20;

    mPtBdOmg[0] = 0.7;
    mPtBdOmg[1] = 1.2;
    mPtBdOmg[2] = 1.6;
    mPtBdOmg[3] = 2.0;
    mPtBdOmg[4] = 2.4;
    mPtBdOmg[5] = 2.8;
    mPtBdOmg[6] = 3.6;

    mOmgYields010[0] = 0.0305707;
    mOmgYields010[1] = 0.0216205;
    mOmgYields010[2] = 0.0114018;
    mOmgYields010[3] = 0.00390573;
    mOmgYields010[4] = 0.00120093;
    mOmgYields010[5] = 0.000509963;

    mAntiOmgYields010[0] = 0.0130111;
    mAntiOmgYields010[1] = 0.008431;
    mAntiOmgYields010[2] = 0.00412566;
    mAntiOmgYields010[3] = 0.00228306;
    mAntiOmgYields010[4] = 0.000983013;
    mAntiOmgYields010[5] = 0.000379698;

    for(int i = 0; i < mKCentBin; i++){
        for(int j = 0; j < mKPtBin; j++)
	    mXCorrSpectra[i][j] = mXRawSpectra[j];
    }
// BES centrality binning initialization
    for(int i = 0; i < 6; i++){
        for(int j = 0; j < mKPtBin; j++)
	    mXCorrSpectraBES[i][j] = mXRawSpectra[j];
    }
}

StPhiDownMaker::~StPhiDownMaker(){}

void StPhiDownMaker::Init(std::string dat_filename){
    std::cout << CYAN << ">> Initializing " << mParticleType << " Analysis Maker..." << RESET << std::endl;
    // Initialize TFile pointers 
    mDatFile = new TFile(dat_filename.c_str(), "read");

    // Initialize branching ratio, phi->K-K+ 
    mBr = 0.491;

    // Initialize signal counting range
    mSigRangeLeft = 0.9;
    mSigRangeRight = 1.05;

    // Initialize inv-mass fitting function BW
    bwFuncInit();

    // Get initial efficiency
    effInit();

    // Initialize Levy Function
    levyInit(); 

    // Get the event number of weighted and unweighted 
    nEventsInit();

    //Initialize rotational background parameters
    mixBgAnalysisInit();
}

void StPhiDownMaker::bwFuncInit(){
    std::cout << YELLOW << ".... BW function initialization..." << RESET << std::endl;
    mTotal = new TF1("total_func", "[0] + [1] * x + 1/(2 * 3.1415926) * [2] * [3] / ((x - [4]) * (x - [4]) + [3] * [3] / 4)", 0.99, 1.05);
    mBW = new TF1("BW_func", "1/(2 * 3.1415926) * [0] * [1] / ((x - [2]) * (x - [2]) + [1] * [1] / 4)",  0.99, 1.05);
    mPolyBg = new TF1("bg_func", "[0]+[1]*x", 0.99, 1.05);

    mTotal->SetParName(0, "p0");
    mTotal->SetParName(1, "p1");
    mTotal->SetParName(2, "BW Area");
    mTotal->SetParName(3, "#Gamma");
    mTotal->SetParName(4, "M_{0}");

    mTotal->SetParameter(3, 0.007);
    mTotal->SetParameter(4, 1.0195);
    mTotal->SetParLimits(3, 0.006, 0.007);
    mTotal->SetParLimits(4, 1.015, 1.023);
}

void StPhiDownMaker::effInit(){
    std::cout << YELLOW << ".... Efficiency initialization..." << std::endl;
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
	mFpEffError[centbin][ptbin] = efferr;
    }
    while(infile_exp >> centbin){
	infile_exp >> ptbin >> dummy >> dummy >> dummy >> dummy >> eff >> efferr;
	mExpEff[centbin][ptbin] = eff;
	mExpEffError[centbin][ptbin] = efferr;
    }

    for(int i = 0; i < mKCentBin; i++){
        for(int j = 0; j < mKPtBin; j++){
            if(j < 3){ //TODO: For first 3 pt bin to use exp dist. eff. 
                mEff[i][j] = mExpEff[i][j];
                mEffError[i][j] = mExpEffError[i][j];
	    }
            else{
                mEff[i][j] = mFpEff[i][j];
                mEffError[i][j] = mFpEffError[i][j];
	    }
	}
    }
}


void StPhiDownMaker::levyInit(){
    std::cout << YELLOW << ".... levy function initialization" << RESET << std::endl;
    mLevy = new TF1("levy", "[0]*pow(1+(sqrt(x*x+1.67245*1.67245)-1.67245)/([1]*[2]),-[1])*([1]-1)*([1]-2)/(2*3.14159265*[1]*[2]*([1]*[2]+1.67245*([1]-2)))", 0., 5.);
    mLevyPt = new TF1("levyPt", "x*[0]*pow(1+(sqrt(x*x+1.67245*1.67245)-1.67245)/([1]*[2]),-[1])*([1]-1)*([1]-2)/(2*3.14159265*[1]*[2]*([1]*[2]+1.67245*([1]-2)))", 0., 8.);
    mLevyPt2 = new TF1("levyPt2", "x*x*[0]*pow(1+(sqrt(x*x+1.67245*1.67245)-1.67245)/([1]*[2]),-[1])*([1]-1)*([1]-2)/(2*3.14159265*[1]*[2]*([1]*[2]+1.67245*([1]-2)))", 0., 8.);

    mLevy->SetParName(0, "dN/dy");
    mLevy->SetParName(1, "a");
    mLevy->SetParName(2, "T");
    mLevy->SetLineColor(4);
    mLevy->SetLineStyle(2);
    for(int i = 0; i < mKCentBin; i++){
	mLevyPar[i][0] = 0.3*pow(10, -mKCentBin+i+1);   
	mLevyPar[i][1] = 2.6e+07;
	mLevyPar[i][2] = 0.25;
    }

    for(int i = 0; i < 6; i++){
	mLevyParBES[i][0] = 1.0*pow(10, -6+i+1);   
	mLevyParBES[i][1] = 2.6e+07;
	mLevyParBES[i][2] = 0.3;
    }
}

void StPhiDownMaker::nEventsInit(){
    //TODO: should get the overview soon
    std::cout << YELLOW << ".... Nevents initialization..." << RESET << std::endl;
    TH1F* h_centbin9_unweighted = (TH1F*)mDatFile->Get("hCentBin9_after");
    TH1F* h_centbin9_weighted = (TH1F*)mDatFile->Get("hCentBin9_after");

    for(int i = 0; i < mKCentBin; i++){
	mNEventsUnweighted[i] = h_centbin9_unweighted->GetBinContent(i+2);
	mNEventsWeighted[i] = h_centbin9_weighted->GetBinContent(i+2);
    }

    //BES centrality bins
    mNEventsWeightedBES[0] = mNEventsWeighted[0] + mNEventsWeighted[1]; // 60-80%
    mNEventsWeightedBES[1] = mNEventsWeighted[2] + mNEventsWeighted[3]; // 40-60%
    mNEventsWeightedBES[2] = mNEventsWeighted[4]; // 30-40%
    mNEventsWeightedBES[3] = mNEventsWeighted[5]; // 20-30%
    mNEventsWeightedBES[4] = mNEventsWeighted[6]; // 10-20%
    mNEventsWeightedBES[5] = mNEventsWeighted[7] + mNEventsWeighted[8]; // 0-10%

    for(int i = 0; i < 6; i++)
	std::cout << BLUE << "........ nEvents in cen" << i << ": " << mNEventsWeightedBES[i] << RESET << std::endl;
}

void StPhiDownMaker::mixBgAnalysisInit(){
    std::cout << YELLOW << ".... Mix-event background initialization..." << RESET << std::endl;
    // Initialize normalization range
    mMixNormLeftLowB = 0; //TODO:1.625;//pdgmass_xi - 0.05;
    mMixNormLeftHighB = 0; //TODO 1.655;//pdgmass_xi - 0.015;

    mMixNormRightLowB = 1.05; //TODO pdgmass_xi + 0.015;
    mMixNormRightHighB = 1.07; //TODO pdgmass_xi + 0.05;
}

Double_t StPhiDownMaker::compMixNormFactor(Int_t centbin, Int_t ptbin,  TH1F* hdat, TH1F* hmix){
    //std::cout << YELLOW << ".... Compute mixing-event norm factor..." << RESET << std::endl;

    Int_t ratio_l1 = hmix->FindBin(mMixNormLeftLowB);
    Int_t ratio_l2 = hmix->FindBin(mMixNormRightLowB);
    Int_t ratio_u1 = hmix->FindBin(mMixNormLeftHighB);
    Int_t ratio_u2 = hmix->FindBin(mMixNormRightHighB);
 
    mMixScale_ratio[centbin][ptbin] = (hmix->Integral(ratio_l1, ratio_u1) + hmix->Integral(ratio_l2, ratio_u2)) / (hdat->Integral(ratio_l1, ratio_u1) + hdat->Integral(ratio_l2, ratio_u2));    

    std::cout << BLUE << "........ Mix-event background normalization factor for cent " << centbin << " ptbin " << ptbin << " is " << mMixScale_ratio[centbin][ptbin] << RESET << std::endl;

    return mMixScale_ratio[centbin][ptbin];
}

void StPhiDownMaker::plotInvMassAfterBgSubtraction(Int_t centbin, Int_t ptbin, TH1F* hdat, TH1F* hmix, Double_t scale){
    //std::cout << YELLOW << ".... Plot inv-mass after bg subtraction..." << RESET << std::endl;
    TCanvas* c = new TCanvas();

    TH1F* hdat_copy = (TH1F*)hdat->Clone();
    TH1F* hmix_copy = (TH1F*)hmix->Clone();
    //hdat_copy->Sumw2();
    hdat_copy->Add(hmix_copy, -1./scale);
    hdat_copy->SetMarkerStyle(24);
    hdat_copy->SetMarkerColor(1);
    hdat_copy->SetLineColor(1);
    gPad->SetTicks(1, 1);
    hdat_copy->Draw("PE");

    hdat_copy->Fit(mTotal, "QREM"); 
    
    Double_t par[5]; 
    Double_t* parerr = NULL;
    mTotal->GetParameters(par);
    parerr = mTotal->GetParErrors();
    mInvMassPar[0][centbin][ptbin] = par[2];
    mInvMassPar[1][centbin][ptbin] = par[3];
    mInvMassPar[2][centbin][ptbin] = par[4];
    mInvMassParError[0][centbin][ptbin] = parerr[2];
    mInvMassParError[1][centbin][ptbin] = parerr[3];
    mInvMassParError[2][centbin][ptbin] = parerr[4];

    mBW->SetParameter(0, mInvMassPar[0][centbin][ptbin]);
    mBW->SetParameter(1, mInvMassPar[1][centbin][ptbin]);
    mBW->SetParameter(2, mInvMassPar[2][centbin][ptbin]);
    mBW->SetLineColor(3);
    mBW->Draw("same");

    mPolyBg->SetParameter(0, par[0]);
    mPolyBg->SetParameter(1, par[1]);
    mPolyBg->SetLineColor(4);
    mPolyBg->Draw("same");

    TLine* lline = new TLine(mInvMassPar[2][centbin][ptbin] - 3 * mInvMassPar[1][centbin][ptbin], 0, mInvMassPar[2][centbin][ptbin] - 3 * mInvMassPar[1][centbin][ptbin], hdat_copy->GetMaximum());
    TLine* uline = new TLine(mInvMassPar[2][centbin][ptbin] + 3 * mInvMassPar[1][centbin][ptbin], 0, mInvMassPar[2][centbin][ptbin] + 3 * mInvMassPar[1][centbin][ptbin], hdat_copy->GetMaximum());
    lline->SetLineColor(6);
    lline->SetLineStyle(10);
    lline->SetLineWidth(2);
    lline->Draw("same");

    uline->SetLineColor(6);
    uline->SetLineStyle(10);
    uline->SetLineWidth(2);
    uline->Draw("same");
    
    char plotname[50];
    sprintf(plotname, "../%s_plots/pure_%scent%dpt%d.pdf", mParticleType.c_str(), mParticleType.c_str(), centbin+1, ptbin+1);
    gPad->SaveAs(plotname);
    delete c;
}

void StPhiDownMaker::plotMixInvMassWithData(Int_t centbin, Int_t ptbin, TH1F* hdat, TH1F* hmix, Double_t scale){
    //std::cout << "........ Plot Inv Mass" << std::endl;
    TH1F* hdat_copy = (TH1F*)hdat->Clone();

    TCanvas* c = new TCanvas();
    hdat_copy->Draw("PE");
    hdat_copy->GetXaxis()->SetTitle("InvMass(GeV)");
    hdat_copy->GetYaxis()->SetTitle("Counts");

    TH1F* hmix_copy = (TH1F*)hmix->Clone();
    //hmix_copy->Sumw2();
    hmix_copy->Scale(1./scale);
    hmix_copy->SetLineColor(2);
    hmix_copy->SetFillColor(2);
    hmix_copy->SetFillStyle(3354);
    hmix_copy->Draw("PEsames");

    TLine* lline = new TLine(mSigRangeLeft, 0, mSigRangeLeft, hdat->GetMaximum());
    TLine* uline = new TLine(mSigRangeRight, 0, mSigRangeRight, hdat->GetMaximum());
    lline->SetLineColor(4);
    lline->SetLineWidth(2);
    lline->SetLineStyle(10);
    uline->SetLineColor(4);
    uline->SetLineWidth(2);
    uline->SetLineStyle(10);
    lline->Draw("sames");
    uline->Draw("sames");

    TLine* l1line = new TLine(mMixNormLeftLowB, 0, mMixNormLeftLowB, hdat->GetMaximum());
    TLine* u1line = new TLine(mMixNormLeftHighB, 0, mMixNormLeftHighB, hdat->GetMaximum());
    l1line->SetLineColor(6);
    l1line->SetLineWidth(2);
    l1line->SetLineStyle(10);
    u1line->SetLineColor(6);
    u1line->SetLineWidth(2);
    u1line->SetLineStyle(10);
    l1line->Draw("sames");
    u1line->Draw("sames");

    TLine* l2line = new TLine(mMixNormRightLowB, 0, mMixNormRightLowB, hdat->GetMaximum());
    TLine* u2line = new TLine(mMixNormRightHighB, 0, mMixNormRightHighB, hdat->GetMaximum());
    l2line->SetLineColor(6);
    l2line->SetLineWidth(2);
    l2line->SetLineStyle(10);
    u2line->SetLineColor(6);
    u2line->SetLineWidth(2);
    u2line->SetLineStyle(10);
    l2line->Draw("sames");
    u2line->Draw("sames"); 
    gPad->SetTicks(1, 1);
    
    char plotname[50];
    sprintf(plotname, "../%s_plots/mix_cent%dpt%d_mix.pdf", mParticleType.c_str(), centbin+1, ptbin+1);
    gPad->SaveAs(plotname);
    delete c;
}

void StPhiDownMaker::compRawSigCounts(Int_t centbin, Int_t ptbin, Double_t bin_width){
    // Using fit function integration to count the signal number
    mBW->SetParameter(0, mInvMassPar[0][centbin][ptbin]);
    mBW->SetParameter(1, mInvMassPar[1][centbin][ptbin]);
    mBW->SetParameter(2, mInvMassPar[2][centbin][ptbin]);
    mRawSigCounts[centbin][ptbin] = mBW->Integral(mInvMassPar[2][centbin][ptbin] - 3 * mInvMassPar[1][centbin][ptbin], mInvMassPar[2][centbin][ptbin] + 3 * mInvMassPar[1][centbin][ptbin]) / bin_width;
    mRawSigCountsError[centbin][ptbin] = mInvMassParError[0][centbin][ptbin] / mInvMassPar[0][centbin][ptbin] * mRawSigCounts[centbin][ptbin];
    std::cout << BLUE << "........ mRawSigCounts for cent " << centbin << " pt " << ptbin << " is " << mRawSigCounts[centbin][ptbin] << " error = " << mRawSigCountsError[centbin][ptbin] << RESET << std::endl;

}

void StPhiDownMaker::plotInvMassQA(){

     std::cout << YELLOW << ".... Plotting inv-mass QA... " << RESET << std::endl;
     TCanvas* can = new TCanvas();
     TLegend* legend = new TLegend(0.55, 0.15, 0.80, 0.45);
     legend->SetBorderSize(0);
     for(int j = 0; j < mKCentBin; j++){
	 Int_t i = mKCentBin - j - 1;
	 TGraphErrors* gerr = new TGraphErrors(mKPtBin - 1, mXRawSpectra, mInvMassPar[2][i], 0, mInvMassParError[2][i]);
	 gerr->SetMarkerStyle(20);
	 gerr->SetMaximum(1.03);
	 gerr->SetMinimum(0.09);
	 gerr->GetXaxis()->SetTitle("pT(GeV/c)");
	 gerr->GetYaxis()->SetTitle("Mass(GeV/c^{2})");
	 gerr->GetYaxis()->SetTitleOffset(1.3);
	 gerr->SetMarkerColor(i+1);
	 gerr->SetTitle("BW Mass");
	 if(j == 0)
	     gerr->Draw("AP");
	 else
	     gerr->Draw("P same");

	 legend->AddEntry(gerr, mCentString[i].c_str(), "p");
     }
     legend->Draw("same");
     gPad->SetTicks(1, 1);
     gPad->SaveAs("../Phi_plots/InvMassQA_mass.pdf");
}

void StPhiDownMaker::compRawSigCountsBES(){
    std::cout << YELLOW << ".... Computing BES centrality bin raw signal counts..." << RESET << std::endl;
    for(int i = 0; i < 6; i++){
	for(int j = 0; j < mKPtBin; j++){
	    if(i == 0){
		mRawSigCountsBES[i][j] = mRawSigCounts[0][j] + mRawSigCounts[1][j];
		mRawSigCountsBESError[i][j] = sqrt(mRawSigCountsError[0][j]*mRawSigCountsError[0][j] + mRawSigCountsError[1][j]*mRawSigCountsError[1][j]);
	    }
	    else if(i == 1){
		mRawSigCountsBES[i][j] = mRawSigCounts[2][j] + mRawSigCounts[3][j];
		mRawSigCountsBESError[i][j] = sqrt(mRawSigCountsError[2][j]*mRawSigCountsError[2][j] + mRawSigCountsError[3][j]*mRawSigCountsError[3][j]);
                //std::cout << i << " " << j << " " << mRawSigCounts[2][j] << " (((())))) " << mRawSigCounts[3] << std::endl;
	    }
	    else if(i >= 2 && i <=4){
		mRawSigCountsBES[i][j] = mRawSigCounts[i+2][j];
		mRawSigCountsBESError[i][j] = mRawSigCountsError[i+2][j];
	    }
	    else{
		mRawSigCountsBES[i][j] = mRawSigCounts[7][j] + mRawSigCounts[8][j];
		mRawSigCountsBESError[i][j] = sqrt(mRawSigCountsError[7][j]*mRawSigCountsError[7][j] + mRawSigCountsError[8][j]*mRawSigCountsError[8][j]);
	    }

            std::cout << BLUE << "........ BESRawCounts: cen" <<i << " pt" << j << " raw counts = " << mRawSigCountsBES[i][j] << RESET << std::endl;
	}
    }
}

void StPhiDownMaker::compRawSpectra(){
    std::cout << YELLOW << ".... Computing raw spectra for 9 cen bins..." << RESET << std::endl;
    double PI = 3.1415926;
    for(int i = 0; i < mKCentBin; i++){
        for(int j = 0; j < mKPtBin; j++){
            mYRawSpectra[i][j] = 1/(2*PI) * mRawSigCounts[i][j] / mXRawSpectra[j] / mDptSpectra[j] / mNEventsWeighted[i] / mBr;
            mYRawSpectraScale[i][j] = pow(10, i-mKCentBin+1) * mYRawSpectra[i][j];
            mYRawSpectraError[i][j] = 1/(2*PI) * mRawSigCountsError[i][j] / mXRawSpectra[j] / mDptSpectra[j] / mNEventsWeighted[i] / mBr;
            mYRawSpectraErrorScale[i][j] = pow(10, i-mKCentBin+1) * mYRawSpectraError[i][j];
            std::cout << BLUE << "........ RawSpectra: cent" << i << " pt" << j <<" y =  " << mYRawSpectra[i][j] << ", yerr = " << mYRawSpectraError[i][j] << RESET << std::endl;
	}
    }
}

void StPhiDownMaker::compRawSpectraBES(){
    std::cout << YELLOW << ".... Computing raw spectra for BES cen bins..." << RESET << std::endl;
    double PI = 3.1415926; 
    for(int i = 0; i < 6; i++){
        for(int j = 0; j < mKPtBin; j++){
            mYRawSpectraBES[i][j] = 1/(2*PI) * mRawSigCountsBES[i][j] / mXRawSpectra[j] / mDptSpectra[j] / mNEventsWeightedBES[i] / mBr;
            mYRawSpectraBESScale[i][j] = pow(10, i-5) * mYRawSpectraBES[i][j];
            mYRawSpectraBESError[i][j] = 1/(2*PI) * mRawSigCountsBESError[i][j] / mXRawSpectra[j] / mDptSpectra[j] / mNEventsWeightedBES[i] / mBr;
            mYRawSpectraBESErrorScale[i][j] = pow(10, i-6+1) * mYRawSpectraBESError[i][j];
            std::cout << BLUE << "........ RawSpectraBES: cent" << i << "pt" << j <<" y = " << mYRawSpectraBES[i][j] << ", yerr = " << mYRawSpectraBESError[i][j] << RESET << std::endl;
	}
    }
}

void StPhiDownMaker::plotRawSpectra(){ 
    std::cout << YELLOW << ".... Plotting raw spectra for 9 cen bins..." << RESET << std::endl;
    TCanvas* rawspectra_can = new TCanvas("rawspectra_can", "rawspectra_can");
    rawspectra_can->SetLogy();
    rawspectra_can->SetTicks(1, 1);

    TGraphErrors* gr = NULL;

    TLegend* leg = new TLegend(0.60, 0.55, 0.85, 0.85);
    leg->SetBorderSize(1);
    
    Int_t i = mKCentBin;
    while(i > 0){
        i--;
	gr = new TGraphErrors(mKPtBin-1, mXRawSpectra, mYRawSpectraScale[i], 0, mYRawSpectraErrorScale[i]); 
	gr->SetMarkerSize(1.0);
	gr->SetMarkerStyle(34);
	gr->SetMarkerColor(1+i);
	if(i == (mKCentBin-1)){
	    gr->SetMaximum(1); 
	    gr->SetMinimum(10e-17);
	    gr->GetXaxis()->SetLimits(0.0, 5.0);
	    gr->GetXaxis()->SetTitle("Pt(GeV/c)");
	    gr->GetYaxis()->SetTitle("#frac{d^{2}N}{2#piNP_{T}dP_{T}dy}(GeV/c)^{-2}");
	    gr->GetYaxis()->SetTitleOffset(1.2);
	    gr->SetTitle("#phi Spectra, Au+Au 14.5GeV");
	    gr->Draw("AP"); 

            leg->AddEntry(gr, mCentString[i].c_str(), "p");
	}
        else{
            gr->Draw("P same");
            char leg_string[50];
            sprintf(leg_string, "%s#times 10^{-%d}", mCentString[i].c_str(), mKCentBin-i-1);
            leg->AddEntry(gr, leg_string, "p");
	}
    }
    leg->Draw("sames");

    std::string plotname = "../" + mParticleType + "_plots/" + mParticleType + "_rawspectra.pdf";
    rawspectra_can->SaveAs(plotname.c_str());  
}

void StPhiDownMaker::plotRawSpectraBES(){ 
    std::cout << YELLOW << ".... Plotting raw spectra for BES cen bins..." << RESET << std::endl;
    TCanvas* rawspectra_can = new TCanvas("rawspectra_can", "rawspectra_can");
    rawspectra_can->SetLogy();
    rawspectra_can->SetTicks(1, 1);

    TGraphErrors* gr = NULL;

    TLegend* leg = new TLegend(0.60, 0.55, 0.85, 0.85);
    leg->SetBorderSize(1);
    
    Int_t i = 6;
    while(i > 0){
        i--;
	gr = new TGraphErrors(mKPtBin-1, mXRawSpectra, mYRawSpectraBESScale[i], 0, mYRawSpectraBESErrorScale[i]); 
	gr->SetMarkerSize(1.0);
	gr->SetMarkerStyle(34);
	gr->SetMarkerColor(1+i);
	if(i == 5){
	    gr->SetMaximum(1); 
	    gr->SetMinimum(10e-17);
	    gr->GetXaxis()->SetLimits(0.0, 5.0);
	    gr->GetXaxis()->SetTitle("Pt(GeV/c)");
	    gr->GetYaxis()->SetTitle("#frac{d^{2}N}{2#piNP_{T}dP_{T}dy}(GeV/c)^{-2}");
	    gr->GetYaxis()->SetTitleOffset(1.2);
	    gr->SetTitle("#phi Raw Spectra, Au+Au 14.5GeV");
	    gr->Draw("AP"); 

            leg->AddEntry(gr, mCentStringBES[i].c_str(), "p");
	}
        else{
            gr->Draw("P same");
            char leg_string[50];
            sprintf(leg_string, "%s#times 10^{-%d}", mCentStringBES[i].c_str(), 5-i);
            leg->AddEntry(gr, leg_string, "p");
	}
    }

    leg->Draw("sames");

    std::string plotname = "../" + mParticleType + "_plots/" + mParticleType + "_rawspectraBES.pdf";
    rawspectra_can->SaveAs(plotname.c_str());  
}

// No use any more
void StPhiDownMaker::analyzeEff(){
    std::cout << YELLOW << ".... Analyzing efficiency..." << RESET << std::endl;
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
    
    std::cout << "Efficiency analysis done!" << std::endl;
}

std::string StPhiDownMaker::getCentString(Int_t i){
    return mCentString[i];
}

void StPhiDownMaker::plotEffLinear(){
    std::cout << YELLOW << ".... Plotting efficiency on linear scale.." << RESET << std::endl;
    TCanvas* canEff = new TCanvas("canEff", "canEff");
    canEff->SetTicks(1, 1);
    TLegend* leg = new TLegend(0.55, 0.45, 0.75, 0.65);
    leg->SetBorderSize(0);
    for(int j = 0; j < mKCentBin; j++){
        Int_t i =  mKCentBin - 1 - j;
	TGraphErrors* geff = new TGraphErrors(mKPtBin, mXRawSpectra, mEff[i], 0, mEffError[i]);
        geff->SetMarkerSize(1.0);
        geff->SetMarkerStyle(20);
        geff->SetMarkerColor(i+1);
        geff->SetMaximum(0.6);
        geff->SetMinimum(0.0);
        geff->GetXaxis()->SetLimits(0.0, 3.6);
        if(j == 1){
            geff->SetTitle("#Omega^{-} Efficiency, Au+Au 14.5GeV"); 
            geff->GetYaxis()->SetTitle("efficiency");
            geff->GetXaxis()->SetTitle("pT(GeV/c)");
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

// No use any more
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

// No use any more
Double_t StPhiDownMaker::getSpectraWeight(Int_t centbin, Double_t pt){
    Double_t wgt;
    mLevyPt->SetParameters(mLevyPar[centbin]);
    wgt = mLevyPt->Eval(pt);
    return wgt;
}

// No use any more
void StPhiDownMaker::compCorrSpectra(){
    //Compute Data points
    std::cout << YELLOW << ".... Computing corrected spectra data points for 9 cen bins..." << RESET << std::endl;
    for(int i = 0;  i < mKCentBin; i++){
        for(int j = 0; j < mKPtBin; j++){
            mYCorrSpectra[i][j] = mYRawSpectra[i][j] / mEff[i][j]; 
            mYCorrSpectraScale[i][j] = mYRawSpectraScale[i][j] / mEff[i][j]; 

            double relative_rawyerror = mYRawSpectraError[i][j] / mYRawSpectra[i][j];
            double relative_efferror = mEffError[i][j] / mEff[i][j];
            mYCorrSpectraError[i][j] = mYCorrSpectra[i][j] * sqrt(relative_rawyerror*relative_rawyerror + relative_efferror*relative_efferror); 
            mYCorrSpectraErrorScale[i][j] = pow(10, i-mKCentBin+1)*mYCorrSpectraError[i][j];
            //mYCorrSpectraErrorScale[i][j] = mYCorrSpectraScale[i][j] * sqrt(relative_rawyerror*relative_rawyerror + relative_efferror*relative_efferror); 
            //mYCorrSpectraErrorScale[i][j] = pow(10, i+1-mKCentBin) * mYCorrSpectra[i][j] * sqrt(relative_rawyerror*relative_rawyerror); 
            //mYCorrSpectraError[i][j] = mYCorrSpectra[i][j] * sqrt(relative_rawyerror*relative_rawyerror+relative_efferror*relative_efferror); 
            //std::cout << BLUE << "........ Calculation========YCorrSpectraCent" << i << "Pt" << j<< "= " << mYCorrSpectraScale[i][j] << "with error = " << mYCorrSpectraErrorScale[i][j] << " with efficiency being " << "============ re_y = " << relative_rawyerror << "========re_eff" << relative_efferror << std::endl; 
	}
    }

    //Fitting Levy Function, obtain the X position and output the fitting par
    for(int i = 0; i < mKCentBin; i++){
        TGraphErrors* g = new TGraphErrors(mKPtBin-1, mXCorrSpectra[i], mYCorrSpectraScale[i], 0, mYCorrSpectraErrorScale[i]);
        mLevy->SetParameters(mLevyPar[i]);
        g->Fit(mLevy, "QREM0");
        
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

void StPhiDownMaker::compCorrSigCountsBES(){
    std::cout << YELLOW << ".... Computing corrected signal counts for BES cen bins..." << RESET << std::endl;

    for(int i = 0; i < 6; i++){
	for(int j = 0; j < mKPtBin; j++){
	    if(i == 0){// 60-70% + 70-80% = 60-80%
		mCorrSigCountsBES[i][j] = mRawSigCounts[0][j]/mEff[0][j] + mRawSigCounts[1][j]/mEff[1][j];
                Double_t err0 = mRawSigCounts[0][j]/mEff[0][j] * sqrt(pow(mRawSigCountsError[0][j]/mRawSigCounts[0][j], 2) + pow(mEffError[0][j]/mEff[0][j], 2)); 
                Double_t err1 = mRawSigCounts[1][j]/mEff[1][j] * sqrt(pow(mRawSigCountsError[1][j]/mRawSigCounts[1][j], 2) + pow(mEffError[1][j]/mEff[1][j], 2)); 
		mCorrSigCountsBESError[i][j] = sqrt(err0*err0 + err1*err1);
                //std::cout << i << " " << j << " " << mCorrSigCountsBES[i][j] << " <<------------>> " << mCorrSigCountsBESError[i][j] << std::endl;
	    }
	    else if(i == 1){// 50-60% + 40-50% = 40-60%
		mCorrSigCountsBES[i][j] = mRawSigCounts[2][j]/mEff[2][j] + mRawSigCounts[3][j]/mEff[3][j];
                Double_t err0 = mRawSigCounts[2][j]/mEff[2][j] * sqrt(pow(mRawSigCountsError[2][j]/mRawSigCounts[2][j], 2) + pow(mEffError[2][j]/mEff[2][j], 2)); 
                Double_t err1 = mRawSigCounts[3][j]/mEff[3][j] * sqrt(pow(mRawSigCountsError[3][j]/mRawSigCounts[3][j], 2) + pow(mEffError[3][j]/mEff[3][j], 2)); 
		mCorrSigCountsBESError[i][j] = sqrt(err0*err0 + err1*err1);
                //std::cout << i << " " << j << " " << mCorrSigCountsBES[i][j] << " <<------------>> " << mCorrSigCountsBESError[i][j] << std::endl;
	    }
	    else if(i >= 2 && i <=4){// 30-40%, 20-30%, 10-20% 
		mCorrSigCountsBES[i][j] = mRawSigCounts[i+2][j]/mEff[i+2][j];
		mCorrSigCountsBESError[i][j] = mRawSigCounts[i+2][j]/mEff[i+2][j] * sqrt(pow(mRawSigCountsError[i+2][j]/mRawSigCounts[i+2][j], 2) + pow(mEffError[i+2][j]/mEff[i+2][j], 2));
                //std::cout << i << " " << j << " " << mCorrSigCountsBES[i][j] << " <<------------>> " << mCorrSigCountsBESError[i][j] << std::endl;
	    }
	    else{// 0-5% + 5-10% = 0-10%
		mCorrSigCountsBES[i][j] = mRawSigCounts[7][j]/mEff[7][j] + mRawSigCounts[8][j]/mEff[8][j];
                Double_t err0 = mRawSigCounts[7][j]/mEff[7][j] * sqrt(pow(mRawSigCountsError[7][j]/mRawSigCounts[7][j], 2) + pow(mEffError[7][j]/mEff[7][j], 2)); 
                Double_t err1 = mRawSigCounts[8][j]/mEff[8][j] * sqrt(pow(mRawSigCountsError[8][j]/mRawSigCounts[8][j], 2) + pow(mEffError[8][j]/mEff[8][j], 2)); 
		mCorrSigCountsBESError[i][j] = sqrt(err0*err0 + err1*err1);
                //std::cout << i << " " << j << " " << mCorrSigCountsBES[i][j] << " <<------------>> " << mCorrSigCountsBESError[i][j] << std::endl;
	    }

            //std::cout << "===>>>" << i << " " << j << " " << mCorrSigCountsBES[i][j] << "=>BES Corrected Counts == with relative error => " << mCorrSigCountsBESError[i][j]/mCorrSigCountsBES[i][j] << std::endl;
	}
    }
}

void StPhiDownMaker::compCorrSpectraBES(){
    std::cout << YELLOW << ".... Computing corrected spectra data points for BES cen bins..." << RESET << std::endl;
    //Compute Data points
    Double_t PI = 3.1415926;
    for(int i = 0;  i < 6; i++){
        for(int j = 0; j < mKPtBin; j++){
	    mYCorrSpectraBES[i][j] = 1/(2*PI) * mCorrSigCountsBES[i][j] / mXCorrSpectraBES[i][j] / mDptSpectra[j] / mNEventsWeightedBES[i] / mBr;
	    mYCorrSpectraBESScale[i][j] = pow(10, i-5) * mYCorrSpectraBES[i][j];
	    mYCorrSpectraBESError[i][j] = 1/(2*PI) * mCorrSigCountsBESError[i][j] / mXCorrSpectraBES[i][j] / mDptSpectra[j] / mNEventsWeightedBES[i] / mBr;
	    mYCorrSpectraBESErrorScale[i][j] = pow(10, i-5) * mYCorrSpectraBESError[i][j];
	    //std::cout << "mYCorrSpectra for BES cent" << i << "pt" << j <<" is " << mYCorrSpectraBES[i][j] << " mYCorrSpectraError is " << mYCorrSpectraBESError[i][j] << std::endl;
	}
    }

    //Fitting Levy Function, obtain the X position and output the fitting par
    for(int i = 0; i < 6; i++){
        TGraphErrors* g = new TGraphErrors(mKPtBin-2, mXCorrSpectraBES[i], mYCorrSpectraBESScale[i], 0, mYCorrSpectraBESErrorScale[i]); //FIXME: mXRawSpectra
        mLevy->SetParameters(mLevyParBES[i]);
        g->Fit(mLevy, "QREM0");
        
        mLevyParBES[i][0] = mLevy->GetParameter(0);
        mLevyParBESError[i][0] = mLevy->GetParError(0);
        mLevyParBES[i][1] = mLevy->GetParameter(1);
        mLevyParBESError[i][1] = mLevy->GetParError(1);
        mLevyParBES[i][2] = mLevy->GetParameter(2);
        mLevyParBESError[i][2] = mLevy->GetParError(2);
        mLevyPt->SetParameters(mLevyParBES[i]);
        mLevyPt2->SetParameters(mLevyParBES[i]);

	for(int j = 0; j < mKPtBin; j++){
	    mXCorrSpectraBES[i][j] = mLevyPt2->Integral(mPtBd[j], mPtBd[j+1]) / mLevyPt->Integral(mPtBd[j], mPtBd[j+1]);
	}
    }
}

void StPhiDownMaker::compAndPlotRcp(){
    std::cout << YELLOW << ".... Computing and plotting Rcp..." << RESET << std::endl;
    TCanvas* can = new TCanvas("Rcp_can");
    for(int i = 0; i < mKPtBin; i++){
        mRcp[i] = mYCorrSpectra[8][i] * mNCollBES[1] * mXCorrSpectraBES[1][i] / (mXCorrSpectra[8][i] * mNColl[8] * mYCorrSpectraBES[1][i]);
        std::cout << BLUE << "........ Rcp(0-5%/40-60%): pt " << i << " y = " << mRcp[i] << ", mYCorrSpectra[8](0-5%) = " << mYCorrSpectra[8][i] << ", mYCorrSpectra[1](40-60%) = " << mYCorrSpectraBES[1][i] << RESET << std::endl; 
    }

    TGraphErrors* gerr = new TGraphErrors(mKPtBin-2, mXRawSpectra, mRcp, 0, 0);
    gerr->SetMarkerStyle(20);
    gerr->SetTitle("#phi Rcp");
    gerr->GetXaxis()->SetTitle("pT(GeV/c)");
    gerr->GetYaxis()->SetTitle("Rcp(0-5%/40-60%)");
    gerr->Draw("AP");
    gPad->SetTicks(1, 1);

    char plotname[50]; 
    sprintf(plotname, "../%s_plots/Rcp_0_5_40_60.pdf", mParticleType.c_str());
    gPad->SaveAs(plotname);
    delete can;
}

void StPhiDownMaker::compDndy(){
    // Initialize measured and unmeasured pt range
    std::cout << YELLOW << ".... Computing dndy for 9 bins..." << RESET << std::endl;
    Double_t leftlow = 0.;
    Double_t lefthigh = 0.4;
    Double_t rightlow = 5.0;
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
	    mDndy[i] += mYCorrSpectra[i][j] * 2 * 3.1415926 * mDptSpectra[j] * mXCorrSpectra[i][j];
            tmp_dndy_err02 += (mYCorrSpectraError[i][j] * 2 * 3.1415926 * mDptSpectra[j] * mXCorrSpectra[i][j]) * (mYCorrSpectraError[i][j] * 2 * 3.1415926 * mDptSpectra[j] * mXCorrSpectra[i][j]);
	}
        mDndyFit[i] = 2 * 3.1415926 * mLevyPt->Integral(leftlow, righthigh);
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

    TLegend* leg = new TLegend(0.25, 0.05, 0.45, 0.25);
    leg->SetBorderSize(0);

    Int_t i = mKCentBin;
    while(i > 0){
        i--;  
	TGraphErrors* gerr = new TGraphErrors(mKPtBin-1, mXCorrSpectra[i], mYCorrSpectraScale[i], 0, mYCorrSpectraErrorScale[i]); 
        gerr->SetMarkerSize(1.0);
        gerr->SetMarkerStyle(20);
        gerr->SetMarkerColor(i+1);
        char tmp_string[50];
        if(i == (mKCentBin-1)){
	    gerr->SetMaximum(10);
	    gerr->SetMinimum(10E-17);
	    gerr->GetXaxis()->SetLimits(0.0, 3.60);
	    std::string title = "#phi Spectra@AuAu 14.5GeV"; 
	    gerr->SetTitle(title.c_str());
	    gerr->GetYaxis()->SetTitle("#frac{d^{2}N}{2#piNP_{T}dP_{T}dy}(GeV/c)^{-2}");
	    gerr->GetXaxis()->SetTitle("P_{T}(GeV/c)");
	    gerr->GetYaxis()->SetTitleOffset(1.3);
	    gerr->Draw("AP");
	    leg->AddEntry(gerr, mCentString[i].c_str(), "p");
	}
        else{
            gerr->Draw("P same"); 
            sprintf(tmp_string, "%s#times 10^{-%d}", mCentString[i].c_str(), (mKCentBin-i-1)); 
            leg->AddEntry(gerr, tmp_string, "p");
	}
        
        mLevy->SetParameters(mLevyPar[i]);
        TF1* levy_copy = (TF1*)mLevy->Clone();
        levy_copy->SetLineStyle(2);
        levy_copy->SetLineColor(4);
        //levy_copy->Draw("sames");
    }
    leg->Draw("sames");

    char plotname[50]; 
    sprintf(plotname, "../%s_plots/finalCorrSpectra.pdf", mParticleType.c_str());
    canCorrSpectra->SaveAs(plotname);
}

void StPhiDownMaker::plotCorrSpectraBES(){
    std::cout << YELLOW << ".... Start plotting corrected spectra data points for BES centrality bins!" << RESET << std::endl;
    TCanvas* canCorrSpectraBES = new TCanvas("canCorrSpectraBES", "canCorrSpectraBES");
    canCorrSpectraBES->SetLogy();
    canCorrSpectraBES->SetTicks(1, 1);

    TLegend* leg = new TLegend(0.25, 0.15, 0.45, 0.35);
    leg->SetBorderSize(0);

    Int_t i = 6; //mKCentBin;
    while(i > 0){
        i--;  
	TGraphErrors* gerr = new TGraphErrors(mKPtBin-2, mXCorrSpectraBES[i], mYCorrSpectraBESScale[i], 0, mYCorrSpectraBESErrorScale[i]); 
        gerr->SetMarkerSize(1.0);
        gerr->SetMarkerStyle(20);
        gerr->SetMarkerColor(i+1);
        char tmp_string[50];
        if(i == (6-1)){
	    gerr->SetMaximum(10);
	    gerr->SetMinimum(10E-17);
	    gerr->GetXaxis()->SetLimits(0.0, 3.60);
	    std::string title = "#phi Spectra@Au+Au 14.5GeV"; 
	    gerr->SetTitle(title.c_str());
	    gerr->GetYaxis()->SetTitle("#frac{d^{2}N}{2#piNP_{T}dP_{T}dy}(GeV/c)^{-2}");
	    gerr->GetXaxis()->SetTitle("P_{T}(GeV/c)");
	    gerr->GetYaxis()->SetTitleOffset(1.3);
	    gerr->Draw("AP");
	    leg->AddEntry(gerr, mCentStringBES[i].c_str(), "p");
	}
        else{
            gerr->Draw("P same"); 
            sprintf(tmp_string, "%s#times 10^{-%d}", mCentStringBES[i].c_str(), (6-i-1)); 
            leg->AddEntry(gerr, tmp_string, "p");
	}
        
        mLevy->SetParameters(mLevyParBES[i]);
        TF1* levy_copy = (TF1*)mLevy->Clone();
        levy_copy->SetLineStyle(2);
        levy_copy->SetLineColor(4);
        levy_copy->Draw("sames");
    }
    //leg->AddEntry(mLevy, "Levy Function", "l");
    leg->Draw("sames");

    char plotname[50]; 
    sprintf(plotname, "../%s_plots/finalCorrSpectraBES.pdf", mParticleType.c_str());
    canCorrSpectraBES->SaveAs(plotname);
}

void StPhiDownMaker::printYieldsBES(){
    std::cout << YELLOW << ".... Printing yields for BES Centralities!" << RESET << std::endl;

    for(int i = 0; i < 6; i++){
        for(int j = 0; j < mKPtBin; j++){
            //std::cout << BLUE << "........ Yields for cen" <<i << " pt" << j << ": " << mCorrSigCountsBES[i][j] / mNEventsWeightedBES[i] << " " << mCorrSigCountsBESError[i][j] << RESET << std::endl;
	}
    }

    // Print Data for omega ratio
    TF1* levypt_copy = (TF1*) mLevyPt->Clone();
    TF1* levypt2_copy = (TF1*) mLevyPt2->Clone();
    levypt_copy->SetParameters(mLevyParBES[5]);
    levypt2_copy->SetParameters(mLevyParBES[5]);
    
    // Compute the Phi Yields for 0-10% only
    for(int i = 0; i < 6; i++){
        double x_pos = levypt2_copy->Integral(mPtBdOmg[i], mPtBdOmg[i+1])/levypt_copy->Integral(mPtBdOmg[i], mPtBdOmg[i+1]);
        mPhiYields010[i] = 2 * 3.1415926 * levypt_copy->Integral(mPtBdOmg[i], mPtBdOmg[i+1]);
	std::cout << BLUE << "........ Printing phi yields for 0-10%: pt" << i << ", x_position = " << x_pos << ", mPhiYieds010 = " << mPhiYields010[i] << ", yield010_counting = " <<  mCorrSigCountsBES[5][i] / mNEventsWeightedBES[5] << RESET <<std::endl;
    }
}
/*
void StPhiDownMaker::compYieldsBES(){
    for(int i = 0; i < 6; i++){
        for(int j = 0; j < mKPtBin; j++){
            mYields[i][j] = m
	}
    }
}
*/
void StPhiDownMaker::Analyze(){
    std::cout << "Load infile_dat/rot successfully!" << std::endl;
    for(int i = 0; i < mKCentBin; i++){
    //for(int i = 0; i < 1; i++){
	for(int j = 0; j < mKPtBin; j++){
	//for(int j = 0; j < 2; j++){
	    char hist_name_dat[50];
	    char hist_name_mix[50]; 
	    sprintf(hist_name_dat, "sig_phipt%dcent%d", j+1, i+1); 
	    sprintf(hist_name_mix, "phipt%dcent%d", j+1, i+1);

	    TH1F* hdat = (TH1F*)mDatFile->Get(hist_name_dat);
	    TH1F* hmix = (TH1F*)mDatFile->Get(hist_name_mix);

	    hmix->Sumw2();
	    hdat->Sumw2();

	    double mix_scale = compMixNormFactor(i, j, hdat, hmix);
            double bin_width = hdat->GetBinWidth(50);

	    plotMixInvMassWithData(i, j, hdat, hmix, mix_scale);
	    plotInvMassAfterBgSubtraction(i, j, hdat, hmix, mix_scale);

	    // Count the raw number of phi candidates
	    compRawSigCounts(i, j, bin_width); 
	}
    }
    
    // Compute the raw spectra of phi
    compRawSpectra(); 

    // Plot the raw spectra of phi
    plotRawSpectra();

    // Compute the corrected spectra of phi with efficiency correction
    compCorrSpectra();
    compCorrSpectra();
    compCorrSpectra();

    // Plot the corrected spectra of phi
    plotCorrSpectra();

    // Compute the dndy of phi 
    compDndy();

/*
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
*/
}

void StPhiDownMaker::plotOmegaPhiRatio(){
    std::cout << YELLOW << ".... Plotting omega/phi ratio" << RESET << std::endl;
    TCanvas* can = new TCanvas();
    for(int i = 0; i < 6; i++) {
	float tempOmgYield = 0.5*(mOmgYields010[i] + mAntiOmgYields010[i]);
        mOmgPhiRatio010[i] = 0.5*(mOmgYields010[i] + mAntiOmgYields010[i]) / mPhiYields010[i];
	std::cout << BLUE << "........ In 0-10% most central collision, omega(omegabar) yield for pt" << i << " is " << tempOmgYield << ", phi yield is " << mPhiYields010[i] << RESET << std::endl;
    }

    TGraph* graph = new TGraph(6, mXRawSpectraOmg, mOmgPhiRatio010);
    graph->SetMarkerStyle(20);
    graph->Draw("AP");

    char plotname[50];
    sprintf(plotname, "../%s_plots/omg_phi_ratio.pdf", mParticleType.c_str());
    gPad->SaveAs(plotname);
}

void StPhiDownMaker::plotOmgPhiSpectra010(){
    Double_t omgCorrSpectra010[6] = {0.0102435, 0.00654648, 0.00260305, 0.000729941, 0.00018738, 3.37409e-05};
    Double_t omgbarCorrSpectra010[6] = {0.00435631, 0.00251424, 0.00091688, 0.000424534, 0.000155178, 2.48972e-05};
    Double_t omgCorrSpectraError010[6] = {0.00557996, 0.000921012, 0.000285038, 0.000106882, 4.09707e-05, 8.40746e-06};
    Double_t omgbarCorrSpectraError010[6] = {0.0018596, 0.000409823, 0.000134689, 5.66398e-05, 2.56832e-05, 5.58363e-06};

    TGraphErrors* gerr = new TGraphErrors(6, mXRawSpectraOmg, omgCorrSpectra010, 0, omgCorrSpectraError010);
    gPad->SetLogy();
    gerr->SetTitle("Omg/AntiOmg(Red/Green) vs Phi(Black)0-10%");
    gerr->SetMaximum(10);
    gerr->SetMinimum(10E-11);
    gerr->GetXaxis()->SetLimits(0, 4.0);
    gerr->GetXaxis()->SetTitle("pT(GeV/c)");
    gerr->GetYaxis()->SetTitle("#frac{d^{2}N}{2#piNP_{T}dP_{T}dy}(GeV/c)^{-2}");
    gerr->SetMarkerColor(2);
    gerr->SetMarkerStyle(4);
    gerr->Draw("AP");

    TGraphErrors* gerr0 = new TGraphErrors(6, mXRawSpectraOmg, omgbarCorrSpectra010, 0, omgbarCorrSpectraError010);
    gerr0->SetMarkerColor(3);
    gerr0->SetMarkerStyle(4);
    gerr0->Draw("P same");

    TGraphErrors* gerr1 = new TGraphErrors(mKPtBin, mXCorrSpectraBES[5], mYCorrSpectraBES[5], 0, mYCorrSpectraBESError[5]);
    gerr1->SetMarkerColor(1);
    gerr1->SetMarkerStyle(8);
    gerr1->Draw("P same");
   
    gPad->SetTicks(1, 1);
    gPad->SaveAs("../Phi_plots/omgphiCorrSpectraComparison010.pdf");
}

void StPhiDownMaker::compare11GeVRawSpectra010(){
    Double_t nEvents11GeV010 = 2.3034924e+6;
    Double_t phiRawYield11GeV010[11] = {646.518, 3357.29, 9750.72, 15671.5, 50887.3, 27077.5, 60889.9, 36305.7, 10070.9, 9152.78, 2189.22};
    Double_t phiPt11GeV010[11] = {.35, .45, .55, .65, .8, .95, 1.15, 1.5, 1.85, 2.25, 3.0};
    Double_t phiDpt11GeV010[11] = {.1, .1, .1, .1, .2, .1, .3, .4, .3, .5, 1.0};
    Double_t phiYRawSpectra11GeV010[11] = {0, };

    for(int i = 0; i < 11; i++){
        phiYRawSpectra11GeV010[i] = phiRawYield11GeV010[i] / (2 * 3.1415926 * phiPt11GeV010[i] * phiDpt11GeV010[i] * nEvents11GeV010);
    }

}
// The main analysis member function 
void StPhiDownMaker::AnalyzeBES(){
    std::cout << CYAN << ">> Analyzing data..." << std::endl;

    for(int i = 0; i < mKCentBin; i++){
	for(int j = 0; j < mKPtBin; j++){
	    char hist_name_dat[50];
	    char hist_name_mix[50];
	    sprintf(hist_name_dat, "sig_phipt%dcent%d", j+1, i+1);
	    sprintf(hist_name_mix, "phipt%dcent%d", j+1, i+1);

	    TH1F* hdat = (TH1F*)mDatFile->Get(hist_name_dat);
	    TH1F* hmix = (TH1F*)mDatFile->Get(hist_name_mix);

	    hmix->Sumw2();
	    hdat->Sumw2();

	    double mix_scale = compMixNormFactor(i, j, hdat, hmix);
            double bin_width = hdat->GetBinWidth(50);

	    plotMixInvMassWithData(i, j, hdat, hmix, mix_scale);
	    plotInvMassAfterBgSubtraction(i, j, hdat, hmix, mix_scale);

	    // Compute the raw counting number of phi candidates  
	    compRawSigCounts(i, j, bin_width); 
	}
    }

    // Plot invariant mass QA plots
    plotInvMassQA();

    // Plot raw spectra for 6 centrality bins
    compRawSigCountsBES();

    // Compute raw spectra for 6 centrality bins
    compRawSpectraBES(); 

    // Plot raw spectra for 6 centrality bins
    plotRawSpectraBES();

    // Compute the raw spectra for 10 centrality bins
    compRawSpectra();
     
    // Compute the corrected spectra for 10 centrality bins 
    compCorrSpectra();
    compCorrSpectra();
    compCorrSpectra();

    for(int i = 0; i < mKCentBin; i++){
        for(int j = 0; j < mKPtBin; j++) { 
	std::cout << BLUE << "........ CorrSpectra: cen" << i << " pt" << j<< " y = " << mYCorrSpectraScale[i][j] << ", yerr = " << mYCorrSpectraErrorScale[i][j] << RESET << std::endl; 
	}
    }
    // Compute the corrected phi candidates number with efficiency correction for 6 centrality bins;
    compCorrSigCountsBES(); 

    // Compute corrected spectra for 6 centrality bins for 3 times to make the fitting coverge
    compCorrSpectraBES();
    compCorrSpectraBES();
    compCorrSpectraBES();

    for(int i = 0; i < 6; i++){
        for(int j = 0; j < mKPtBin; j++){
	std::cout << BLUE << "........ CorrSPectraBES: cen" << i << " pt" << j<< " y = " << mYCorrSpectraScale[i][j] << ", yerr = " << mYCorrSpectraErrorScale[i][j] << RESET << std::endl; 
	}
    }
    // Plot the efficiency data of Phi
    plotEffLinear();

    // Plot the corrected spectra for 6 centrality bins
    plotCorrSpectraBES();

    // Print out the phi yields for 6 centrality bins
    printYieldsBES();

    // Compute and plot the Rcp of Phi
    compAndPlotRcp();

    // Plot the Omega/Phi ratio
    plotOmegaPhiRatio();

    plotOmgPhiSpectra010();
}
