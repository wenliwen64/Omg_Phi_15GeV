#ifndef StPhiDownMaker_H
#define StPhiDownMaker_H
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TObject.h"
#include "TProfile.h" 
#include <string>
class StPhiDownMaker: public TObject{
    TFile* mDatFile;
 
    //Fitting Function
    TF1* mLevy;
    TF1* mLevyPt;
    TF1* mLevyPt2;
    TF1* mTotal;
    TF1* mBW;
    TF1* mPolyBg;

    TH1F* mHFpEffFine[9];
    TH1F* mHExpEffFine[9];

    std::string mCentString[9]; 
   
    std::string mParticleType;
    Double_t pdgmass_phi;
    Int_t mKCentBin;
    Int_t mKCentBin2; // This no of centrality bin is used for Physics Plots
    Int_t mKPtBin;
    Double_t mNEventsUnweighted[9];
    Double_t mNEventsWeighted[9];
    Double_t mBr;
    Double_t mSigRangeLeft;
    Double_t mSigRangeRight;

    Double_t mMixNormLeftLowB;
    Double_t mMixNormLeftHighB;
    Double_t mMixNormRightLowB;
    Double_t mMixNormRightHighB;

    Double_t mPtBd[12];

    Double_t mFpEff[9][11];
    Double_t mFpEffError[9][11];
    Double_t mExpEff[9][11];
    Double_t mExpEffError[9][11];
    Double_t mEff[9][11];
    Double_t mEffError[9][11];

    Double_t mMixScale_ratio[9][11];

    Double_t mRawSigCounts[9][11];
    Double_t mRawSigCountsError[9][11];
    Double_t mXRawSpectra[11];
    Double_t mXRawSpectraError[11];
    Double_t mYRawSpectra[9][11];
    Double_t mYRawSpectraScale[9][11];
    Double_t mYRawSpectraError[9][11];
    Double_t mYRawSpectraErrorScale[9][11];

    Double_t mDptSpectra[11];

    Double_t mXCorrSpectra[9][11];
    Double_t mYCorrSpectra[9][11];
    Double_t mYCorrSpectraScale[9][11];
    Double_t mYCorrSpectraError[9][11];
    Double_t mYCorrSpectraErrorScale[9][11];

    Double_t mLevyPar[9][3]; 
    Double_t mLevyParError[9][3]; 
    Double_t mInvMassPar[9][11][3]; 
    Double_t mInvMassParError[9][11][3]; 
    Double_t mDndy[9];
    Double_t mDndyError[9];
    Double_t mDndyFit[9];

    void effInit();
    void levyInit();
    void bwFuncInit();
    void nEventsInit();
    void mixBgAnalysisInit();
    Double_t compMixNormFactor(Int_t centbin, Int_t ptbin, TH1F* hdat, TH1F* hrot);
    void plotMixInvMassWithData(Int_t centbin, Int_t ptbin, TH1F* hdat, TH1F* hrot, Double_t scale); 
    void plotInvMassAfterBgSubtraction(Int_t centbin, Int_t ptbin, TH1F* hdat, TH1F* hrot, Double_t scale); 
    void compRawSigCounts(Int_t centbin, Int_t ptbin, Double_t bin_width);
    void compRawSpectra();
    void plotRawSpectra();
    void analyzeEff(); // input: eff, levy; Load efficiency raw data file and apply the cuts; output efficiency data; iteratively compute the data points
    std::string getCentString(Int_t);
    void plotEff();
    Double_t getSpectraWeight(Int_t, Double_t);
    void compCorrSpectra(); // input: eff, xpos, levy
    void compDndy(); // Integrate the unmeasured range and add up the measured range.
    Double_t getDndy(Int_t);
    Double_t getDndyError(Int_t);
    void plotCorrSpectra();
public:
    StPhiDownMaker(std::string par_type);
    ~StPhiDownMaker();
    void Init(std::string datfile);
    void Analyze(); // call analyzeEff(); compCorrSpectra(); analyzeEff()

    ClassDef(StPhiDownMaker, 1)
};
#endif
