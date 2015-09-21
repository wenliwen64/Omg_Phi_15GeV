#ifndef StrAnalyMaker_H
#define StrAnalyMaker_H
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TObject.h"
#include "TProfile.h" 
#include <string>
class StrAnalyMaker: public TObject{
    TFile* mOverviewFile;
    TFile* mDatFile;
    TFile* mRotBgFile;
    TFile* mFpEffFile;
    TFile* mExpEffFile;
    TF1* mLevy;
    TF1* mLevyPt;
    TF1* mLevyPt2;

    TH1F* mHFpEffFine[2];
    TH1F* mHExpEffFine[2];

    std::string mCentString[2]; 
   
    std::string mParticleType;
    Double_t pdgmass_xi;
    Int_t mKCentBin;
    Int_t mKPtBin;
    Double_t mNEventsUnweighted[2];
    Double_t mNEventsWeighted[2];
    Double_t mBr;
    Double_t mSigRangeLeft;
    Double_t mSigRangeRight;

    Double_t mRotNormLeftLowB;
    Double_t mRotNormLeftHighB;
    Double_t mRotNormRightLowB;
    Double_t mRotNormRightHighB;

    Double_t mPtBd[7];

    Double_t mFpEff[2][6];
    Double_t mFpEffError[2][6];
    Double_t mExpEff[2][6];
    Double_t mExpEffError[2][6];
    Double_t mEff[2][6];
    Double_t mEffError[2][6];

    Double_t mRotScale_ratio[2][6];

    Double_t mRawSigCounts[2][6];
    Double_t mRawSigCountsError[2][6];
    Double_t mXRawSpectra[6];
    Double_t mXRawSpectraError[6];
    Double_t mYRawSpectra[2][6];
    Double_t mYRawSpectraError[2][6];

    Double_t mDptSpectra[6];
    Double_t mXCorrSpectra[2][6];
    Double_t mYCorrSpectra[2][6];
    Double_t mYCorrSpectraError[2][6];
    Double_t mYields[2][6];
    Double_t mYieldsError[2][6];
//11GeV comparison data
    Double_t mDptSpectra11GeV[6];
    Double_t mXCorrSpectra11GeV[2][6];
    Double_t mYCorrSpectra11GeV[2][6];
    Double_t mYCorrSpectra11GeVError[2][6];
    Double_t mYields11GeV[2][6];
    Double_t mYields11GeVError[2][6];

    Double_t mLevyPar[2][3]; 
    Double_t mLevyParError[2][3]; 
    Double_t mDndy[2];
    Double_t mDndyError[2];
    Double_t mDndyFit[2];

    void effInit();
    void levyInit();
    void nEventsInit();
    void rotBgAnalysisInit();
    Double_t compRotNormFactor(Int_t centbin, Int_t ptbin, TH1F* hdat, TH1F* hrot);
    void plotRotInvMassWithData(Int_t centbin, Int_t ptbin, TH1F* hdat, TH1F* hrot, Double_t scale); 
    void compRawSigCounts(Int_t centbin, Int_t ptbin, TH1F* hdat, TH1F* hrot, Double_t scale);
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
    void compYields();
    void compare11GeV();
public:
    StrAnalyMaker(std::string par_type);
    ~StrAnalyMaker();
    void Init(std::string overveiwfile, std::string datfile, std::string rotfile, std::string fpefffile, std::string expefffile);
    void Analyze(); // call analyzeEff(); compCorrSpectra(); analyzeEff()

    ClassDef(StrAnalyMaker, 1)
};
#endif
