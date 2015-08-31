#ifndef StrAnalyMaker_H
#define StrAnalyMaker_H
#include "TH1F.h"
#include "TFile.h"
#include "TObject.h"
#include <string>
class StrAnalyMaker: public TObject{
    TFile* mOverviewFile;
    TFile* mDatFile;
    TFile* mRotBgFile;
    TFile* mFpEffFile;
    TFile* mExpEffFile;
    TFile* mLevy;
    TFile* mLevyPt;
    TFile* mLevyPt2;

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

    Double_t mPtBD[7];

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

    Double_t mLevyPar[2][3]; 
    Double_t mDndy[2][6];

    void rotBgAnalysisInit();
    Double_t compRotNormFactor(Int_t centbin, Int_t ptbin, TH1F* hdat, TH1F* hrot);
public:
    StrAnalyMaker(std::string par_type);
    ~StrAnalyMaker();
    void Init(std::string overveiwfile, std::string datfile, std::string rotfile, std::string fpefffile, std::string expefffile);
    void plotRotInvMassWithData(Int_t centbin, Int_t ptbin, TH1F* hdat, TH1F* hrot, Double_t scale); 
    void compRawSigCounts(Int_t centbin, Int_t ptbin, TH1F* hdat, TH1F* hrot, Double_t scale);
    void compRawSpectra();
    void plotRawSpectra();
    void analyzeEff(); // input: eff, levy; Load efficiency raw data file and apply the cuts; output efficiency data; iteratively compute the data points
    void compCorrSpectra(); // input: eff, xpos, levy
    void compDndy(); // Integrate the unmeasured range and add up the measured range.
    void Analyze(); // call analyzeEff(); compCorrSpectra(); analyzeEff()

    Double_t getDndy(Int_t centbin);
    Double_t getDndyError(Int_t centbin);
    ClassDef(StrAnalyMaker, 1)
};
#endif
