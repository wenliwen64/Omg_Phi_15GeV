{
    ifstream scalerot_dat("/Users/lwen/Documents/Omg_Phi_14GeV/Omega/plot_scripts/scale_rot.dat");
    gStyle->SetOptStat(0);
    const Float_t pdgV0Mass = 1.11568;
    const Float_t pdgXiMass = 1.67245;
    const Float_t masswidth = 0.07 ;
    const Int_t kCentBin = 2;
    const Int_t kPtBin = 6;
    const Float_t ptbd[kPtBin+1] = {0.7, 1.2, 1.6, 2.0, 2.4, 2.8, 3.6};
    const Float_t centbd[kCentBin+1] = {1.5, 6.5, 8.5};

    //TFile* mc_file = new TFile("mcomg_fp0.coarse_ptbin_test.cuts.histo.root", "read");
    TFile* mc_file = new TFile("mcomg_exp0.coarse_ptbin_test.cuts.histo.root", "read");
    TFile* data_file = new TFile("test_hist_embedding_15GeV.root", "read");
    TFile* rot_file = new TFile("testrot_hist_embedding_15GeV.root", "read");
    
    //==== Read Data Hist. ====
    TH1F* hmDau1nHits_data[kCentBin][kPtBin];
    TH1F* hmDau2nHits_data[kCentBin][kPtBin];
    TH1F* hmBachnHits_data[kCentBin][kPtBin];
    TH1F* hmDau1Pt_data[kCentBin][kPtBin];
    TH1F* hmDau2Pt_data[kCentBin][kPtBin];
    TH1F* hmBachPt_data[kCentBin][kPtBin];

    TH1F* hmDau1Dca_data[kCentBin][kPtBin];
    TH1F* hmDau2Dca_data[kCentBin][kPtBin];
    TH1F* hmDca1to2_data[kCentBin][kPtBin];
    TH1F* hmV0Dca_data[kCentBin][kPtBin]; 
    TH1F* hmV0DecLen_data[kCentBin][kPtBin];
    TH1F* hmV0InvMass_data[kCentBin][kPtBin];

    TH1F* hmBachDca_data[kCentBin][kPtBin];
    TH1F* hmDcav0tobach_data[kCentBin][kPtBin];
    TH1F* hmXiDca_data[kCentBin][kPtBin];
    TH1F* hmXiDecLen_data[kCentBin][kPtBin];
    TH1F* hmDecLenDiff_data[kCentBin][kPtBin];
    TH1F* hmXiInvMass_data[kCentBin][kPtBin];
    TH1F* hmXiSinth_data[kCentBin][kPtBin];

    //==== Read Rot Hist. ====
    TH1F* hmDau1nHits_rot[kCentBin][kPtBin];
    TH1F* hmDau2nHits_rot[kCentBin][kPtBin];
    TH1F* hmBachnHits_rot[kCentBin][kPtBin];
    TH1F* hmDau1Pt_rot[kCentBin][kPtBin];
    TH1F* hmDau2Pt_rot[kCentBin][kPtBin];
    TH1F* hmBachPt_rot[kCentBin][kPtBin];

    TH1F* hmDau1Dca_rot[kCentBin][kPtBin];
    TH1F* hmDau2Dca_rot[kCentBin][kPtBin];
    TH1F* hmDca1to2_rot[kCentBin][kPtBin];
    TH1F* hmV0Dca_rot[kCentBin][kPtBin]; 
    TH1F* hmV0DecLen_rot[kCentBin][kPtBin];
    TH1F* hmV0InvMass_rot[kCentBin][kPtBin];

    TH1F* hmBachDca_rot[kCentBin][kPtBin];
    TH1F* hmDcav0tobach_rot[kCentBin][kPtBin];
    TH1F* hmXiDca_rot[kCentBin][kPtBin];
    TH1F* hmXiDecLen_rot[kCentBin][kPtBin];
    TH1F* hmDecLenDiff_rot[kCentBin][kPtBin];
    TH1F* hmXiInvMass_rot[kCentBin][kPtBin];
    TH1F* hmXiSinth_rot[kCentBin][kPtBin];

    //==== MC Data Hist. ====
    TH1F* hmDau1nHits_mc[kCentBin][kPtBin];
    TH1F* hmDau2nHits_mc[kCentBin][kPtBin];
    TH1F* hmBachnHits_mc[kCentBin][kPtBin];
    TH1F* hmDau1Pt_mc[kCentBin][kPtBin];
    TH1F* hmDau2Pt_mc[kCentBin][kPtBin];
    TH1F* hmBachPt_mc[kCentBin][kPtBin];

    TH1F* hmDau1Dca_mc[kCentBin][kPtBin];
    TH1F* hmDau2Dca_mc[kCentBin][kPtBin];
    TH1F* hmDca1to2_mc[kCentBin][kPtBin];
    TH1F* hmV0Dca_mc[kCentBin][kPtBin]; 
    TH1F* hmV0DecLen_mc[kCentBin][kPtBin];
    TH1F* hmV0InvMass_mc[kCentBin][kPtBin];

    TH1F* hmBachDca_mc[kCentBin][kPtBin];
    TH1F* hmDcav0tobach_mc[kCentBin][kPtBin];
    TH1F* hmXiDca_mc[kCentBin][kPtBin];
    TH1F* hmXiDecLen_mc[kCentBin][kPtBin];
    TH1F* hmDecLenDiff_mc[kCentBin][kPtBin];
    TH1F* hmXiInvMass_mc[kCentBin][kPtBin];
    TH1F* hmXiSinth_mc[kCentBin][kPtBin];

    for(Int_t iCent = 0; iCent < kCentBin; iCent++){
	for(Int_t iPt = 0; iPt < kPtBin; iPt++){
            Float_t scale_ratio = 0.; 
            Float_t dummy = 0.;
            scalerot_dat >> dummy >> dummy >> scale_ratio; 
            cout << scale_ratio << endl;
	    TString hDau1nHitsName;
	    TString hDau2nHitsName;
            TString hBachnHitsName;
	    TString hDau1PtName;
	    TString hDau2PtName;
            TString hBachPtName;

	    TString hDau1DcaName; 
	    TString hDau2DcaName; 
	    TString hDca1to2Name; 
	    TString hV0DcaName; 
	    TString hV0DecLenName;
	    TString hV0InvMassName;
	    //TString hV0ColinearName;

	    TString hBachDcaName; 
	    TString hDcav0tobachName;
	    TString hXiDcaName;
	    TString hXiDecLenName;
	    TString hDecLenDiffName;
	    TString hXiInvMassName;
	    //TString hXiColinearName;
            TString hXiSinthName;

	    hDau1nHitsName.Form("hmDau1nHits_cen%dpt%d", iCent, iPt);
	    hDau2nHitsName.Form("hmDau2nHits_cen%dpt%d", iCent, iPt);
            hBachnHitsName.Form("hmBachnHits_cen%dpt%d", iCent, iPt);
	    hDau1PtName.Form("hmDau1Pt_cen%dpt%d", iCent, iPt);
	    hDau2PtName.Form("hmDau2Pt_cen%dpt%d", iCent, iPt);
            hBachPtName.Form("hmBachPt_cen%dpt%d", iCent, iPt);

	    hDau1DcaName.Form("hmDau1Dca_cen%dpt%d", iCent, iPt);
	    hDau2DcaName.Form("hmDau2Dca_cen%dpt%d", iCent, iPt);
	    hDca1to2Name.Form("hmDca1to2_cen%dpt%d", iCent, iPt);
	    hV0DcaName.Form("hmV0Dca_cen%dpt%d", iCent, iPt);
	    hV0DecLenName.Form("hmV0DecLen_cen%dpt%d", iCent, iPt);
	    hV0InvMassName.Form("hmV0Invmass_cen%dpt%d", iCent, iPt);
	    //hV0ColinearName.Form("hmV0Colinear_cen%dpt%d", iCent, iPt);

	    hBachDcaName.Form("hmBachDca_cen%dpt%d", iCent, iPt);
	    hDcav0tobachName.Form("hmDcav0tobach_cen%dpt%d", iCent, iPt);
	    hXiDcaName.Form("hmXiDca_cen%dpt%d", iCent, iPt);
	    hXiDecLenName.Form("hmXiDecLen_cen%dpt%d", iCent, iPt);
	    hDecLenDiffName.Form("hmDecLenDiff_cen%dpt%d", iCent, iPt);
	    hXiInvMassName.Form("hmXiInvMass_cen%dpt%d", iCent, iPt);
	    //hXiColinearName.Form("hmXiColinear_cen%dpt%d", iCent, iPt);
	    hXiSinthName.Form("hmXisinth_cen%dpt%d", iCent, iPt);
	    cout<<"happy ipt = "<<iPt<<endl;

            //==== Get Real Data Hist. ====
	    hmDau1nHits_data[iCent][iPt] = (TH1F*)data_file->Get(hDau1nHitsName->Data());
	    hmDau2nHits_data[iCent][iPt] = (TH1F*)data_file->Get(hDau2nHitsName->Data()); 
            hmBachnHits_data[iCent][iPt] = (TH1F*)data_file->Get(hBachnHitsName->Data());
	    hmDau1Pt_data[iCent][iPt] = (TH1F*)data_file->Get(hDau1PtName->Data()); 
	    hmDau2Pt_data[iCent][iPt] = (TH1F*)data_file->Get(hDau2PtName->Data()); 
            hmBachPt_data[iCent][iPt] = (TH1F*)data_file->Get(hBachPtName->Data()); 

	    hmDau1Dca_data[iCent][iPt] = (TH1F*)data_file->Get(hDau1DcaName->Data());
	    hmDau2Dca_data[iCent][iPt] = (TH1F*)data_file->Get(hDau2DcaName->Data());
	    hmDca1to2_data[iCent][iPt] = (TH1F*)data_file->Get(hDca1to2Name->Data());
	    hmV0Dca_data[iCent][iPt] = (TH1F*)data_file->Get(hV0DcaName->Data());
	    hmV0DecLen_data[iCent][iPt] = (TH1F*)data_file->Get(hV0DecLenName->Data()); 
	    hmV0InvMass_data[iCent][iPt] = (TH1F*)data_file->Get(hV0InvMassName->Data());

	    hmBachDca_data[iCent][iPt] = (TH1F*)data_file->Get(hBachDcaName->Data());
	    hmDcav0tobach_data[iCent][iPt] = (TH1F*)data_file->Get(hDcav0tobachName->Data());
	    hmXiDca_data[iCent][iPt] = (TH1F*)data_file->Get(hXiDcaName->Data());
	    hmXiDecLen_data[iCent][iPt] = (TH1F*)data_file->Get(hXiDecLenName->Data()); 
	    hmDecLenDiff_data[iCent][iPt] = (TH1F*)data_file->Get(hDecLenDiffName->Data());
	    hmXiInvMass_data[iCent][iPt] = (TH1F*)data_file->Get(hXiInvMassName->Data());
	    hmXiSinth_data[iCent][iPt] = (TH1F*)data_file->Get(hXiSinthName->Data());

	    hmDau1nHits_data[iCent][iPt]->Sumw2();
	    hmDau2nHits_data[iCent][iPt]->Sumw2();
            hmBachnHits_data[iCent][iPt]->Sumw2();
	    hmDau1Pt_data[iCent][iPt]->Sumw2();
	    hmDau2Pt_data[iCent][iPt]->Sumw2();
            hmBachPt_data[iCent][iPt]->Sumw2();

	    hmDau1Dca_data[iCent][iPt]->Sumw2();
	    hmDau2Dca_data[iCent][iPt]->Sumw2();
	    hmDca1to2_data[iCent][iPt]->Sumw2();
	    hmV0Dca_data[iCent][iPt]->Sumw2();
	    hmV0DecLen_data[iCent][iPt]->Sumw2();
	    hmV0InvMass_data[iCent][iPt]->Sumw2();

	    hmBachDca_data[iCent][iPt]->Sumw2();
	    hmDcav0tobach_data[iCent][iPt]->Sumw2();
	    hmXiDca_data[iCent][iPt]->Sumw2();
	    hmXiDecLen_data[iCent][iPt]->Sumw2();
	    hmDecLenDiff_data[iCent][iPt]->Sumw2();
	    hmXiInvMass_data[iCent][iPt]->Sumw2();
	    hmXiSinth_data[iCent][iPt]->Sumw2();
            
            //==== Get Rot Data Hist. ====
	    hmDau1nHits_rot[iCent][iPt] = (TH1F*)rot_file->Get(hDau1nHitsName->Data());
	    hmDau2nHits_rot[iCent][iPt] = (TH1F*)rot_file->Get(hDau2nHitsName->Data()); 
            hmBachnHits_rot[iCent][iPt] = (TH1F*)rot_file->Get(hBachnHitsName->Data());
	    hmDau1Pt_rot[iCent][iPt] = (TH1F*)rot_file->Get(hDau1PtName->Data()); 
	    hmDau2Pt_rot[iCent][iPt] = (TH1F*)rot_file->Get(hDau2PtName->Data()); 
            hmBachPt_rot[iCent][iPt] = (TH1F*)rot_file->Get(hBachPtName->Data()); 

	    hmDau1Dca_rot[iCent][iPt] = (TH1F*)rot_file->Get(hDau1DcaName->Data());
	    hmDau2Dca_rot[iCent][iPt] = (TH1F*)rot_file->Get(hDau2DcaName->Data());
	    hmDca1to2_rot[iCent][iPt] = (TH1F*)rot_file->Get(hDca1to2Name->Data());
	    hmV0Dca_rot[iCent][iPt] = (TH1F*)rot_file->Get(hV0DcaName->Data());
	    hmV0DecLen_rot[iCent][iPt] = (TH1F*)rot_file->Get(hV0DecLenName->Data()); 
	    hmV0InvMass_rot[iCent][iPt] = (TH1F*)rot_file->Get(hV0InvMassName->Data());

	    hmBachDca_rot[iCent][iPt] = (TH1F*)rot_file->Get(hBachDcaName->Data());
	    hmDcav0tobach_rot[iCent][iPt] = (TH1F*)rot_file->Get(hDcav0tobachName->Data());
	    hmXiDca_rot[iCent][iPt] = (TH1F*)rot_file->Get(hXiDcaName->Data());
	    hmXiDecLen_rot[iCent][iPt] = (TH1F*)rot_file->Get(hXiDecLenName->Data()); 
	    hmDecLenDiff_rot[iCent][iPt] = (TH1F*)rot_file->Get(hDecLenDiffName->Data());
	    hmXiInvMass_rot[iCent][iPt] = (TH1F*)rot_file->Get(hXiInvMassName->Data());
	    hmXiSinth_rot[iCent][iPt] = (TH1F*)rot_file->Get(hXiSinthName->Data());

	    hmDau1nHits_rot[iCent][iPt]->Sumw2();
	    hmDau2nHits_rot[iCent][iPt]->Sumw2();
            hmBachnHits_rot[iCent][iPt]->Sumw2();
	    hmDau1Pt_rot[iCent][iPt]->Sumw2();
	    hmDau2Pt_rot[iCent][iPt]->Sumw2();
            hmBachPt_rot[iCent][iPt]->Sumw2();

	    hmDau1Dca_rot[iCent][iPt]->Sumw2();
	    hmDau2Dca_rot[iCent][iPt]->Sumw2();
	    hmDca1to2_rot[iCent][iPt]->Sumw2();
	    hmV0Dca_rot[iCent][iPt]->Sumw2();
	    hmV0DecLen_rot[iCent][iPt]->Sumw2();
	    hmV0InvMass_rot[iCent][iPt]->Sumw2();

	    hmBachDca_rot[iCent][iPt]->Sumw2();
	    hmDcav0tobach_rot[iCent][iPt]->Sumw2();
	    hmXiDca_rot[iCent][iPt]->Sumw2();
	    hmXiDecLen_rot[iCent][iPt]->Sumw2();
	    hmDecLenDiff_rot[iCent][iPt]->Sumw2();
	    hmXiInvMass_rot[iCent][iPt]->Sumw2();
	    hmXiSinth_rot[iCent][iPt]->Sumw2();
            
/*
	    hmDau1nHits_rot[iCent][iPt]->Scale(1./scale_ratio);
	    hmDau2nHits_rot[iCent][iPt]->Scale(1./scale_ratio);
            hmBachnHits_rot[iCent][iPt]->Scale(1./scale_ratio);
	    hmDau1Pt_rot[iCent][iPt]->Scale(1./scale_ratio);
	    hmDau2Pt_rot[iCent][iPt]->Scale(1./scale_ratio);
            hmBachPt_rot[iCent][iPt]->Scale(1./scale_ratio);

	    hmDau1Dca_rot[iCent][iPt]->Scale(1./scale_ratio);
	    hmDau2Dca_rot[iCent][iPt]->Scale(1./scale_ratio);
	    hmDca1to2_rot[iCent][iPt]->Scale(1./scale_ratio);
	    hmV0Dca_rot[iCent][iPt]->Scale(1./scale_ratio);
	    hmV0DecLen_rot[iCent][iPt]->Scale(1./scale_ratio);
	    hmV0InvMass_rot[iCent][iPt]->Scale(1./scale_ratio);

	    hmBachDca_rot[iCent][iPt]->Scale(1./scale_ratio);
	    hmDcav0tobach_rot[iCent][iPt]->Scale(1./scale_ratio);
	    hmXiDca_rot[iCent][iPt]->Scale(1./scale_ratio);
	    hmXiDecLen_rot[iCent][iPt]->Scale(1./scale_ratio);
	    hmDecLenDiff_rot[iCent][iPt]->Scale(1./scale_ratio);
	    hmXiInvMass_rot[iCent][iPt]->Scale(1./scale_ratio);
	    hmXiSinth_rot[iCent][iPt]->Scale(1./scale_ratio);
*/            
            cout<<"****************"<<hmDau1nHits_data[iCent][iPt]->GetEntries()<<"************************"<<hmDau1nHits_rot[iCent][iPt]->GetEntries()<<endl;

            hmDau1nHits_data[iCent][iPt]->Draw();
	    hmDau1nHits_data[iCent][iPt]->Add(hmDau1nHits_rot[iCent][iPt], -1./scale_ratio);

            cout<<"**"<<iCent<<"**"<<iPt<<"************"<<hmDau1nHits_data[iCent][iPt]->GetEntries()<<endl;
            //gStyle->SetOptStat(1);
            //break; 

	    hmDau2nHits_data[iCent][iPt]->Add(hmDau2nHits_rot[iCent][iPt], -1./scale_ratio);
            hmBachnHits_data[iCent][iPt]->Add(hmBachnHits_rot[iCent][iPt], -1./scale_ratio);
	    hmDau1Pt_data[iCent][iPt]->Add(hmDau1Pt_rot[iCent][iPt], -1./scale_ratio);
	    hmDau2Pt_data[iCent][iPt]->Add(hmDau2Pt_rot[iCent][iPt], -1./scale_ratio);
            hmBachPt_data[iCent][iPt]->Add(hmBachPt_rot[iCent][iPt], -1./scale_ratio);

	    hmDau1Dca_data[iCent][iPt]->Add(hmDau1Dca_rot[iCent][iPt], -1./scale_ratio);
	    hmDau2Dca_data[iCent][iPt]->Add(hmDau2Dca_rot[iCent][iPt], -1./scale_ratio);
	    hmDca1to2_data[iCent][iPt]->Add(hmDca1to2_rot[iCent][iPt], -1./scale_ratio);
	    hmV0Dca_data[iCent][iPt]->Add(hmV0Dca_rot[iCent][iPt], -1./scale_ratio);
	    hmV0DecLen_data[iCent][iPt]->Add(hmV0DecLen_rot[iCent][iPt], -1./scale_ratio);
	    hmV0InvMass_data[iCent][iPt]->Add(hmV0InvMass_rot[iCent][iPt], -1./scale_ratio);

	    hmBachDca_data[iCent][iPt]->Add(hmBachDca_rot[iCent][iPt], -1./scale_ratio);
	    hmDcav0tobach_data[iCent][iPt]->Add(hmDcav0tobach_rot[iCent][iPt], -1./scale_ratio);
	    hmXiDca_data[iCent][iPt]->Add(hmXiDca_rot[iCent][iPt], -1./scale_ratio);
	    hmXiDecLen_data[iCent][iPt]->Add(hmXiDecLen_rot[iCent][iPt], -1./scale_ratio);
	    hmDecLenDiff_data[iCent][iPt]->Add(hmDecLenDiff_rot[iCent][iPt], -1./scale_ratio);
	    hmXiInvMass_data[iCent][iPt]->Add(hmXiInvMass_rot[iCent][iPt], -1./scale_ratio);
	    hmXiSinth_data[iCent][iPt]->Add(hmXiSinth_rot[iCent][iPt], -1./scale_ratio);

            //==== Get MC Data Hist. ====
	    hmDau1nHits_mc[iCent][iPt] = (TH1F*)mc_file->Get(hDau1nHitsName->Data());
	    hmDau2nHits_mc[iCent][iPt] = (TH1F*)mc_file->Get(hDau2nHitsName->Data()); 
            hmBachnHits_mc[iCent][iPt] = (TH1F*)mc_file->Get(hBachnHitsName->Data());
	    hmDau1Pt_mc[iCent][iPt] = (TH1F*)mc_file->Get(hDau1PtName->Data()); 
	    hmDau2Pt_mc[iCent][iPt] = (TH1F*)mc_file->Get(hDau2PtName->Data()); 
            hmBachPt_mc[iCent][iPt] = (TH1F*)mc_file->Get(hBachPtName->Data()); 

	    hmDau1Dca_mc[iCent][iPt] = (TH1F*)mc_file->Get(hDau1DcaName->Data());
	    hmDau2Dca_mc[iCent][iPt] = (TH1F*)mc_file->Get(hDau2DcaName->Data());
	    hmDca1to2_mc[iCent][iPt] = (TH1F*)mc_file->Get(hDca1to2Name->Data());
	    hmV0Dca_mc[iCent][iPt] = (TH1F*)mc_file->Get(hV0DcaName->Data());
	    hmV0DecLen_mc[iCent][iPt] = (TH1F*)mc_file->Get(hV0DecLenName->Data()); 
	    hmV0InvMass_mc[iCent][iPt] = (TH1F*)mc_file->Get(hV0InvMassName->Data());

	    hmBachDca_mc[iCent][iPt] = (TH1F*)mc_file->Get(hBachDcaName->Data());
	    hmDcav0tobach_mc[iCent][iPt] = (TH1F*)mc_file->Get(hDcav0tobachName->Data());
	    hmXiDca_mc[iCent][iPt] = (TH1F*)mc_file->Get(hXiDcaName->Data());
	    hmXiDecLen_mc[iCent][iPt] = (TH1F*)mc_file->Get(hXiDecLenName->Data()); 
	    hmDecLenDiff_mc[iCent][iPt] = (TH1F*)mc_file->Get(hDecLenDiffName->Data());
	    hmXiInvMass_mc[iCent][iPt] = (TH1F*)mc_file->Get(hXiInvMassName->Data());
	    hmXiSinth_mc[iCent][iPt] = (TH1F*)mc_file->Get(hXiSinthName->Data());
            
	    hmDau1nHits_mc[iCent][iPt]->Sumw2();
	    hmDau2nHits_mc[iCent][iPt]->Sumw2();
            hmBachnHits_mc[iCent][iPt]->Sumw2();
	    hmDau1Pt_mc[iCent][iPt]->Sumw2();
	    hmDau2Pt_mc[iCent][iPt]->Sumw2();
            hmBachPt_mc[iCent][iPt]->Sumw2();

	    hmDau1Dca_mc[iCent][iPt]->Sumw2();
	    hmDau2Dca_mc[iCent][iPt]->Sumw2();
	    hmDca1to2_mc[iCent][iPt]->Sumw2();
	    hmV0Dca_mc[iCent][iPt]->Sumw2();
	    hmV0DecLen_mc[iCent][iPt]->Sumw2();
	    hmV0InvMass_mc[iCent][iPt]->Sumw2();

	    hmBachDca_mc[iCent][iPt]->Sumw2();
	    hmDcav0tobach_mc[iCent][iPt]->Sumw2();
	    hmXiDca_mc[iCent][iPt]->Sumw2();
	    hmXiDecLen_mc[iCent][iPt]->Sumw2();
	    hmDecLenDiff_mc[iCent][iPt]->Sumw2();
	    hmXiInvMass_mc[iCent][iPt]->Sumw2();
	    hmXiSinth_mc[iCent][iPt]->Sumw2();
            
            //std::pair <TH1F*, TH1F*> pair_datamc(hm);
            //==== Plot Comparison ==== 
            Float_t ratio_datatomc = 0;
            Float_t canvas_w = 800; Float_t canvas_l = 600;

            TCanvas* c_dau1nhits = new TCanvas("c1", "c1", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_dau1nhits->Print("dau1nhits_mcdata.pdf[");
            hmDau1nHits_data[iCent][iPt]->SetLineColor(1); 
            hmDau1nHits_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmDau1nHits_data[iCent][iPt]->Integral()/hmDau1nHits_mc[iCent][iPt]->Integral();
            //cout<<hmDau1nHits_data[iCent][iPt]->GetEntries()<<endl;
            //cout<<hmDau1nHits_mc[iCent][iPt]->GetEntries()<<endl;
            //if((hmDau1nHits_data[iCent][iPt]->GetEntries() - 1./scale_ratio*hmDau1nHits_rot[iCent][iPt]->GetEntries()) < 0.) break;
            hmDau1nHits_mc[iCent][iPt]->Scale(ratio_datatomc);
            gPad->SetTicks(1, 1);
            hmDau1nHits_data[iCent][iPt]->Draw("E");
            hmDau1nHits_mc[iCent][iPt]->Draw("E same");
            c_dau1nhits->Print("dau1nhits_mcdata.pdf");

            TCanvas* c_dau2nhits = new TCanvas("c2", "c2", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_dau2nhits->Print("dau2nhits_mcdata.pdf[");
            hmDau2nHits_data[iCent][iPt]->SetLineColor(1); 
            hmDau2nHits_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmDau2nHits_data[iCent][iPt]->Integral()/hmDau2nHits_mc[iCent][iPt]->Integral();
            hmDau2nHits_mc[iCent][iPt]->Scale(ratio_datatomc);
            gPad->SetTicks(1, 1);
            hmDau2nHits_data[iCent][iPt]->Draw("E");
            hmDau2nHits_mc[iCent][iPt]->Draw("E same");
            c_dau2nhits->Print("dau2nhits_mcdata.pdf");

            TCanvas* c_bachnhits = new TCanvas("c3", "c3", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_bachnhits->Print("bachnhits_mcdata.pdf[");
            hmBachnHits_data[iCent][iPt]->SetLineColor(1); 
            hmBachnHits_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmBachnHits_data[iCent][iPt]->Integral()/hmBachnHits_mc[iCent][iPt]->Integral();
            hmBachnHits_mc[iCent][iPt]->Scale(ratio_datatomc);
            gPad->SetTicks(1, 1);
            hmBachnHits_data[iCent][iPt]->Draw("E");
            hmBachnHits_mc[iCent][iPt]->Draw("E same");
            c_bachnhits->Print("bachnhits_mcdata.pdf");

            TCanvas* c_dau1dca = new TCanvas("c4", "c4", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_dau1dca->Print("dau1dca_mcdata.pdf[");
            hmDau1Dca_data[iCent][iPt]->SetLineColor(1); 
            hmDau1Dca_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmDau1Dca_data[iCent][iPt]->Integral()/hmDau1Dca_mc[iCent][iPt]->Integral();
            hmDau1Dca_mc[iCent][iPt]->Scale(ratio_datatomc);
            gPad->SetTicks(1, 1);
            hmDau1Dca_data[iCent][iPt]->Draw("E");
            hmDau1Dca_mc[iCent][iPt]->Draw("E same");
            c_dau1dca->Print("dau1dca_mcdata.pdf");

            TCanvas* c_dau2dca = new TCanvas("c5", "c5", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_dau2dca->Print("dau2dca_mcdata.pdf[");
            hmDau2Dca_data[iCent][iPt]->SetLineColor(1); 
            hmDau2Dca_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmDau2Dca_data[iCent][iPt]->Integral()/hmDau2Dca_mc[iCent][iPt]->Integral();
            hmDau2Dca_mc[iCent][iPt]->Scale(ratio_datatomc);
            gPad->SetTicks(1, 1);
            hmDau2Dca_data[iCent][iPt]->Draw("E");
            hmDau2Dca_mc[iCent][iPt]->Draw("E same");
            cout<<ratio_datatomc<<" ratio of data to mc" <<endl;
            c_dau2dca->Print("dau2dca_mcdata.pdf");

            TCanvas* c_dca1to2 = new TCanvas("c6", "c6", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_dca1to2->Print("dca1to2_mcdata.pdf[");
            hmDca1to2_data[iCent][iPt]->SetLineColor(1); 
            hmDca1to2_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmDca1to2_data[iCent][iPt]->Integral()/hmDca1to2_mc[iCent][iPt]->Integral();
            hmDca1to2_mc[iCent][iPt]->Scale(ratio_datatomc);
            gPad->SetTicks(1, 1);
            hmDca1to2_data[iCent][iPt]->Draw("E");
            hmDca1to2_mc[iCent][iPt]->Draw("E same");
            cout<<ratio_datatomc<<" ratio of data to mc" <<endl;
            c_dca1to2->Print("dca1to2_mcdata.pdf");

            TCanvas* c_v0dca = new TCanvas("c7", "c7", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_v0dca->Print("v0dca_mcdata.pdf[");
            hmV0Dca_data[iCent][iPt]->SetLineColor(1); 
            hmV0Dca_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmV0Dca_data[iCent][iPt]->Integral()/hmV0Dca_mc[iCent][iPt]->Integral();
            hmV0Dca_mc[iCent][iPt]->Scale(ratio_datatomc);
            gPad->SetTicks(1, 1);
            hmV0Dca_data[iCent][iPt]->Draw("E");
            hmV0Dca_mc[iCent][iPt]->Draw("E same");
            cout<<ratio_datatomc<<" ratio of data to mc" <<endl;
            c_v0dca->Print("v0dca_mcdata.pdf");

            TCanvas* c_v0declen = new TCanvas("c8", "c8", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_v0declen->Print("v0declen_mcdata.pdf[");
            hmV0DecLen_data[iCent][iPt]->SetLineColor(1); 
            hmV0DecLen_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmV0DecLen_data[iCent][iPt]->Integral()/hmV0DecLen_mc[iCent][iPt]->Integral();
            hmV0DecLen_mc[iCent][iPt]->Scale(ratio_datatomc);
            gPad->SetTicks(1, 1);
            hmV0DecLen_data[iCent][iPt]->Draw("E");
            hmV0DecLen_mc[iCent][iPt]->Draw("E same");
            cout<<ratio_datatomc<<" ratio of data to mc"<<endl;
            c_v0declen->Print("v0declen_mcdata.pdf");

            TCanvas* c_v0mass = new TCanvas("c9", "c9", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_v0mass->Print("v0mass_mcdata.pdf[");
            hmV0InvMass_data[iCent][iPt]->SetLineColor(1); 
            hmV0InvMass_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmV0InvMass_data[iCent][iPt]->Integral()/hmV0InvMass_mc[iCent][iPt]->Integral();
            hmV0InvMass_mc[iCent][iPt]->Scale(ratio_datatomc);
            gPad->SetTicks(1, 1);
            hmV0InvMass_data[iCent][iPt]->Draw("E");
            hmV0InvMass_mc[iCent][iPt]->Draw("E same");
            cout<<ratio_datatomc<<" ratio of data to mc" <<endl;
            c_v0mass->Print("v0mass_mcdata.pdf");

            TCanvas* c_bachdca = new TCanvas("c10", "c10", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_bachdca->Print("bachdca_mcdata.pdf[");
            hmBachDca_data[iCent][iPt]->SetLineColor(1); 
            hmBachDca_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmBachDca_data[iCent][iPt]->Integral()/hmBachDca_mc[iCent][iPt]->Integral();
            hmBachDca_mc[iCent][iPt]->Scale(ratio_datatomc);
            gPad->SetTicks(1, 1);
            hmBachDca_data[iCent][iPt]->Draw("E");
            hmBachDca_mc[iCent][iPt]->Draw("E same");
            cout<<ratio_datatomc<<" ratio of data to mc" <<endl;
            c_bachdca->Print("bachdca_mcdata.pdf");

            TCanvas* c_dcav0tobach = new TCanvas("c11", "c11", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_dcav0tobach->Print("dcav0tobach_mcdata.pdf[");
            hmDcav0tobach_data[iCent][iPt]->SetLineColor(1); 
            hmDcav0tobach_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmDcav0tobach_data[iCent][iPt]->Integral()/hmDcav0tobach_mc[iCent][iPt]->Integral();
            hmDcav0tobach_mc[iCent][iPt]->Scale(ratio_datatomc);
            gPad->SetTicks(1, 1);
            hmDcav0tobach_data[iCent][iPt]->Draw("E");
            hmDcav0tobach_mc[iCent][iPt]->Draw("E same");
            cout<<ratio_datatomc<<" ratio of data to mc" <<endl;
            c_dcav0tobach->Print("dcav0tobach_mcdata.pdf");

            TCanvas* c_xidca = new TCanvas("c12", "c12", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_xidca->Print("xidca_mcdata.pdf[");
            hmXiDca_data[iCent][iPt]->SetLineColor(1); 
            hmXiDca_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmXiDca_data[iCent][iPt]->Integral()/hmXiDca_mc[iCent][iPt]->Integral();
            hmXiDca_mc[iCent][iPt]->Scale(ratio_datatomc);
            gPad->SetTicks(1, 1);
            hmXiDca_data[iCent][iPt]->Draw("E");
            hmXiDca_mc[iCent][iPt]->Draw("E same");
            cout<<ratio_datatomc<<" ratio of data to mc" <<endl;
            c_xidca->Print("xidca_mcdata.pdf");

            TCanvas* c_xideclen = new TCanvas("c13", "c13", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_xideclen->Print("xideclen_mcdata.pdf[");
            hmXiDecLen_data[iCent][iPt]->SetLineColor(1); 
            hmXiDecLen_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmXiDecLen_data[iCent][iPt]->Integral()/hmXiDecLen_mc[iCent][iPt]->Integral();
            hmXiDecLen_mc[iCent][iPt]->Scale(ratio_datatomc);
            gPad->SetTicks(1, 1);
            hmXiDecLen_data[iCent][iPt]->Draw("E");
            hmXiDecLen_mc[iCent][iPt]->Draw("E same");
            cout<<ratio_datatomc<<" ratio of data to mc" <<endl;
            c_xideclen->Print("xideclen_mcdata.pdf");

            TCanvas* c_ximass = new TCanvas("c14", "c14", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_ximass->Print("ximass_mcdata.pdf[");
            hmXiInvMass_data[iCent][iPt]->SetLineColor(1); 
            hmXiInvMass_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmXiInvMass_data[iCent][iPt]->Integral()/hmXiInvMass_mc[iCent][iPt]->Integral();
            hmXiInvMass_mc[iCent][iPt]->Scale(ratio_datatomc);
            gPad->SetTicks(1, 1);
            hmXiInvMass_data[iCent][iPt]->Draw("E");
            hmXiInvMass_mc[iCent][iPt]->Draw("E same");
            cout<<ratio_datatomc<<" ratio of data to mc" <<endl;
            c_ximass->Print("ximass_mcdata.pdf");

	}
    }
    c_dau1nhits->Print("dau1nhits_mcdata.pdf]");
    c_dau2nhits->Print("dau2nhits_mcdata.pdf]");
    c_bachnhits->Print("bachnhits_mcdata.pdf]");
    c_dau1dca->Print("dau1dca_mcdata.pdf]");
    c_dau2dca->Print("dau2dca_mcdata.pdf]");
    c_dca1to2->Print("dca1to2_mcdata.pdf]");
    c_v0dca->Print("v0dca_mcdata.pdf]");
    c_v0declen->Print("v0declen_mcdata.pdf]");
    c_v0mass->Print("v0mass_mcdata.pdf]");
    c_bachdca->Print("bachdca_mcdata.pdf]");
    c_dcav0tobach->Print("dcav0tobach_mcdata.pdf]");
    c_xidca->Print("xidca_mcdata.pdf]");
    c_xideclen->Print("xideclen_mcdata.pdf]");
    c_ximass->Print("ximass_mcdata.pdf]");
}
