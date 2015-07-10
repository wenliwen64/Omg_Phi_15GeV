{
    gStyle->SetOptStat(0);
    const Float_t pdgV0Mass = 1.11568;
    const Float_t pdgXiMass = 1.67245;
    const Float_t masswidth = 0.07 ;
    const Int_t kCentBin = 2;
    const Int_t kPtBin = 6;
    const Float_t ptbd[kPtBin+1] = {0.7, 1.2, 1.6, 2.0, 2.4, 2.8, 3.6};
    const Float_t centbd[kCentBin+1] = {1.5, 6.5, 8.5};

    TFile* mc_file = new TFile("mcomg_fp0.coarse_ptbin_test.cuts.histo.root", "read");
    TFile* data_file = new TFile("test_hist_embedding_15GeV.root", "read");
    
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
            ratio_datatomc = hmDau1nHits_data[iCent][iPt]->GetEntries()/hmDau1nHits_mc[iCent][iPt]->GetEntries();
            hmDau1nHits_mc[iCent][iPt]->Scale(ratio_datatomc);
            hmDau1nHits_data[iCent][iPt]->Draw("E");
            hmDau1nHits_mc[iCent][iPt]->Draw("E same");
            c_dau1nhits->Print("dau1nhits_mcdata.pdf");

            TCanvas* c_dau2nhits = new TCanvas("c2", "c2", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_dau2nhits->Print("dau2nhits_mcdata.pdf[");
            hmDau2nHits_data[iCent][iPt]->SetLineColor(1); 
            hmDau2nHits_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmDau2nHits_data[iCent][iPt]->GetEntries()/hmDau2nHits_mc[iCent][iPt]->GetEntries();
            hmDau2nHits_mc[iCent][iPt]->Scale(ratio_datatomc);
            hmDau2nHits_data[iCent][iPt]->Draw("E");
            hmDau2nHits_mc[iCent][iPt]->Draw("E same");
            c_dau2nhits->Print("dau2nhits_mcdata.pdf");

            TCanvas* c_bachnhits = new TCanvas("c3", "c3", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_bachnhits->Print("bachnhits_mcdata.pdf[");
            hmBachnHits_data[iCent][iPt]->SetLineColor(1); 
            hmBachnHits_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmDau1nHits_data[iCent][iPt]->GetEntries()/hmDau1nHits_mc[iCent][iPt]->GetEntries();
            hmBachnHits_mc[iCent][iPt]->Scale(ratio_datatomc);
            hmBachnHits_data[iCent][iPt]->Draw("E");
            hmBachnHits_mc[iCent][iPt]->Draw("E same");
            c_bachnhits->Print("bachnhits_mcdata.pdf");

            TCanvas* c_dau1dca = new TCanvas("c4", "c4", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_dau1dca->Print("dau1dca_mcdata.pdf[");
            hmDau1Dca_data[iCent][iPt]->SetLineColor(1); 
            hmDau1Dca_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmDau1nHits_data[iCent][iPt]->GetEntries()/hmDau1nHits_mc[iCent][iPt]->GetEntries();
            hmDau1Dca_mc[iCent][iPt]->Scale(ratio_datatomc);
            hmDau1Dca_data[iCent][iPt]->Draw("E");
            hmDau1Dca_mc[iCent][iPt]->Draw("E same");
            c_dau1dca->Print("dau1dca_mcdata.pdf");

            TCanvas* c_dau2dca = new TCanvas("c5", "c5", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_dau2dca->Print("dau2dca_mcdata.pdf[");
            hmDau2Dca_data[iCent][iPt]->SetLineColor(1); 
            hmDau2Dca_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmDau2Dca_data[iCent][iPt]->GetEntries()/hmDau2Dca_mc[iCent][iPt]->GetEntries();
            hmDau2Dca_mc[iCent][iPt]->Scale(ratio_datatomc);
            hmDau2Dca_data[iCent][iPt]->Draw("E");
            hmDau2Dca_mc[iCent][iPt]->Draw("E same");
            cout<<ratio_datatomc<<" ratio of data to mc" <<endl;
            c_dau2dca->Print("dau2dca_mcdata.pdf");

            TCanvas* c_dca1to2 = new TCanvas("c6", "c6", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_dca1to2->Print("dca1to2_mcdata.pdf[");
            hmDca1to2_data[iCent][iPt]->SetLineColor(1); 
            hmDca1to2_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmDca1to2_data[iCent][iPt]->GetEntries()/hmDca1to2_mc[iCent][iPt]->GetEntries();
            hmDca1to2_mc[iCent][iPt]->Scale(ratio_datatomc);
            hmDca1to2_data[iCent][iPt]->Draw("E");
            hmDca1to2_mc[iCent][iPt]->Draw("E same");
            cout<<ratio_datatomc<<" ratio of data to mc" <<endl;
            c_dca1to2->Print("dca1to2_mcdata.pdf");

            TCanvas* c_v0dca = new TCanvas("c7", "c7", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_v0dca->Print("v0dca_mcdata.pdf[");
            hmV0Dca_data[iCent][iPt]->SetLineColor(1); 
            hmV0Dca_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmV0Dca_data[iCent][iPt]->GetEntries()/hmV0Dca_mc[iCent][iPt]->GetEntries();
            hmV0Dca_mc[iCent][iPt]->Scale(ratio_datatomc);
            hmV0Dca_data[iCent][iPt]->Draw("E");
            hmV0Dca_mc[iCent][iPt]->Draw("E same");
            cout<<ratio_datatomc<<" ratio of data to mc" <<endl;
            c_v0dca->Print("v0dca_mcdata.pdf");

            TCanvas* c_v0declen = new TCanvas("c8", "c8", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_v0declen->Print("v0declen_mcdata.pdf[");
            hmV0DecLen_data[iCent][iPt]->SetLineColor(1); 
            hmV0DecLen_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmV0DecLen_data[iCent][iPt]->GetEntries()/hmV0DecLen_mc[iCent][iPt]->GetEntries();
            hmV0DecLen_mc[iCent][iPt]->Scale(ratio_datatomc);
            hmV0DecLen_data[iCent][iPt]->Draw("E");
            hmV0DecLen_mc[iCent][iPt]->Draw("E same");
            cout<<ratio_datatomc<<" ratio of data to mc"<<endl;
            c_v0declen->Print("v0declen_mcdata.pdf");

            TCanvas* c_v0mass = new TCanvas("c9", "c9", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_v0mass->Print("v0mass_mcdata.pdf[");
            hmV0InvMass_data[iCent][iPt]->SetLineColor(1); 
            hmV0InvMass_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmV0InvMass_data[iCent][iPt]->GetEntries()/hmV0InvMass_mc[iCent][iPt]->GetEntries();
            hmV0InvMass_mc[iCent][iPt]->Scale(ratio_datatomc);
            hmV0InvMass_data[iCent][iPt]->Draw("E");
            hmV0InvMass_mc[iCent][iPt]->Draw("E same");
            cout<<ratio_datatomc<<" ratio of data to mc" <<endl;
            c_v0mass->Print("v0mass_mcdata.pdf");

            TCanvas* c_bachdca = new TCanvas("c10", "c10", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_bachdca->Print("bachdca_mcdata.pdf[");
            hmBachDca_data[iCent][iPt]->SetLineColor(1); 
            hmBachDca_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmBachDca_data[iCent][iPt]->GetEntries()/hmBachDca_mc[iCent][iPt]->GetEntries();
            hmBachDca_mc[iCent][iPt]->Scale(ratio_datatomc);
            hmBachDca_data[iCent][iPt]->Draw("E");
            hmBachDca_mc[iCent][iPt]->Draw("E same");
            cout<<ratio_datatomc<<" ratio of data to mc" <<endl;
            c_bachdca->Print("bachdca_mcdata.pdf");

            TCanvas* c_dcav0tobach = new TCanvas("c11", "c11", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_dcav0tobach->Print("dcav0tobach_mcdata.pdf[");
            hmDcav0tobach_data[iCent][iPt]->SetLineColor(1); 
            hmDcav0tobach_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmDcav0tobach_data[iCent][iPt]->GetEntries()/hmDcav0tobach_mc[iCent][iPt]->GetEntries();
            hmDcav0tobach_mc[iCent][iPt]->Scale(ratio_datatomc);
            hmDcav0tobach_data[iCent][iPt]->Draw("E");
            hmDcav0tobach_mc[iCent][iPt]->Draw("E same");
            cout<<ratio_datatomc<<" ratio of data to mc" <<endl;
            c_dcav0tobach->Print("dcav0tobach_mcdata.pdf");

            TCanvas* c_xidca = new TCanvas("c12", "c12", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_xidca->Print("xidca_mcdata.pdf[");
            hmXiDca_data[iCent][iPt]->SetLineColor(1); 
            hmXiDca_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmXiDca_data[iCent][iPt]->GetEntries()/hmXiDca_mc[iCent][iPt]->GetEntries();
            hmXiDca_mc[iCent][iPt]->Scale(ratio_datatomc);
            hmXiDca_data[iCent][iPt]->Draw("E");
            hmXiDca_mc[iCent][iPt]->Draw("E same");
            cout<<ratio_datatomc<<" ratio of data to mc" <<endl;
            c_xidca->Print("xidca_mcdata.pdf");

            TCanvas* c_xideclen = new TCanvas("c13", "c13", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_xideclen->Print("xideclen_mcdata.pdf[");
            hmXiDecLen_data[iCent][iPt]->SetLineColor(1); 
            hmXiDecLen_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmXiDecLen_data[iCent][iPt]->GetEntries()/hmXiDecLen_mc[iCent][iPt]->GetEntries();
            hmXiDecLen_mc[iCent][iPt]->Scale(ratio_datatomc);
            hmXiDecLen_data[iCent][iPt]->Draw("E");
            hmXiDecLen_mc[iCent][iPt]->Draw("E same");
            cout<<ratio_datatomc<<" ratio of data to mc" <<endl;
            c_xideclen->Print("xideclen_mcdata.pdf");

            TCanvas* c_ximass = new TCanvas("c14", "c14", canvas_w, canvas_l); 
            if(iCent == 0 && iPt == 0) c_ximass->Print("ximass_mcdata.pdf[");
            hmXiInvMass_data[iCent][iPt]->SetLineColor(1); 
            hmXiInvMass_mc[iCent][iPt]->SetLineColor(2);
            ratio_datatomc = hmXiInvMass_data[iCent][iPt]->GetEntries()/hmXiInvMass_mc[iCent][iPt]->GetEntries();
            hmXiInvMass_mc[iCent][iPt]->Scale(ratio_datatomc);
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
