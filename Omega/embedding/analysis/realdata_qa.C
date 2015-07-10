void realdata_qa(){
    TFile* mc_file = new TFile("mcomg_fp0.coarse_ptbin_test.cuts.histo.root", "read");
    TFile* real_file = new TFile("0709_2015_omg.local_analysis.root", "read");
    TFile* histoutput_file = new TFile("test_hist_embedding_15GeV.root", "recreate");

    //==== Book Histograms ====
    const Float_t pdgV0Mass = 1.11568;
    const Float_t pdgXiMass = 1.67245;
    const Float_t masswidth = 0.07 ;
    const Int_t kCentBin = 2;
    const Int_t kPtBin = 6;
    const Float_t ptbd[kPtBin+1] = {0.7, 1.2, 1.6, 2.0, 2.4, 2.8, 3.6};
    const Float_t centbd[kCentBin+1] = {1.5, 6.5, 8.5};
    TH1F* hmPtBinFinder = new TH1F("hmPtBinFinder", "hmPtBinFinder", kPtBin, ptbd);
    TH1F* hmCentBinFinder = new TH1F("hmCentBinFinder", "hmCentBinFinder", kCentBin, centbd);

    //Declare histograms
    histoutput_file->cd();
    TH1F* hmDau1nHits[kCentBin][kPtBin];
    TH1F* hmDau2nHits[kCentBin][kPtBin];
    TH1F* hmBachnHits[kCentBin][kPtBin];
    TH1F* hmDau1Pt[kCentBin][kPtBin];
    TH1F* hmDau2Pt[kCentBin][kPtBin];
    TH1F* hmBachPt[kCentBin][kPtBin];

    TH1F* hmDau1Dca[kCentBin][kPtBin];
    TH1F* hmDau2Dca[kCentBin][kPtBin];
    TH1F* hmDca1to2[kCentBin][kPtBin];
    TH1F* hmV0Dca[kCentBin][kPtBin]; 
    TH1F* hmV0DecLen[kCentBin][kPtBin];
    TH1F* hmV0InvMass[kCentBin][kPtBin];
    //TH1F* V0Colinear[kCentBin][kPtBin]   new TH1F(hV0ColinearName.Data(), hV0ColinearTitle.Data(), 100, 0, 1);

    cout<<"happy start"<<kPtBin<<endl;
    TH1F* hmBachDca[kCentBin][kPtBin];
    TH1F* hmDcav0tobach[kCentBin][kPtBin];
    TH1F* hmXiDca[kCentBin][kPtBin];
    TH1F* hmXiDecLen[kCentBin][kPtBin];
    TH1F* hmDecLenDiff[kCentBin][kPtBin];
    TH1F* hmXiInvMass[kCentBin][kPtBin];
    //TH1F* XiColinear[kCentBin][kPtBin]   new TH1F(hXiColinearName.Data(), hXiColinearTitle.Data(), 100, 0, 1);
    TH1F* hmXiSinth[kCentBin][kPtBin];

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

	    TString hDau1nHitsTitle;
	    TString hDau2nHitsTitle;
            TString hBachnHitsTitle;
	    TString hDau1PtTitle;
	    TString hDau2PtTitle;
            TString hBachPtTitle;

	    TString hDau1DcaTitle; 
	    TString hDau2DcaTitle; 
	    TString hDca1to2Title; 
	    TString hV0DcaTitle; 
	    TString hV0DecLenTitle;
	    TString hV0InvMassTitle;
	    //TString hV0ColinearTitle;

	    TString hBachDcaTitle; 
	    TString hDcav0tobachTitle;
	    TString hXiDcaTitle;
	    TString hXiDecLenTitle;
	    TString hDecLenDiffTitle;
	    TString hXiInvMassTitle;
            //TString hXiColinearTitle;
            TString hXiSinthTitle;

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

	    hDau1nHitsTitle.Form("Dau1_nHits cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]);
	    hDau2nHitsTitle.Form("Dau2_nHits cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]);
            hBachnHitsTitle.Form("Bach_nHits cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]);
	    hDau1PtTitle.Form("Dau1_Pt cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]);
	    hDau2PtTitle.Form("Dau2_Pt cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]);
            hBachPtTitle.Form("Bach_Pt cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]);

	    hDau1DcaTitle.Form("Dau1_Dca cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]);
	    hDau2DcaTitle.Form("Dau2_Dca cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]);
	    hDca1to2Title.Form("Dca1to2 cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]);
	    hV0DcaTitle.Form("V0_Dca cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]);
	    hV0DecLenTitle.Form("V0_DecLen cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]);
	    hV0InvMassTitle.Form("V0_Invmass cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]);
	    //hV0ColinearTitle.Form("V0_Colinear cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]);

	    hBachDcaTitle.Form("Bach_Dca cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]);
	    hDcav0tobachTitle.Form("Dca_v0tobach cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]);
	    hXiDcaTitle.Form("Xi_Dca cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]);
	    hXiDecLenTitle.Form("Xi_DecLen cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]);
	    hDecLenDiffTitle.Form("DecLenDiff cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]);
	    hXiInvMassTitle.Form("Xi_InvMass cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]);
	    //hXiColinearTitle.Form("Xi_Colinear cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]); 
	    hXiSinthTitle.Form("Xi_sinth cen%d(%3.1f < pt < %3.1f)", iCent, ptbd[iPt], ptbd[iPt+1]); 

	    hmDau1nHits[iCent][iPt] = new TH1F(hDau1nHitsName.Data(), hDau1nHitsTitle.Data(), 50, -0.5, 49.5);//new
	    hmDau2nHits[iCent][iPt] = new TH1F(hDau2nHitsName.Data(), hDau2nHitsTitle.Data(), 50, -0.5, 49.5);//new
            hmBachnHits[iCent][iPt] = new TH1F(hBachnHitsName.Data(), hBachnHitsTitle.Data(), 50, -0.5, 49.5);
	    hmDau1Pt[iCent][iPt] = new TH1F(hDau1PtName.Data(), hDau1PtTitle.Data(), 100, 0, 10);
	    hmDau2Pt[iCent][iPt] = new TH1F(hDau2PtName.Data(), hDau2PtTitle.Data(), 100, 0, 10);
            hmBachPt[iCent][iPt] = new TH1F(hBachPtName.Data(), hBachPtTitle.Data(), 100, 0, 10);

	    hmDau1Dca[iCent][iPt] = new TH1F(hDau1DcaName.Data(), hDau1DcaTitle.Data(), 100, 0, 10.0);
	    hmDau2Dca[iCent][iPt] = new TH1F(hDau2DcaName.Data(), hDau2DcaTitle.Data(), 100, 0, 10.0);
	    hmDca1to2[iCent][iPt] = new TH1F(hDca1to2Name.Data(), hDca1to2Title.Data(), 200, 0, 1.6);
	    hmV0Dca[iCent][iPt] = new TH1F(hV0DcaName.Data(), hV0DcaTitle.Data(), 200, 0, 5.5);
	    hmV0DecLen[iCent][iPt] = new TH1F(hV0DecLenName.Data(), hV0DecLenTitle.Data(), 200, 0, 100); 
	    hmV0InvMass[iCent][iPt] = new TH1F(hV0InvMassName.Data(), hV0InvMassTitle.Data(), 300, pdgV0Mass-masswidth, pdgV0Mass+masswidth);
	    //hmV0Colinear[iCent][iPt] = new TH1F(hV0ColinearName.Data(), hV0ColinearTitle.Data(), 100, 0, 1);

	    cout<<"happy mid pt = "<<iPt<<endl;
	    hmBachDca[iCent][iPt] = new TH1F(hBachDcaName.Data(), hBachDcaTitle.Data(), 200, 0, 25);
	    hmDcav0tobach[iCent][iPt] = new TH1F(hDcav0tobachName.Data(), hDcav0tobachTitle.Data(), 200, 0, 2);
	    hmXiDca[iCent][iPt] = new TH1F(hXiDcaName.Data(), hXiDcaTitle.Data(), 100, 0, 1.2);
	    hmXiDecLen[iCent][iPt] = new TH1F(hXiDecLenName.Data(), hXiDecLenTitle.Data(), 200, 0, 80);
	    hmDecLenDiff[iCent][iPt] = new TH1F(hDecLenDiffName.Data(), hDecLenDiffTitle.Data(), 200, 0, 80);
	    hmXiInvMass[iCent][iPt] = new TH1F(hXiInvMassName.Data(), hXiInvMassTitle.Data(), 300, pdgXiMass - masswidth, pdgXiMass + masswidth);
	    //hmXiColinear[iCent][iPt] = new TH1F(hXiColinearName.Data(), hXiColinearTitle.Data(), 100, 0, 1);
	    hmXiSinth[iCent][iPt] = new TH1F(hXiSinthName.Data(), hXiSinthTitle.Data(), 200, 0, 0.3);
             
	}
    }
    // process real data 
    TTree* fermidst = (TTree*)real_file->Get("Xi_FermiDst");

    Int_t nentries = fermidst->GetEntries();
    cout<<"in total there are nentries"<<nentries<<endl;
    for(int ientry = 0; ientry < nentries; ientry++){
        if(ientry % 10000 == 0) cout << ientry << " entries processed!/" << nentries << endl;
	fermidst->GetEntry(ientry);

        TLeaf* leaf_nrefmult = fermidst->GetLeaf("nrefmult");
	TLeaf* leaf_cenbin9 = fermidst->GetLeaf("centbin9");
	TLeaf* leaf_nxi = fermidst->GetLeaf("nxi");
	TLeaf* leaf_v0mass = fermidst->GetLeaf("V0Mass");
	TLeaf* leaf_v0declen = fermidst->GetLeaf("V0DecLen");
	TLeaf* leaf_v0dca = fermidst->GetLeaf("V0Dca");
	//TLeaf* leaf_v0x = fermidst->GetLeaf("v0x");
	//TLeaf* leaf_v0y = fermidst->GetLeaf("v0y");
	//TLeaf* leaf_v0z = fermidst->GetLeaf("v0z");
	//TLeaf* leaf_v0px = fermidst->GetLeaf("v0px");
	//TLeaf* leaf_v0py = fermidst->GetLeaf("v0py");
	//TLeaf* leaf_v0pz = fermidst->GetLeaf("v0pz");
	//TLeaf* leaf_dau1id = fermidst->GetLeaf("Dau1id");
	TLeaf* leaf_dau1dca = fermidst->GetLeaf("Dau1Dca");
	TLeaf* leaf_dau1nhits = fermidst->GetLeaf("Dau1nHits");
	TLeaf* leaf_dau1pt = fermidst->GetLeaf("Dau1Pt");
	//TLeaf* leaf_dau2id = fermidst->GetLeaf("Dau2id");
	TLeaf* leaf_dau2dca = fermidst->GetLeaf("Dau2Dca");
	TLeaf* leaf_dau2nhits = fermidst->GetLeaf("Dau2nHits");
	TLeaf* leaf_dau2pt = fermidst->GetLeaf("Dau2Pt");
        TLeaf* leaf_bachpt = fermidst->GetLeaf("BachPt");
	TLeaf* leaf_dca1to2 = fermidst->GetLeaf("Dca1to2");
	TLeaf* leaf_bachnhits = fermidst->GetLeaf("BachnHits");
	//TLeaf* leaf_bachid = fermidst->GetLeaf("bachid");
	TLeaf* leaf_bachdca = fermidst->GetLeaf("BachDca");
	TLeaf* leaf_dcav0tobach = fermidst->GetLeaf("Dcav0tobach");
	TLeaf* leaf_xipt = fermidst -> GetLeaf("XiPt");
	//TLeaf* leaf_xix = fermidst -> GetLeaf("xix");
	//TLeaf* leaf_xiy = fermidst -> GetLeaf("xiy");
	//TLeaf* leaf_xiz = fermidst -> GetLeaf("xiz");
	//TLeaf* leaf_xipx = fermidst -> GetLeaf("xipx");
	//TLeaf* leaf_xipy = fermidst -> GetLeaf("xipy");
	//TLeaf* leaf_xipz = fermidst -> GetLeaf("xipz");
	TLeaf* leaf_xidca = fermidst->GetLeaf("XiDca");
	TLeaf* leaf_xipt = fermidst->GetLeaf("XiPt");
	TLeaf* leaf_xideclen = fermidst->GetLeaf("XiDecLen");
	TLeaf* leaf_xirapidity = fermidst->GetLeaf("XiRapidity");
	//TLeaf* leaf_ximcid = fermidst->GetLeaf("ximcid");
	TLeaf* leaf_ximass = fermidst->GetLeaf("XiMass");
	TLeaf* leaf_xisinth = fermidst->GetLeaf("XiSinth");

	int nxi = leaf_nxi->GetValue(0);
        //cout<<"nxi = "<<nxi<<" "<<leaf_nrefmult->GetValue(0)<<endl;
        int cenbin = leaf_cenbin9->GetValue(0);
        cenbin = hmCentBinFinder->FindBin(cenbin);
        if(cenbin > kCentBin) cenbin = -1;
        if(cenbin <= 0) continue;
        cenbin = cenbin - 1;
        for(int ixi = 0; ixi < nxi; ixi++){
            int dau1nhits = leaf_dau1nhits->GetValue(ixi);	   
            int dau2nhits = leaf_dau2nhits->GetValue(ixi);
            int bachnhits = leaf_bachnhits->GetValue(ixi);
            float dau1pt = leaf_dau1pt->GetValue(ixi);
            float dau2pt = leaf_dau2pt->GetValue(ixi);
            float bachpt = leaf_bachpt->GetValue(ixi);
            float dau1dca = leaf_dau1dca->GetValue(ixi);
            float dau2dca = leaf_dau2dca->GetValue(ixi);
            float dca1to2 = leaf_dca1to2->GetValue(ixi);
            float v0dca = leaf_v0dca->GetValue(ixi);
            float v0declen = leaf_v0declen->GetValue(ixi);
            float v0mass = leaf_v0mass->GetValue(ixi);
            float bachdca = leaf_bachdca->GetValue(ixi);
            float dcav0tobach = leaf_dcav0tobach->GetValue(ixi);
            float xipt = leaf_xipt->GetValue(ixi);
            float xidca = leaf_xidca->GetValue(ixi); 
            float xideclen = leaf_xideclen->GetValue(ixi);
            float xirapidity = leaf_xirapidity->GetValue(ixi);
            float ximass = leaf_ximass->GetValue(ixi);
            float xisinth = leaf_xisinth->GetValue(ixi);

            int ptbin = hmPtBinFinder->FindBin(xipt); 
            if(ptbin > kPtBin) ptbin = -1; 
            if(ptbin <= 0) continue;
            ptbin = ptbin - 1;

            hmDau1nHits[cenbin][ptbin]->Fill(dau1nhits); 
            hmDau2nHits[cenbin][ptbin]->Fill(dau2nhits); 
            hmBachnHits[cenbin][ptbin]->Fill(bachnhits); 
            hmDau1Pt[cenbin][ptbin]->Fill(dau1pt);
            hmDau2Pt[cenbin][ptbin]->Fill(dau2pt);
            hmBachPt[cenbin][ptbin]->Fill(bachpt);
            hmDau1Dca[cenbin][ptbin]->Fill(dau1dca);
            hmDau2Dca[cenbin][ptbin]->Fill(dau2dca);
            hmDca1to2[cenbin][ptbin]->Fill(dca1to2);
            hmV0Dca[cenbin][ptbin]->Fill(v0dca);
            hmV0DecLen[cenbin][ptbin]->Fill(v0declen);
            hmV0InvMass[cenbin][ptbin]->Fill(v0mass);

            hmBachDca[cenbin][ptbin]->Fill(bachdca);
            hmDcav0tobach[cenbin][ptbin]->Fill(dcav0tobach);
            hmXiDca[cenbin][ptbin]->Fill(xidca);
            hmXiDecLen[cenbin][ptbin]->Fill(xideclen);
            hmXiInvMass[cenbin][ptbin]->Fill(ximass);
            //hmXiRapidity[cenbin][ptbin]->Fill(xirapdity); 
            hmXiSinth[cenbin][ptbin]->Fill(xisinth);
	}
    }
    histoutput_file->Write();
}
