Float_t findEff(TH1F* eff_hist, Float_t pt);
void correct_eff(){
    const Int_t kPtBin = 6;
    const Float_t pdgmass_xi = 1.67245;
    TH1F* eff_hist[2];
    const Float_t ptbd[kPtBin+1] = {1.0, 1.2, 1.6, 2.0, 2.4, 2.8, 3.6};
    TH1F* hmptbin =new TH1F("hmptbin", "pt bin finder", kPtBin, ptbd);
    TFile* fermidst_file = new TFile("0628_2015_omg.local_analysis.root", "read");
    TTree* fermidst = (TTree*)fermidst_file->Get("Xi_FermiDst");
    TFile* eff_file = new TFile("mcomg_fp0.cuts.histo.root ", "read");
    eff_hist[0] = (TH1F*)eff_file->Get("hrcxiptcen0");
    eff_hist[1] = (TH1F*)eff_file->Get("hrcxiptcen1");

    TFile* file = new TFile("corrected_eff.root", "recreate");
    file->cd();
    TH1F* h_010[kPtBin];
    TH1F* h_1060[kPtBin];
    for(int iptbin; iptbin < kPtBin; iptbin++){
        TString hname_010;
        TString hname_1060;
        hname_010.Form("sig_xipt%dcent_010", iptbin+1);
        hname_1060.Form("sig_xipt%dcent_1060", iptbin+1);
        h_010[iptbin] = new TH1F(hname_010.Data(), hname_010.Data(), 100, pdgmass_xi - 0.07, pdgmass_xi + 0.07);
        h_1060[iptbin] = new TH1F(hname_1060.Data(), hname_1060.Data(), 100, pdgmass_xi - 0.07, pdgmass_xi + 0.07);
        h_010[iptbin]->Sumw2();
        h_1060[iptbin]->Sumw2();
    }

    Int_t nentries = fermidst->GetEntries();
    for(int ientry = 0; ientry < nentries; ientry++){
	if(ientry % 10000 == 0) cout<<ientry<< "events processed/"<<nentries<<endl;
        fermidst -> GetEntry(ientry);
        TLeaf* leaf_centbin9 = fermidst->GetLeaf("centbin9");
        TLeaf* leaf_nxi = fermidst->GetLeaf("nxi");
        TLeaf* leaf_ximass = fermidst->GetLeaf("XiMass");
        TLeaf* leaf_xipt = fermidst->GetLeaf("XiPt");
        
        Int_t nxi = leaf_nxi->GetValue(0);
	Int_t centbin = leaf_centbin9->GetValue(0);
        //cout<<nxi<<"xi paricles"<<endl;
        for(int ixi = 0; ixi < nxi; ixi++){
            Float_t ximass = leaf_ximass->GetValue(ixi);
            Float_t xipt = leaf_xipt->GetValue(ixi);
    
            Int_t ptbin = hmptbin->FindBin(xipt);
            if(ptbin >= kPtBin) ptbin = -1;
            if(ptbin <= 0 || centbin < 0) continue;
            if(centbin == 7 || centbin == 8){//h_010
                icent = 0; 
		Float_t eff = findEff(eff_hist[icent], xipt);
                if(eff <= 0) continue;
                h_010[ptbin-1]->Fill(ximass, 1/eff);
	    }
            else if(centbin == 2 || centbin == 3 || centbin== 4 || centbin== 5 || centbin == 6){
                icent = 1;
		Float_t eff = findEff(eff_hist[icent], xipt);
                if(eff <= 0) continue;
                h_1060[ptbin-1]->Fill(ximass, 1/eff);
	    }
            else{}
	} 
    } 
    file->Write();
}
 
Float_t findEff(TH1F* eff_hist, Float_t pt){
    Int_t bin = eff_hist->FindBin(pt);
    Float_t eff = eff_hist->GetBinContent(bin);
    return eff;
}
