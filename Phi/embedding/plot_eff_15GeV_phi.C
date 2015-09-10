{
   std::ifstream infile_flat("weight_phi_fp.txt");
   std::ifstream infile_exp("weight_phi_exp.txt");

   const Int_t kCentBin = 9;
   const Int_t kPtBin = 11;
   Float_t  PtBinCenter[kPtBin]={0.45,0.55,0.65,0.75,0.9,1.15,1.5,1.85,2.25,3.0,4.25};

   // Read Eff data
   Int_t centbin = 0;
   Int_t ptbin = 0;
   Float_t dummy = 0;
   Float_t eff = 0;
   Float_t efferr = 0;
   Float_t eff_flat[kCentBin][kPtBin]; 
   Float_t efferr_flat[kCentBin][kPtBin]; 
   Float_t eff_exp[kCentBin][kPtBin]; 
   Float_t efferr_exp[kCentBin][kPtBin]; 

   while(infile_flat >> centbin){
       infile_flat >> ptbin >> dummy >> dummy >> dummy >> dummy >> eff >> efferr;         
       eff_flat[centbin][ptbin] = eff*pow(10, -(8-centbin));
       efferr_flat[centbin][ptbin] = efferr*pow(10, -(8-centbin));
   }

   while(infile_exp >> centbin){
       infile_exp >> ptbin >> dummy >> dummy >> dummy >> dummy >> eff >> efferr;
       eff_exp[centbin][ptbin] = eff*pow(10, -(8-centbin));
       efferr_exp[centbin][ptbin] = efferr*pow(10, -(8-centbin));
   }

   TCanvas* can_eff = new TCanvas();
   gPad->SetLogy();
   gPad->SetTicks(1, 1);

   TGraphErrors* gr_flat = NULL;
   TGraphErrors* gr_exp = NULL;
   TLegend* leg = new TLegend(0.5, 0.2, 0.7, 0.3);
   leg->SetBorderSize(0);
   for(int i = 0; i < kCentBin; i++){
       gr_flat = new TGraphErrors(11, PtBinCenter, eff_flat[i], 0, efferr_flat[i]);
       gr_flat->SetMarkerStyle(20);
       gr_flat->SetMarkerSize(1.3);
       gr_flat->SetLineColor(1);
       gr_flat->SetMarkerColor(1);
       if(i == 0){
           gr_flat->SetTitle("");
           gr_flat->SetMaximum(10.);
           gr_flat->SetMinimum(10e-14);
           gr_flat->GetXaxis()->SetTitle("Pt(GeV/c)");
           gr_flat->GetYaxis()->SetTitle("efficiency");
	   gr_flat->Draw("AP");
           leg->AddEntry(gr_flat, "flat", "p");
       }
       else
           gr_flat->Draw("P");

       gr_exp = new TGraphErrors(11, PtBinCenter, eff_exp[i], 0, efferr_exp[i]);
       gr_exp->SetMarkerStyle(20);
       gr_exp->SetMarkerSize(1.3);
       gr_exp->SetLineColor(2);
       gr_exp->SetMarkerColor(2);
       gr_exp->Draw("P");
       if(i == 0)
           leg->AddEntry(gr_exp, "exp", "p");
   }
   leg->Draw("same");
   gPad->SaveAs("15GeV_phieff_exp_flat.pdf");   
}
