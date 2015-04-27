void test_error_after_scale(){

     gRandom -> SetSeed(0);
     float a = 0;
     float b = 0;
     TH1F* h1 = new TH1F("h1", "h1", 100, -1, 1);
     //h1 -> Sumw2();
     TH1F* h2 = new TH1F("h2", "h2", 100, -1, 1);
     //h2 -> Sumw2();
     for(int i = 0; i < 10000; i++){
         gRandom -> Rannor(a, b);
         h1 -> Fill(a);
         h2 -> Fill(b);
        
     }

     TH1F* h1new = (TH1F*)h1 -> Clone("h1new");
     TH1F* h2new = (TH1F*)h2 -> Clone("h2new");

//================================================================

TCanvas* c1 = new TCanvas("c1");
     c1 -> Divide(2,2);
     c1 -> cd(1);
     h1 -> Draw("E");
     c1 -> cd(2);
     TH1F* h1new1 = (TH1F*)h1 -> Clone();
     h1new1 -> Sumw2();
     h1new1 -> Scale(0.1);
     h1new1 -> Draw("E");
     //h1 -> Scale(10);
     
     //h1 -> Draw("E");

     //h1 -> Add(h2);
     c1 -> cd(3);
     TH1F* h1new2 = (TH1F*) h1 -> Clone();
     //h1new2 -> Add(h2, 1);
     h1new2 -> Scale(0.1);
     h1new2 -> SetLineColor(2);
     h1new2 -> Draw("E same");
     c1 -> cd(4);
     TH1F* h1new3 = (TH1F*) h1 -> Clone();
     h1new3 -> Sumw2();
     h1new3 -> Add(h2, 1);
     h1new3 -> Draw("E"); 
c1 -> SaveAs("h1_test_error.eps");
//===============================================================
/*
TCanvas* c2 = new TCanvas("c2");
     c2 -> Divide(2, 2);
     c2 -> cd(1);
     h1new -> Draw("E");
     c2 -> cd(2);
     //h1new -> Scale(0.5);
     h1new -> Draw("E");
  
     c2 -> cd(3);
     //h1new -> Add(h2new);
     h1new -> Draw("E");
     c2 -> cd(4);
     h2new -> Draw("E");
c2 -> SaveAs("h1new_test_error.eps");
*/
}
