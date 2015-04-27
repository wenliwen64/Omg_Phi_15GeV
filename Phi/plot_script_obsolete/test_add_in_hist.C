void test_add_in_hist(){
    gRandom -> SetSeed(0);
    TFile* f1 = new TFile("test_add.root", "recreate");
    TH1F* h1 = new TH1F("h1", "h1", 100, -1, 1);  
    TH1F* h2 = new TH1F("h2", "h2", 100, -1, 1);  
    
    float a = 0; 
    float b = 0; 
    for(int i = 0; i< 10000; i++){
        gRandom -> Rannor(a, b);
        h1 -> Fill(a);
        h2 -> Fill(b);
    }
    f1 -> Write();
    f1 -> Close();
    
    //h1 -> Draw();
    //delete h2;
    //delete h1;

    TFile* f2 = new TFile("test_add.root", "r");
    TH1F* new_h2 = (TH1F*)f2 -> Get("h2");
    TH1F* new_h1 = (TH1F*)f2 -> Get("h1");
    //new_h1 -> Add(new_h2);
    delete new_h2;
   
    //new_h2 -> Draw();
    new_h1 -> Draw();
}
