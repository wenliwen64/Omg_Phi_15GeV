void run_spectra(){
    gSystem->Load("./StPhiDownMaker_cc.so");
    StPhiDownMaker* aMaker = new StPhiDownMaker("Phi");
    //aMaker->Init("./0818_overview.reweight.histo.root");
//    aMaker->Init("./0818_overview_11GeV.reweight.histo.root", "./0826_11GeV_omg.local_analysis.root", "./0826_11GeV_omgrot.local_analysis.root");// feng's upstream files
    //std::string overviewfile = "./0820_15GeV_overview.histo.root";
    //std::string overviewfile = "./0826_15GeV_overview.reweight.histo.root";
    //std::string datfile = "./0826_omg_15GeV.local_analysis.root";
    //std::string rotfile = "./0826_omgrot_15GeV.local_analysis.root";
    std::string datfile = "./phi_15GeV.histo.root";
    aMaker->Init(datfile);// feng's upstream files
    aMaker->AnalyzeBES();
}
