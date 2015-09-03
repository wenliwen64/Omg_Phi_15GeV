void run_rawspectra(){
    gSystem->Load("./StrAnalyMaker.so");
    StrAnalyMaker* aMaker = new StrAnalyMaker("Omega");
    //aMaker->Init("./0818_overview.reweight.histo.root");
//    aMaker->Init("./0818_overview_11GeV.reweight.histo.root", "./0826_11GeV_omg.local_analysis.root", "./0826_11GeV_omgrot.local_analysis.root");// feng's upstream files
    //std::string overviewfile = "./0820_15GeV_overview.histo.root";
    //std::string overviewfile = "./0826_15GeV_overview.reweight.histo.root";
    std::string overviewfile = "./0901_15GeV_overview.reweight.histo.root";
    //std::string datfile = "./0826_omg_15GeV.local_analysis.root";
    std::string datfile = "./0901_omg_15GeV.local_analysis.root";
    //std::string rotfile = "./0826_omgrot_15GeV.local_analysis.root";
    std::string rotfile = "./0901_omgrot_15GeV.local_analysis.root";
    std::string fpefffile = "./mcomg_fp0.manyeff.histo.root";
    std::string expefffile = "./mcomg_exp0.manyeff.histo.root";
    aMaker->Init(overviewfile, datfile, rotfile, fpefffile, expefffile);// feng's upstream files
    aMaker->Analyze();
}
