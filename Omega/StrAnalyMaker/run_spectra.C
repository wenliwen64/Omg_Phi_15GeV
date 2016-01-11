void run_spectra(){
    //gSystem->Load("./StrAnalyMaker_cc.so");
    gROOT->LoadMacro("./StrAnalyMaker.cc++");
    StrAnalyMaker* aMaker = new StrAnalyMaker("omg");
    //aMaker->Init("./0818_overview.reweight.histo.root");
//    aMaker->Init("./0818_overview_11GeV.reweight.histo.root", "./0826_11GeV_omg.local_analysis.root", "./0826_11GeV_omgrot.local_analysis.root");// feng's upstream files
    //std::string overviewfile = "./0820_15GeV_overview.histo.root";
    //std::string overviewfile = "./0826_15GeV_overview.reweight.histo.root";
    std::string overviewfile = "./0901_15GeV_overview.reweight.histo.root";
    //std::string datfile = "./0826_omg_15GeV.local_analysis.root";
    std::string datfile = "./0901_omg_15GeV.local_analysis.root";
    std::string rotfile = "./0826_omgrot_15GeV.local_analysis.root";
    //std::string rotfile = "./0901_omgrot_15GeV.local_analysis.root";
    std::string rotfile1 = "./1130_omgrot1_15GeV.mid_analysis.root";
    std::string rotfile2 = "./1201_omgrot2_15GeV.mid_analysis.root";
    std::string rotfile3 = "./1202_omgrot3_15GeV.mid_analysis.root";
    std::string rotfile4 = "./1203_omgrot4_15GeV.mid_analysis.root";
    std::string fpefffile = "./mcomg_fp.manyeff.histo.root";
    std::string expefffile = "./mcomg_exp.manyeff.histo.root";
    aMaker->Init(overviewfile, datfile, rotfile, rotfile1, rotfile2, rotfile3, rotfile4, fpefffile, expefffile);// feng's upstream files
    aMaker->AnalyzeII();
    //aMaker->Analyze();
}
