#include "TH1F.h"
int plot_omg(){
    TFile* output_histo = new TFile("output_histo.root", "recreate");
    double mV0 = 1.67245;
    ifstream fList_sig("./total_omg_sig.list");
    ifstream fList_bg("./total_omg_bg.list");
    ofstream output_file("./output_omg");

    float pt_bin[7] = {0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.6}; 
    int cent_bin[3] = {66, 97, 198}; 
    char linefromlist[250];
    int no_centbin = 4;
    int no_ptbin = 6;
    //int nEvents[4] = {2970564, 3340670, 2136072, 10037274};

//===========Initialize Total Invariant Mass Histogram==============================
    TH1F* hist_total_sig = new TH1F("hist_total_sig", "hist_total_sig", 100, mV0 - 0.07, mV0 + 0.07);
    TH1F* hist_total_bg = new TH1F("hist_total_bg", "hist_total_bg", 100, mV0 - 0.07, mV0 + 0.07);
//===========Initialize TObjArray to Store Histograms=============================
    TObjArray HList(0);
    
    for(int i = 0; i < no_centbin; i++){
	for(int j = 0; j < no_ptbin; j++){
	    char hist_name_sig[100];
	    char hist_des_sig[100];
	    char hist_name_bg[100];
	    char hist_des_bg[100];

	    sprintf(hist_name_sig, "sig_omgpt%dcent%d", j+1, i+1);
	    sprintf(hist_des_sig, "sig_omg(%.1f< pt<%.1f)cent%d", pt_bin[j], pt_bin[j+1], i+1);
	    sprintf(hist_name_bg, "bg_omgpt%dcent%d", j+1, i+1);
	    sprintf(hist_des_bg, "bg_omg(%.1f<pt<%.1f)cent%d", pt_bin[j], pt_bin[j+1], i+1);

	    TH1F* new_histo_sig = new TH1F(hist_name_sig, hist_des_sig, 100, mV0 - 0.07, mV0 + 0.07);
            new_histo_sig -> Sumw2();
	    TH1F* new_histo_bg = new TH1F(hist_name_bg, hist_des_bg, 100, mV0 - 0.07, mV0 + 0.07);
            new_histo_bg -> Sumw2();

	    HList.Add(new_histo_sig);
	    HList.Add(new_histo_bg);
	}
    }

//===========Input files listed in filelist===================
//===========QA Plots========================================
/*
    TH1F* hNRefMult  = new TH1F("RefMult", "Reference Multiplicity", 1000, 0.0, 1000.0 ) ;
    TH1F* hSelectedMult  = new TH1F("SelectedRefMult", "Selected Reference Multiplicity", 1000, 0.0, 1000.0 ) ;
    TH1F* hVz_before = new TH1F("vz_before", "vz_before_cut", 400, -100, 100);
    TH1F* hVz_after = new TH1F("vz_after", "vz_after_cut", 400, -100, 100);
    TH2F* hVr_before = new TH2F("vr_before", "vr_before_cut", 200, -5.0, 5.0, 200, -5.0, 5.0);
    TH2F* hVr_after = new TH2F("vr_after", "vr_after_cut", 200, -5.0, 5.0, 200, -5.0, 5.0);
    TH2F* hDeDx = new TH2F("dedxP", "DeDx_P", 300, -3, 3, 300, 0, 3e-5);
*/
//===================Add up signal invariant mass histograms===============================
    int count_sig = 0;
    while(fList_sig.getline(linefromlist, 250)){
        TFile* file = new TFile(linefromlist, "read");
	count_sig++;
        if(count_sig%100==0) std::cout<< count_sig<<" files processed"<<std::endl;
        //===========Invariant Mass Histograms===================================
	for(int i = 0; i < no_centbin; i++){
	    for(int j = 0; j < no_ptbin; j++){
		char hist_name_sig[100];
		sprintf(hist_name_sig, "sig_omgpt%dcent%d", j+1, i+1);

		char cur_hist_name_sig[100];
		sprintf(cur_hist_name_sig, "sig_omgpt%dcent%d", j+1, i+1);

		TH1F* cur_histo_sig = (TH1F*)file->Get(cur_hist_name_sig);
   

                TH1F* hist_sig = (TH1F*)HList.FindObject(hist_name_sig);
                if(cur_histo_sig) hist_sig -> Add(cur_histo_sig);
                if(cur_histo_sig) hist_total_sig -> Add(cur_histo_sig);

                //std::cout<<"sofarsogood in signal part"<<std::endl;
                delete cur_histo_sig;
	    }
	}
        delete file;
    }
//======================Add up background histograms==============================
    int count_bg = 0;
    while(fList_bg.getline(linefromlist, 250)){
        TFile* file = new TFile(linefromlist, "read");
        count_bg++;
        if(count_bg%100==0) std::cout<< count_bg<<" files processed"<<std::endl;
        //===========Invariant Mass Histograms===================================
	for(int i = 0; i < 4; i++){
	    for(int j = 0; j < 5; j++){
		char hist_name_bg[100];
		sprintf(hist_name_bg, "bg_omgpt%dcent%d", j+1, i+1);

		char cur_hist_name_bg[100];
		sprintf(cur_hist_name_bg, "sig_omgpt%dcent%d", j+1, i+1);

		TH1F* cur_histo_bg = (TH1F*)file -> Get(cur_hist_name_bg);

                TH1F* hist_bg = (TH1F*)HList.FindObject(hist_name_bg);
                
                if(cur_histo_bg) hist_bg -> Add(cur_histo_bg);
                if(cur_histo_bg) hist_total_bg -> Add(cur_histo_bg);
                //std::cout<<"sofarsogood in background part"<<std::endl;
                delete cur_histo_bg;
	    }
	}
        delete file;
    }


//=============Define scale factor=====================
    int low_bin = 21;
    int high_bin = 42;
//=============Scale Background and Signal==============
//==================Plot the Sig with Bg and Bg============================
    for(int i = 0; i < 4; i++){
	for(int j = 0; j < 5; j++){
	    char hist_name_sig[100];
	    char hist_name_bg[100];
	    sprintf(hist_name_sig, "sig_omgpt%dcent%d", j+1, i+1);
	    sprintf(hist_name_bg, "bg_omgpt%dcent%d", j+1, i+1);

	    TH1F* hist_sig = (TH1F*)HList.FindObject(hist_name_sig);

	    TH1F* hist_bg = (TH1F*)HList.FindObject(hist_name_bg);
     
            //hist_sig -> Sumw2();
            //hist_bg -> Sumw2();

            TCanvas* new_canvas = new TCanvas("c1", "c1", 200, 100, 600, 400);
            new_canvas -> cd();
             
            double int_sig = hist_sig -> Integral(low_bin, high_bin); 
            double int_bg = hist_bg -> Integral(low_bin, high_bin); 

            double scale_factor = int_sig / int_bg;
            hist_bg -> Scale(scale_factor);
            hist_bg -> SetLineColor(2);
            
            hist_sig -> Draw("HIST"); 
            hist_bg -> Draw("sames HIST");

            char figure_name_1[100];
            char figure_name_2[100];
            sprintf(figure_name_1, "../../plots/omg/omgpt%dcent%d.png", j+1, i+1);
            sprintf(figure_name_2, "../../plots/omg/omgpt%dcent%d.eps", j+1, i+1);
            new_canvas -> SaveAs(figure_name_1);
            new_canvas -> SaveAs(figure_name_2);
            
            delete new_canvas;
	}
    }
//=============Plot total result =======================================================
    double int_total_sig = hist_total_sig -> Integral(low_bin, high_bin); 
    double int_total_bg = hist_total_bg -> Integral(low_bin, high_bin); 
    double scale_total_factor = int_total_sig / int_total_bg;

    TCanvas* canvas_total = new TCanvas("c1", "c1", 200, 100, 600, 400);
    hist_total_bg -> Scale(scale_total_factor);
    hist_total_bg -> SetLineColor(2);
    hist_total_sig -> Draw("HIST");
    hist_total_bg -> Draw("sames HIST");
    canvas_total -> SaveAs("../../plots/omg/total_mix.png");
    
    hist_total_sig -> Add(hist_total_bg, -1);
    hist_total_sig -> SetMarkerStyle(4);
    TCanvas* canvas_total_pure = new TCanvas("c1", "c1", 200, 100, 600, 400);
    hist_total_sig -> Draw("HIST");
    canvas_total_pure -> SaveAs("../../plots/omg/total_pure.png");


//===============Get the pure signal by subtracting the background===============
//===============Plot the Invariant Mass Figures(After Subtraction)=================================
    float sig_count[4][5] = {};
    float sig_count_cal[4][5] = {};
    for(int i = 0; i < 4; i++){
	for(int j = 0; j < 5; j++){
	    char hist_name_sig[100];
	    char hist_name_bg[100];
	    sprintf(hist_name_sig, "sig_omgpt%dcent%d", j+1, i+1);
	    sprintf(hist_name_bg, "bg_omgpt%dcent%d", j+1, i+1);

	    TH1F* hist_sig = (TH1F*)HList.FindObject(hist_name_sig);

	    TH1F* hist_bg = (TH1F*)HList.FindObject(hist_name_bg);

            //hist_sig -> Sumw2(); 
            //hist_bg -> Sumw2();

	    TCanvas* new_canvas = new TCanvas("c1", "c1", 200, 100, 600, 400);
	    new_canvas -> cd();

	    double int_sig = hist_sig -> Integral(low_bin, high_bin); 
	    double int_bg = hist_bg -> Integral(low_bin, high_bin); 

	    double scale_factor = int_sig / int_bg;
	    hist_bg -> Scale(scale_factor);
	    hist_bg -> SetLineColor(2);

	    hist_sig -> Add(hist_bg, -1);
	    hist_sig -> SetMarkerStyle(4); 
            hist_sig -> Draw("HIST");
               
            //=================Fitting Ends=====================================
            
	    char figure_name_1[100];
	    char figure_name_2[100];
	    sprintf(figure_name_1, "../../plots/omg/pure_omgpt%dcent%d.png", j+1, i+1);
	    sprintf(figure_name_2, "../../plots/omg/pure_omgpt%dcent%d.eps", j+1, i+1);
	    new_canvas -> SaveAs(figure_name_1);
	    new_canvas -> SaveAs(figure_name_2);

            //=================Some Calculation==================================
            /*
            sig_count[i][j] = sig_func -> Integral(1.660, 1.684);  
            sig_count_cal[i][j] = par[2] / 0.0014;  
            */
	    delete new_canvas;
	}
    }

    //plot_pt_spectra(hSelectedMult, sig_count);
    output_histo -> cd();
    hNRefMult -> Write();
    hSelectedMult -> Write();
    hVz_before -> Write();
    hVz_after -> Write();
    hVr_before -> Write();
    hVr_after -> Write();
    hDeDx -> Write();
    hNRefMult -> Write();
    hSelectedMult -> Write();
    HList.Write();
    output_histo -> Close();
    output_file.close();
}


void plot_pt_spectra(TH1F* hSelectedMult, float** sig_count){
    float PI = 3.1415926;
    float pt_bin[12] = {0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.3, 1.5, 2.0, 2.5, 3.5, 5.0}; 
    float x_pt_spectra[11] = {0.45, 0.55, 0.65, 0.75, 0.9, 1.15, 1.5, 1.65, 2.25, 3.0, 4.25};
    float dpt_spectra[11] = {0.1, 0.1, 0.1, 0.1, 0.2, 0.3, 0.4, 0.3, 0.5, 1.0, 1.5};
    float y_pt_spectra[9][11] = {};
//calculate the y coordinates of points
    TCanvas* cpt_phi = new TCanvas("cpt_phi", "cpt_phi", 200, 10, 600, 400);
    cpt_phi -> SetLogy();

    std::cout<<"calculation started"<<std::endl;
    for(int i = 0; i < 9; i++){
        for(int j = 0; j < 11; j++){
            int nEvents_total = get_nevents_cent(i+1, hSelectedMult);
            double scale_factor = TMath::Power(10, -i);
            y_pt_spectra[i][j] = scale_factor*sig_count[i][j]/2/PI/x_pt_spectra[j]/dpt_spectra[j]/nEvents_total/2; 
        }
    }

    std::cout<<"calculation done"<<std::endl;
    for(int i = 0; i < 9; i++){
        std::cout<<"graph no."<<i<<std::endl;
        TGraph* cur_g = new TGraph(11, x_pt_spectra, y_pt_spectra[i]);
	cur_g -> SetMarkerSize(1.5);
	cur_g -> SetMarkerStyle(20+i);
	cur_g -> SetMarkerColor(2+i);
	cur_g -> SetMaximum(10.0);
        double scale_factor = Power(10, -i);
	cur_g -> SetMinimum(1E-14);
        cpt_phi -> cd();
	if(i == 0){ 
	    cur_g -> Draw("AP");
	}
        else{
            cur_g -> Draw("P");
        }
    }
    cpt -> SaveAs("phi_pt_spectra.eps");
    cpt -> SaveAs("phi_pt_spectra.png");
}

int get_nevents_cent(int nCent, TH1F* hSelectedMult){
    
    int low_b = 0;
    int high_b = 0;
    switch(nCent){
        case 1:
	    low_b = 198; 
	    high_b = 1000;
	    break;
	case 2:
	    low_b = 97;
	    high_b = 66;
	    break;
	case 3:
	    low_b = 66;
	    high_b = 0;
	    break;
	default:
	    break;
    } 

    int low_bin = hSelectedMult -> GetBin(low_b); 
    int high_bin = hSelectedMult -> GetBin(high_b);

    int nEvents_cent = hSelectedMult -> Integral(low_bin, high_bin);
    return nEvents_cent;
}
