int plot_phi(){
    gStyle -> SetOptStat();
    gStyle -> SetOptFit();
    TFile* output_histo = new TFile("phi_output.histo.root", "recreate");
    float mV0 = 1.01945;
    TFile* file = NULL;
    //ifstream fList("./total_histo_list_phi");
    ifstream fList("./total_list_1001");
    ofstream output_file("./output_1001");

    const float pt_bin[] = {0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.3, 1.7, 2.0, 2.5};//, 3.5, 5.0}; 
    const int cent_no_bin[] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 100};
    char linefromlist[250];

    const int centbin_no = 7;
    const int ptbin_no = 9;
//===========Initialize TObjArray to Store Histograms========
    TObjArray HList(0);
    
    for(int i = 0; i < centbin_no; i++){
	for(int j = 0; j < ptbin_no; j++){
	    char hist_name_sig[100];
	    char hist_des_sig[100];
	    char hist_name_bg[100];
	    char hist_des_bg[100];

	    sprintf(hist_name_sig, "sig_phipt%dcent%d", j+1, i+1);
	    sprintf(hist_des_sig, "sig_phi(%.1f< pt<%.1f)cent%d", pt_bin[j], pt_bin[j+1], i+1);
	    sprintf(hist_name_bg, "phipt%dcent%d", j+1, i+1);
	    sprintf(hist_des_bg, "phi(%.1f<pt<%.1f)cent%d", pt_bin[j], pt_bin[j+1], i+1);

	    TH1F* new_histo_sig = new TH1F(hist_name_sig, hist_des_sig, 100, mV0 - 0.07, mV0+ 0.07);
            new_histo_sig -> Sumw2();
	    TH1F* new_histo_bg = new TH1F(hist_name_bg, hist_des_bg, 100, mV0 - 0.07, mV0+ 0.07);
            new_histo_bg -> Sumw2();

	    HList.Add(new_histo_sig);
	    HList.Add(new_histo_bg);
	}
    }

//===========Input files listed in filelist===================
//===========QA Plots========================================
    TH1F* hNRefMult = new TH1F("RefMult", "Reference Multiplicity", 1000, 0.0, 1000.0 ) ;
    TH1F* hSelectedMult = new TH1F("SelectedRefMult", "Selected Reference Multiplicity", 1000, 0.0, 1000.0 ) ;
    TH1F* hVz_before = new TH1F("vz_before", "vz_before_cut", 400, -100, 100);
    TH1F* hVz_after = new TH1F("vz_after", "vz_after_cut", 400, -100, 100);
    TH2F* hVr_before = new TH2F("vr_before", "vr_before_cut", 200, -5.0, 5.0, 200, -5.0, 5.0);
    TH2F* hVr_after = new TH2F("vr_after", "vr_after_cut", 200, -5.0, 5.0, 200, -5.0, 5.0);
    TH2F* hDeDx = new TH2F("dedxP", "DeDx_P", 300, -3, 3, 300, 0, 3e-5);

    //while(fList.getline(linefromlist, 250)){
    file = new TFile("phi_BBCMONTOF.histo.root", "read"); // Summarized Histogram ROOT File
    TH1F* hNRefMult  =  (TH1F*)file -> Get("RefMult");
    TH1F* hSelectedMult = (TH1F*)file -> Get("SelectRefMult");
    TH1F* hVz_before = (TH1F*)file -> Get("VertexZ");
    TH1F* hVz_after = (TH1F*)file -> Get("VertexZ_after");
    TH2F* hVr_before = (TH2F*)file -> Get("VertexR");
    TH2F* hVr_after = (TH2F*)file -> Get("VertexR_after");
    TH2F* hDeDx = (TH2F*)file -> Get("DedxP");

    //===========Invariant Mass Histograms===================================
    for(int i = 0; i < centbin_no; i++){
	for(int j = 0; j < ptbin_no; j++){
	    char hist_name_sig[100];
	    char hist_name_bg[100];
	    sprintf(hist_name_sig, "sig_phipt%dcent%d", j+1, i+1);
	    sprintf(hist_name_bg, "phipt%dcent%d", j+1, i+1);

	    TH1F* cur_histo_sig = (TH1F*)file -> Get(hist_name_sig);
	    TH1F* cur_histo_bg = (TH1F*)file -> Get(hist_name_bg);

	    TH1F* hist_sig = (TH1F*)HList.FindObject(hist_name_sig);
	    hist_sig -> Add(cur_histo_sig);

	    TH1F* hist_bg = (TH1F*)HList.FindObject(hist_name_bg);
	    hist_bg -> Add(cur_histo_bg);
	    delete cur_histo_sig;
	    delete cur_histo_bg;
	}
    }
        //delete file;
    //}

//=============Define scale factor=====================
    int low_b = 1.04; //0.99;
    int high_b = 1.06; //1.01;
//=============Scale Background and Signal==============
//==================Plot the Sig with Bg and Bg============================
    for(int i = 0; i < centbin_no; i++){ // Centrality Bin
	for(int j = 0; j < ptbin_no; j++){ // Pt Bin
	    char hist_name_sig[100];
	    char hist_name_bg[100];
	    sprintf(hist_name_sig, "sig_phipt%dcent%d", j+1, i+1);
	    sprintf(hist_name_bg, "phipt%dcent%d", j+1, i+1);

	    TH1F* hist_sig = (TH1F*)HList.FindObject(hist_name_sig);
	    TH1F* hist_bg = (TH1F*)HList.FindObject(hist_name_bg);

            TCanvas* new_canvas = new TCanvas("c1", "c1", 200, 100, 600, 400);
            new_canvas -> cd();
             
            int low_bin = 70; //hist_sig -> FindBin(low_b);//TODO:
            int high_bin = 90; //hist_sig -> FindBin(high_b);
       
            double int_sig = hist_sig -> Integral(low_bin, high_bin); 
            double int_bg = hist_bg -> Integral(low_bin, high_bin); 

            double scale_factor = int_sig / int_bg;
            hist_bg -> Scale(scale_factor);
            hist_bg -> SetLineColor(2);
            
            hist_sig -> Draw("E"); 
            hist_bg -> Draw("sames E");

            char figure_name_png[100];
            char figure_name_eps[100];
            sprintf(figure_name_png, "../phi_plots/phipt%dcent%d.png", j+1, i+1);
            sprintf(figure_name_eps, "../phi_plots/phipt%dcent%d.eps", j+1, i+1);
            new_canvas -> SaveAs(figure_name_png);
            new_canvas -> SaveAs(figure_name_eps);
            
            delete new_canvas;
	}
    }

//===============Get the pure signal by subtracting the background===============
//===============Plot the Invariant Mass Figures(After Subtraction)=================================
    float sig_count[centbin_no][ptbin_no] = {};
    float sig_count_byfit[centbin_no][ptbin_no] = {};
    float bg_count[centbin_no][ptbin_no] = {};
    float sig_count_cal[centbin_no][ptbin_no] = {};
    for(int i = 0; i < centbin_no; i++){
	for(int j = 0; j < ptbin_no; j++){
	    char hist_name_sig[100];
	    char hist_name_bg[100];
	    sprintf(hist_name_sig, "sig_phipt%dcent%d", j+1, i+1);
	    sprintf(hist_name_bg, "phipt%dcent%d", j+1, i+1);

	    TH1F* hist_sig = (TH1F*)HList.FindObject(hist_name_sig);
	    TH1F* hist_bg = (TH1F*)HList.FindObject(hist_name_bg);

	    TCanvas* new_canvas = new TCanvas("c1", "c1", 200, 100, 600, 400);
	    new_canvas -> cd();

	    int low_bin = hist_sig -> FindBin(low_b);//TODO:
            int high_bin = hist_sig -> FindBin(high_b);

	    double int_sig = hist_sig -> Integral(low_bin, high_bin); 
	    double int_bg = hist_bg -> Integral(low_bin, high_bin); 

	    double scale_factor = int_sig / int_bg;
	    hist_bg -> Scale(scale_factor);
	    hist_bg -> SetLineColor(2);

	    hist_sig -> Add(hist_bg, -1);

            //=================Fitting Starts===============================================
	    TF1* total_func = new TF1("total_func", "[0] + [1] * x + 1/(2 * 3.1415926) * [2] * [3] / ((x - [4]) * (x - [4]) + [3] * [3] / 4)", 0.99, 1.05);
	    TF1* sig_func = new TF1("BW_func", "1/(2 * 3.1415926) * [0] * [1] / ((x - [2]) * (x - [2]) + [1] * [1] / 4)",  0.99, 1.05);
	    TF1* bg_func = new TF1("bg_func", "[0]+[1]*x", 0.99, 1.05);

	    total_func -> SetParName(0, "p0");
	    total_func -> SetParName(1, "p1");
	    total_func -> SetParName(2, "BW Area");
	    total_func -> SetParName(3, "#Gamma");
	    total_func -> SetParName(4, "M_{0}");
	    total_func -> SetParameter(3, 0.006);
	    total_func -> SetParameter(4, 1.0195);
            total_func -> SetParLimits(3, 0.004, 0.007);
            total_func -> SetParLimits(4, 1.015, 1.023);
	    total_func -> SetLineColor(2);
            hist_sig -> Fit("total_func", "REM"); 
	    hist_sig -> SetMarkerStyle(4); 
	    hist_sig -> Draw("E");

            double par[5];//to accept fitting parameters
  
            total_func -> GetParameters(par);

            sig_func -> SetParameter(0, par[2]);
            sig_func -> SetParameter(1, par[3]);
            sig_func -> SetParameter(2, par[4]);
            sig_func -> SetLineColor(4);
            sig_func -> Draw("same");

            bg_func -> SetParameter(0, par[0]);
            bg_func -> SetParameter(1, par[1]);
            bg_func -> SetLineColor(3);
            bg_func -> Draw("same");
            //=================Fitting Ends=====================================
            
	    char figure_name_png[100];
	    char figure_name_eps[100];
	    sprintf(figure_name_png, "../phi_plots/pure_phipt%dcent%d.png", j+1, i+1);
	    sprintf(figure_name_eps, "../phi_plots/pure_phipt%dcent%d.eps", j+1, i+1);
	    new_canvas -> SaveAs(figure_name_png);
	    new_canvas -> SaveAs(figure_name_eps);

            //=================Draw Integration Region==========================//TODO:
            float x_low_line = mV0 + 0.007; 
            float x_up_line = mV0 - 0.007; 
            float peak_value = hist_sig -> GetBinContent(hist_sig -> GetMaximumBin()); 
            float valley_value = hist_sig -> GetBinContent(hist_sig -> GetMinimumBin()); 
           
            TLine* l1 = new TLine(x_low_line, valley_value, x_low_line, 0.5*peak_value); 
            TLine* l2 = new TLine(x_up_line, valley_value, x_up_line, 0.5*peak_value); 
            l1 -> SetLineColor(4);
            l1 -> SetLineWidth(3);
            l1 -> SetLineStyle(4);
	    l2 -> SetLineColor(4);
            l2 -> SetLineWidth(3);
            l2 -> SetLineStyle(4);

            l1 -> Draw("sames");
            l1 -> Draw("sames");

            low_bin = hist_sig -> FindBin(mV0 - 0.007);
            high_bin = hist_sig -> FindBin(mV0 + 0.007);
            sig_count[i][j] = hist_sig -> Integral(low_bin, high_bin);  //sig_func -> Integral(0.99, 1.05);
            sig_count_cal[i][j] = hist_sig -> Integral(low_bin, high_bin) - bg_func -> Integral(mV0 - 0.007, mV0 + 0.007)/0.0014;
            bg_count[i][j] = bg_func -> Integral(mV0 - 0.007, mV0 + 0.007)/0.0014;
            sig_count_byfit[i][j] = sig_func -> Integral(mV0 - 0.007, mV0 + 0.007)/0.0014;
	    delete new_canvas;
	}
    }
//==============Plot Spectra========================================
    float PI = 3.1415926;
    //const float pt_bin[10] = {0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.3, 1.7, 2.0, 2.5};//, 3.5, 5.0}; 
    float x_pt_spectra[] = {0.45, 0.55, 0.65, 0.75, 0.9, 1.15, 1.5, 1.85, 2.25}; //, 3.0, 4.25};
    float dpt_spectra[] = {0.1, 0.1, 0.1, 0.1, 0.2, 0.3, 0.4, 0.3, 0.5}; //, 1.0, 1.5};
    float y_pt_spectra[centbin_no][ptbin_no] = {};
    float y_pt_spectra_cal[centbin_no][ptbin_no] = {};

//calculate the y coordinates of points
    std::cout<<"====Calculation Starts===="<<std::endl;
    for(int i = 0; i < centbin_no; i++){
            int nEvents_total = get_nevents_cent(i+1, hSelectedMult);
            output_file << nEvents_total << " found in cent " << i+1 << "\n";
        for(int j = 0; j < ptbin_no; j++){
            double scale_factor = TMath::Power(10, -i);
            y_pt_spectra[i][j] = scale_factor * sig_count_byfit[i][j] / 2 / PI / x_pt_spectra[j] / dpt_spectra[j] / nEvents_total;
            y_pt_spectra_cal[i][j] = scale_factor * sig_count_cal[i][j] / 2 / PI / x_pt_spectra[j] / dpt_spectra[j] / nEvents_total;
            output_file << i+1 << " " << j+1 << " " << y_pt_spectra_cal[i][j] << " " << y_pt_spectra[i][j] << " " << sig_count_cal[i][j] << " " << sig_count[i][j] << " " << bg_count[i][j] << " " << sig_count_byfit[i][j] << " " << nEvents_total << std::endl;
        }
    }

    std::cout<<"====Calculation Done===="<<std::endl;

    TCanvas* cpt_phi = new TCanvas("cpt_phi", "cpt_phi", 200, 10, 600, 400);
    cpt_phi -> SetLogy();
    TLegend* leg1 = new TLegend(0.65, 0.65, 0.89, 0.89);
    for(int i = 0; i < centbin_no; i++){
        std::cout<<"graph no."<<i<<std::endl;
        TGraph* cur_g = new TGraph(ptbin_no, x_pt_spectra, y_pt_spectra[i]);
	cur_g -> SetMarkerSize(1.5);
	cur_g -> SetMarkerStyle(20);
	cur_g -> SetMarkerColor(2+i);
	cur_g -> SetMaximum(10.0);
	cur_g -> SetMinimum(1E-20);
        cpt_phi -> cd();
        cur_g -> SetTitle("#phi Meson Raw Spectra@14.5GeV");
        cur_g -> GetXaxis() -> SetTitle("P_{t}(GeV/c)");
        cur_g -> GetYaxis() -> SetTitle("#frac{d^{2}N}{2#piPtNdPtdy}(GeV/c)^{2}");
        cur_g -> GetYaxis() -> SetTitleOffset(1.4);
        cur_g -> GetYaxis() -> CenterTitle();
	if(i == 0){ 
	    cur_g -> Draw("AP");
	}
        else{
            cur_g -> Draw("P same");
        }
	char leg_text[100];
        sprintf(leg_text, "%d%%~%d%%(10E-%d)", cent_no_bin[i], cent_no_bin[i+1], i);
        leg1 -> AddEntry(cur_g, leg_text, "p");
    }
    leg1 -> Draw("same");
    cpt_phi -> SaveAs("../phi_plots/phi_pt_spectra.eps");
    cpt_phi -> SaveAs("../phi_plots/phi_pt_spectra.png");

    TCanvas* cpt_phi_cal = new TCanvas("cpt_phi_cal", "cpt_phi_cal", 200, 10, 600, 400);
    cpt_phi_cal -> SetLogy();
    TLegend* leg = new TLegend(0.65, 0.65, 0.89, 0.89);
    for(int i = 0; i < centbin_no; i++){
        if(i == 8) continue;
        std::cout<<"graph no."<<i<<std::endl;
        TGraph* cur_g = new TGraph(ptbin_no, x_pt_spectra, y_pt_spectra_cal[i]);
	cur_g -> SetMarkerSize(1.5);
	cur_g -> SetMarkerStyle(20);
	cur_g -> SetMarkerColor(2+i);
	cur_g -> SetMaximum(10.0);
	cur_g -> SetMinimum(1E-20);
        cpt_phi_cal -> cd();
	cur_g -> SetTitle("#phi Meson Raw Spectra@14.5GeV");
        cur_g -> GetXaxis() -> SetTitle("P_{t}(GeV/c)");
        cur_g -> GetYaxis() -> SetTitle("#frac{d^{2}N}{2#piPtNdPtdy}(GeV/c)^{2}");
        cur_g -> GetYaxis() -> SetTitleOffset(1.4);
        cur_g -> GetYaxis() -> CenterTitle();
	if(i == 0){ 
	    cur_g -> Draw("AP");
	}
        else{
            cur_g -> Draw("P");
        }
        char leg_text[100];
        sprintf(leg_text, "%d%%~%d%%(10E-%d)", cent_no_bin[i], cent_no_bin[i+1], i);
        leg -> AddEntry(cur_g, leg_text, "p");
    }
    leg -> Draw("same");
    cpt_phi_cal -> SaveAs("../phi_plots/phi_pt_spectra_cal.eps");
    cpt_phi_cal -> SaveAs("../phi_plots/phi_pt_spectra_cal.png");

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

/*
void plot_pt_spectra(TH1F* hSelectedMult, float** sig_count){
    float PI = 3.1415926;
    float pt_bin[12] = {0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.3, 1.7, 2.0, 2.5, 3.5, 5.0}; 
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
*/

int get_nevents_cent(int nCent, TH1F* hSelectedMult){
    
    int low_b = 0;
    int high_b = 0;
    switch(nCent){ // Based on http://www.star.bnl.gov/protected/lfspectra/jdb/run14/AuAu15/RefMultCorr/Run14AuAu15_Centrality_Systematics_May_15_15.pdf
        case 1:
	    low_b = 239; 
	    high_b = 1000;
	    break;
	case 2:
	    low_b = 200;
	    high_b = 239;
	    break;
	case 3:
	    low_b = 139;
	    high_b = 200;
	    break;
	case 4:
	    low_b = 94;
	    high_b =139;
	    break;
	case 5:
	    low_b = 61;
	    high_b = 94;
	    break;
	case 6:
	    low_b = 37;
	    high_b = 61;
	    break;
	case 7:
	    low_b = 21;
	    high_b = 37;
	    break;
	case 8:
	    low_b = 11;
	    high_b = 21;
	    break;
	case 9:
	    low_b = 0;
	    high_b = 11;
	    break;
	default:
	    break;
    } 

    int low_bin = hSelectedMult -> GetBin(low_b); 
    int high_bin = hSelectedMult -> GetBin(high_b);

    int nEvents_cent = hSelectedMult -> Integral(low_bin, high_bin);
    return nEvents_cent;
}

