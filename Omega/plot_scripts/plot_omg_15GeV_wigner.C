int plot_omg_15GeV_wigner(std::string particle){
    gStyle -> SetOptFit(111);
    double pdgmass_omg = 1.67245;
    TFile* infile;
    if(particle == "omg"){
	infile = new TFile("omg_21M_BBC_BBCMONTOF.local.root", "read");
    }
    else if(particle == "antiomg"){
	infile = new TFile("antiomg_21M_BBC_BBCMONTOF.local.root", "read");
    }
    else{
        std::cout<<"Wrong Particle Type!"<<std::endl;
        return -1;
    }
    

    Float_t pt_bin[8] = {0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.6};
    const Int_t no_centbin = 9;
    const Int_t no_ptbin = 7;

    //Float_t lb = 1.65;
    //Float_t ub = 1.70;
    /*
    for(int i = 0; i < 9; i++){
        for(int j = 0; j < 6; j++){
	    char hist_name_sig[200];
            char can_name[200];
            sprintf(hist_name_sig, "sig_xipt%dcent%d", j+1, i+1);
            sprintf(can_name, "../plots/sig_xipt%dcent%d.eps", j+1, i+1);

	    TH1F* h = (TH1F*) infile -> Get(hist_name_sig);
	    TCanvas* c = new TCanvas("c");

	    h -> Draw();
	    h -> GetXaxis() -> SetTitle("inv_mass(GeV)");
	    h -> GetYaxis() -> SetTitle("counts");

	    f1 -> SetParameter(3, 0.007);
	    f1 -> SetParameter(4, 1.675);
	    f1 -> SetParameter(5, 0.01);

	    h -> Fit("f1", "REM");
            c -> SaveAs(can_name);
            delete c;

            char hist_name_sig[200];
            sprintf(hist_name_sig, "sig_xipt%dcent%d", j+1, i+1);

            TH1F* h = (TH1F*) infile -> Get(hist_name_sig);
            TCanvas* c = new TCanvas("c");
 //     
            h -> Draw();
            h -> GetXaxis() -> SetTitle("inv_mass(GeV)");
            h -> GetYaxis() -> SetTitle("counts");

            int lb_bin = h -> FindBin(lb); 
            int ub_bin = h -> FindBin(ub);
            double int_count = h -> Integral(lb_bin, ub_bin);
            f1 -> SetParameter(2, int_count);
            f1 -> SetParameter(3, pdgmass_omg);
            f1 -> SetParameter(4, 0.03);
    
            h -> Fit("gaus", "EM");
//
        }
    }
*/
    double fit_par_010[6][6];
    double fit_err_010[6][6];
    double fit_par_1060[6][6];
    double fit_err_1060[6][6];
    double fit_mass_010[6];
    double fit_masserr_010[6];
    double fit_mass_1060[6];
    double fit_masserr_1060[6];
    double fit_sigma_010[6];
    double fit_sigmaerr_010[6];
    double fit_sigma_1060[6];
    double fit_sigmaerr_1060[6];
    double sig_counts_010[6];
    double sig_counts_1060[6];
    double int_l = pdgmass_omg - 0.006; 
    double int_u = pdgmass_omg + 0.006;
    double lb = 1.64;
    double ub = 1.704;
    TF1* f1 = new TF1("f1", "1/(2 * 3.1415926) * [0] * [1] / ((x - [2]) * (x - [2]) + [1] * [1] / 4) + [3] + [4]*x + [5]*x*x", lb, ub);//lb, ub);
    TF1* f_sig = new TF1("f_sig", "1/(2 * 3.1415926) * [0] * [1] / ((x - [2]) * (x - [2]) + [1] * [1] / 4)", lb, ub);
    TF1* f_bg = new TF1("f_bg", "[0]+[1]*x+[2]*x*x", lb, ub);

    for(int i = 0; i < 6; i++){
        char hist_name_sig_010[200];
        char hist_name_sig_1060[200];
        sprintf(hist_name_sig_010, "sig_xipt%dcent_010", i+1);
        sprintf(hist_name_sig_1060, "sig_xipt%dcent_1060", i+1);
        char can_name_sig_010[200];
        char can_name_sig_1060[200];
        if(particle == "omg"){
	    sprintf(can_name_sig_010, "../omg_plots/sig_xipt%dcent_010.eps", i+1);
	    sprintf(can_name_sig_1060, "../omg_plots/sig_xipt%dcent_1060.eps", i+1);
	}
        else{
	    sprintf(can_name_sig_010, "../antiomg_plots/sig_xipt%dcent_010.eps", i+1);
	    sprintf(can_name_sig_1060, "../antiomg_plots/sig_xipt%dcent_1060.eps", i+1);
	}
   
        TH1F* h_010 = (TH1F*) infile -> Get(hist_name_sig_010);
        TH1F* h_1060 = (TH1F*) infile -> Get(hist_name_sig_1060);

        TCanvas* c_010 = new TCanvas("c_010"); 
        TCanvas* c_1060 = new TCanvas("c_1060"); 

        c_010 -> cd();
        //h_010 -> Rebin();
        h_010 -> SetMarkerStyle(8);
        h_010 -> Draw("PE");
        h_010 -> GetXaxis() -> SetTitle("inv_mass(GeV)");
	h_010 -> GetYaxis() -> SetTitle("counts");
 
	f1 -> SetParameter(0, 0.007);
	//f1 -> SetParameter(1, 0.001);
	f1 -> SetParameter(1, 0.01);
	f1 -> SetParameter(2, 1.671);
        //if(particle == "omg" && (i == 0 || i== 3)) f1 -> SetParameter(1, 0.01);
        if(particle == "omg" && (i == 0 || i== 3)){
	    f1 -> SetParameter(1, 0.01);
	    f1 -> SetParameter(2, 1.672);
        }
        else if(particle == "omg" && i == 5){
	    f1 -> SetParameter(1, 0.001);
	    f1 -> SetParameter(2, 1.672);
	}
        else if(particle == "antiomg" && i == 5){
            h_010 -> Rebin();
	    f1 -> SetParameter(1, 0.05);
	    f1 -> SetParameter(2, 1.671);
        }
        else if(particle == "antiomg" && i == 0){
            h_010 -> Rebin();
	    f1 -> SetParameter(2, 1.671);
	    f1 -> SetParameter(1, 0.05);
        }
        else if(particle == "antiomg" && i == 2){
	    f1 -> SetParameter(2, 1.671);
	    f1 -> SetParameter(1, 0.05);
        }
        else{
	    f1 -> SetParameter(0, 0.007);
	    //f1 -> SetParameter(1, 0.001);
	    f1 -> SetParameter(1, 0.01);
	    f1 -> SetParameter(2, 1.671);
	}
	h_010 -> Fit("f1", "REM");

    
        f1 -> GetParameters(fit_par_010[i]);
        //fit_err_010[i] = f1 -> GetParErrors();
        fit_mass_010[i] = fit_par_010[i][2];
        fit_masserr_010[i] = f1 -> GetParError(2);//fit_err_010[i][2];
        fit_sigma_010[i] = fit_par_010[i][1];
        fit_sigmaerr_010[i] = f1 -> GetParError(1);//fit_err_010[i][1];
        int_l = fit_par_010[i][2] - 0.006;//2*fit_par_010[i][1];
        int_u = fit_par_010[i][2] + 0.006;//2*fit_par_010[i][1];

	TLine* line_l = new TLine(int_l, 0, int_l, h_010 -> GetMaximum());
        TLine* line_u = new TLine(int_u, 0, int_u, h_010 -> GetMaximum());
        line_l -> SetLineColor(3);
        line_l -> SetLineWidth(2);
        line_l -> SetLineStyle(10);
        line_u -> SetLineColor(3);
        line_u -> SetLineWidth(2);
        line_u -> SetLineStyle(10);
        line_l -> Draw("sames");
        line_u -> Draw("sames");

        f_sig -> SetParameter(0, fit_par_010[i][0]);
        f_sig -> SetParameter(1, fit_par_010[i][1]);
        f_sig -> SetParameter(2, fit_par_010[i][2]);
        f_bg -> SetParameter(0, fit_par_010[i][3]);
        f_bg -> SetParameter(1, fit_par_010[i][4]);
        f_bg -> SetParameter(2, fit_par_010[i][5]);
        float bg_counts_010 = f_bg -> Integral(int_l, int_u)/0.0014;
        if(particle == "antiomg" && (i == 5 || i == 0)){
            bg_counts_010 = f_bg -> Integral(int_l, int_u)/0.0028;
        }
        int lb_bin = h_010 -> FindBin(int_l);
        int ub_bin = h_010 -> FindBin(int_u);
        sig_counts_010[i] = h_010 -> Integral(lb_bin, ub_bin) - bg_counts_010;
        cout<<sig_counts_010[i] << "======sig_counts" <<endl;

        c_010 -> SaveAs(can_name_sig_010);

	//================================================================
        c_1060 -> cd();
        if(particle == "omg" && i == 5){
            h_1060->Rebin();
	    f1 -> SetParameter(2, 1.671);
	    f1 -> SetParameter(1, 0.005);
	}
        else if(particle == "antiomg" && i == 5){
            f1 -> SetParameter(2, 1.671);
            f1 -> SetParameter(1, 0.005 );
        }
        else{
	    f1 -> SetParameter(0, 0.007);
	    //f1 -> SetParameter(1, 0.001);
	    f1 -> SetParameter(1, 0.001);
	    f1 -> SetParameter(2, 1.671);
	}
	h_1060 -> SetMarkerStyle(8);
        h_1060 -> Draw("PE");
	h_1060 -> GetXaxis() -> SetTitle("inv_mass(GeV)");
	h_1060 -> GetYaxis() -> SetTitle("counts");

        h_1060 -> Fit("f1", "REM");
        f1 -> GetParameters(fit_par_1060[i]);
	//fit_err_1060[i] = f1 -> GetParErrors();
        fit_mass_1060[i] = fit_par_1060[i][2];
        fit_masserr_1060[i] = f1 -> GetParError(2);//fit_err_1060[i][2];
        fit_sigma_1060[i] = fit_par_1060[i][1];
        fit_sigmaerr_1060[i] = f1 -> GetParError(1);//fit_err_1060[i][1];

	int_l = fit_par_1060[i][2] - 0.006;//2*fit_par_1060[i][1];
        int_u = fit_par_1060[i][2] + 0.006;//2*fit_par_1060[i][1];
	f_sig -> SetParameter(0, fit_par_1060[i][0]);
        f_sig -> SetParameter(1, fit_par_1060[i][1]);
        f_sig -> SetParameter(2, fit_par_1060[i][2]);
        f_bg -> SetParameter(0, fit_par_1060[i][3]);
        f_bg -> SetParameter(1, fit_par_1060[i][4]);
        f_bg -> SetParameter(2, fit_par_1060[i][5]);

 	TLine* line_l_1060 = new TLine(int_l, 0, int_l, h_1060 -> GetMaximum());
        TLine* line_u_1060 = new TLine(int_u, 0, int_u, h_1060 -> GetMaximum());
        line_l_1060 -> SetLineColor(3);
        line_l_1060 -> SetLineWidth(2);
        line_l_1060 -> SetLineStyle(10);
        line_u_1060 -> SetLineColor(3);
        line_u_1060 -> SetLineWidth(2);
        line_u_1060 -> SetLineStyle(10);
        line_l_1060 -> Draw("sames");
        line_u_1060 -> Draw("sames");

        bg_counts = f_bg -> Integral(int_l, int_u)/0.0014; 
        if(particle=="omg" && i==5){
	    bg_counts = f_bg -> Integral(int_l, int_u)/0.0028; 
        }
	lb_bin = h_1060 -> FindBin(int_l);
        ub_bin = h_1060 -> FindBin(int_u);
        sig_counts_1060[i] = h_1060-> Integral(lb_bin, ub_bin) - bg_counts;
        c_1060 -> SaveAs(can_name_sig_1060);

        delete c_010;
        delete c_1060;
        delete line_l;
        delete line_u;
        delete line_l_1060;
        delete line_u_1060;
    }

//====Plot Spectrum====
    //Calculation

    double PI = 3.1415926;
    int nevents[9] = {1.538113e6, 2.45574e6, 2.623553e6, 2.703106e6, 2.69095e6, 2.725661e6, 2.68739e6, 1.293578e6, 1.333159e6};
    double x_pt_spectra[6] = {0.95, 1.4, 1.8, 2.2, 2.6, 3.2 };
    double x_pterr_spectra[6] = {0, 0, 0, 0, 0, 0};
    double dpt_spectra[] = {0.5, 0.4, 0.4, 0.4, 0.4, 0.8 };
    double y_pt_spectra_010[6] = {};  
    double y_pt_spectra_1060[6] = {};  
    for(int j = 0; j < 6; j++){
	y_pt_spectra_1060[j] = 1/(2*PI) * sig_counts_1060[j] / x_pt_spectra[j] / dpt_spectra[j] / (nevents[6]+nevents[5] + nevents[4] + nevents[3] + nevents[2]); 
	y_pt_spectra_010[j] = 1/(2*PI) * sig_counts_010[j] / x_pt_spectra[j] / dpt_spectra[j]/(nevents[7] + nevents[8]); 
	//cout<<sig_counts_010[j] << "======y_pt_spectra_010 "<<y_pt_spectra_010[j]<< " nevents=6 -> "<<nevents[6]<<endl;
	printf("ypt = %.10f\n", y_pt_spectra_010[j]);
	printf("ypt = %.10f\n", y_pt_spectra_1060[j]);
	printf("sig_count = %.10f<->%.10f\n", sig_counts_010[j], sig_counts_1060[j]);
    }
    //Plotting
    TCanvas* cpt_omg_010 = new TCanvas("cpt_omg_010", "cpt_omg_010", 200, 10, 600, 400);
    cpt_omg_010 -> SetLogy();
    TGraph* cur_g = new TGraph(11, x_pt_spectra, y_pt_spectra_010);
    cur_g -> SetMarkerSize(1.5);
    cur_g -> SetMarkerStyle(20);
    cur_g -> SetMarkerColor(2);
    cur_g -> SetMaximum(10E-3);
    cur_g -> SetMinimum(10E-14);
    cur_g -> GetXaxis() -> SetLimits(0.5, 3.60);
    cur_g -> SetTitle("#Omega^{-} 0-10%@AuAu14.5GeV");
    if(particle == "antiomg")
	cur_g -> SetTitle("#Omega^{+} 0-10%@AuAu14.5GeV");
    cur_g -> GetYaxis() -> SetTitle("dN^{2}/(2#piN_{evnt}P_{t}dP_{t})(GeV/c)^{2}");
    cur_g -> GetXaxis() -> SetTitle("Pt(GeV/c)");
    cur_g -> Draw("AP");
    if(particle == "omg"){
	cpt_omg_010 -> SaveAs("../omg_plots/omg_pt_spectra_010.eps");
	cpt_omg_010 -> SaveAs("../omg_plots/omg_pt_spectra_010.png");
    }
    else{
	cpt_omg_010 -> SaveAs("../antiomg_plots/omg_pt_spectra_010.eps");
	cpt_omg_010 -> SaveAs("../antiomg_plots/omg_pt_spectra_010.png");
    }

    TCanvas* cpt_omg_1060 = new TCanvas("cpt_omg_1060", "cpt_omg_1060", 200, 10, 600, 400);
    cpt_omg_1060 -> SetLogy();
    TGraph* cur_g_1060 = new TGraph(11, x_pt_spectra, y_pt_spectra_1060);
    cur_g_1060 -> SetMarkerSize(1.5);
    cur_g_1060 -> SetMarkerStyle(20);
    cur_g_1060 -> SetMarkerColor(2);
    cur_g_1060 -> SetMaximum(10E-3);
    cur_g_1060 -> SetMinimum(10E-14);
    cur_g_1060 -> GetXaxis() -> SetLimits(0.5, 3.60);
    cur_g_1060 -> SetTitle("#Omega^{-} 10-60%@AuAu14.5GeV");
    if(particle == "antiomg")
	cur_g_1060 -> SetTitle("#Omega^{+} 10-60%@AuAu14.5GeV");
    cur_g_1060 -> GetYaxis() -> SetTitle("dN^{2}/(2#piN_{evnt}P_{t}dP_{t})(GeV/c)^{2}");
    cur_g_1060 -> GetXaxis() -> SetTitle("Pt(GeV/c)");
    cur_g_1060 -> Draw("AP");
    if(particle == "omg"){
	cpt_omg_1060 -> SaveAs("../omg_plots/omg_pt_spectra_1060.eps");
	cpt_omg_1060 -> SaveAs("../omg_plots/omg_pt_spectra_1060.png");
    }
    else{
	cpt_omg_1060 -> SaveAs("../antiomg_plots/omg_pt_spectra_1060.eps");
	cpt_omg_1060 -> SaveAs("../antiomg_plots/omg_pt_spectra_1060.png");
    }


/*
    for(int i = 0; i < 9; i++){
	std::cout<<"graph no."<<i<<std::endl;
	TGraph* cur_g = new TGraph(11, x_pt_spectra, y_pt_spectra[i]);
	cur_g -> SetMarkerSize(1.5);
	cur_g -> SetMarkerStyle(20);
	cur_g -> SetMarkerColor(2+i);
	cur_g -> SetMaximum(10.0);
	cur_g -> SetMinimum(1E-20);
	cpt_phi -> cd();
	if(i == 0){
	    cur_g -> Draw("AP");
	}
	else{
	    cur_g -> Draw("P");
	}
    }
*/

//====Plotting the overview figures====
    TFile* file_overview = new TFile("overview.histo.root", "read");
    TCanvas* c1 = new TCanvas("c1");
    TH1F* h_centbin9 = (TH1F*) file_overview -> Get("h_centbin9_after");
    h_centbin9 -> GetXaxis() -> SetTitle("centrality bin");
    h_centbin9 -> GetYaxis() -> SetTitle("counts");
    h_centbin9 -> Draw();
    c1 -> SaveAs("../omg_plots/centbin9.eps");

    TCanvas* c2 = new TCanvas("c2");
    c2 -> SetLogz();
    TH2F* h_vpr = (TH2F*) file_overview -> Get("h_vpr_after");
    h_vpr -> GetXaxis() -> SetTitle("x(cm)");
    h_vpr -> GetYaxis() -> SetTitle("y(cm)");
    h_vpr -> Draw("colorz"); 
    c2 -> SaveAs("../omg_plots/vpr.eps");

    TCanvas* c3 = new TCanvas("c3");
    TH1F* h_vpz = (TH1F*) file_overview -> Get("h_vpz_after");
    h_vpz -> GetXaxis() -> SetTitle("z(cm)");
    h_vpz -> GetYaxis() -> SetTitle("counts");
    h_vpz -> Draw();
    c3 -> SaveAs("../omg_plots/vpz.eps");
   
    TCanvas* c4 = new TCanvas("c4");
    TH1F* h_ximass = (TH1F*) infile -> Get("ximass");
    h_ximass -> SetTitle("#Omega^{-} Invariant Mass");
    if(particle == "antiomg")
	h_ximass -> SetTitle("#Omega^{+} Invariant Mass");
    h_ximass -> GetXaxis() -> SetTitle("Invariant Mass(GeV)");
    h_ximass -> GetYaxis() -> SetTitle("counts");
    h_ximass -> Draw();
    if(particle == "omg"){
	c4 -> SaveAs("../omg_plots/ximass.eps");
    }
    else{
	c4 -> SaveAs("../antiomg_plots/ximass.eps");
    }

    TCanvas* c5 = new TCanvas("mass_010");
    TGraphErrors* gr_mass_010 = new TGraphErrors(6, x_pt_spectra, fit_mass_010, x_pterr_spectra, fit_masserr_010);
    gr_mass_010 -> SetMarkerColor(2);
    gr_mass_010 -> SetMarkerStyle(21);
    gr_mass_010 -> Draw("AP"); 
    if(particle == "omg"){
	c5 -> SaveAs("../omg_plots/omg_mass_010.eps");
	c5 -> SaveAs("../omg_plots/omg_mass_010.jpg");
    }
    else if(particle == "antiomg"){
	c5 -> SaveAs("../antiomg_plots/antiomg_mass_010.eps");
	c5 -> SaveAs("../antiomg_plots/antiomg_mass_010.jpg");
    }

    TCanvas* c6 = new TCanvas("sigma_010");
    TGraphErrors* gr_sigma_010 = new TGraphErrors(6, x_pt_spectra, fit_sigma_010, x_pterr_spectra, fit_sigmaerr_010);
    gr_sigma_010 -> SetMarkerColor(2);
    gr_sigma_010 -> SetMarkerStyle(21);
    gr_sigma_010 -> Draw("AP"); 
    if(particle == "omg"){
	c6 -> SaveAs("../omg_plots/omg_sigma_010.eps");
	c6 -> SaveAs("../omg_plots/omg_sigma_010.jpg");
    }
    else if(particle == "antiomg"){
	c6 -> SaveAs("../antiomg_plots/antiomg_sigma_010.eps");
	c6 -> SaveAs("../antiomg_plots/antiomg_sigma_010.jpg");
    }

    TCanvas* c7 = new TCanvas("mass_1060");
    TGraphErrors* gr_mass_1060 = new TGraphErrors(6, x_pt_spectra, fit_mass_1060, x_pterr_spectra, fit_masserr_1060);
    gr_mass_1060 -> SetMarkerColor(2);
    gr_mass_1060 -> SetMarkerStyle(21);
    gr_mass_1060 -> Draw("AP"); 
    if(particle == "omg"){
	c7 -> SaveAs("../omg_plots/omg_mass_1060.eps");
	c7 -> SaveAs("../omg_plots/omg_mass_1060.jpg");
    }
    else if(particle == "antiomg"){
	c7 -> SaveAs("../antiomg_plots/antiomg_mass_1060.eps");
	c7 -> SaveAs("../antiomg_plots/antiomg_mass_1060.jpg");
    }

    TCanvas* c8 = new TCanvas("sigma_1060");
    TGraphErrors* gr_sigma_1060 = new TGraphErrors(6, x_pt_spectra, fit_sigma_1060, x_pterr_spectra, fit_sigmaerr_1060);
    gr_sigma_1060 -> SetMarkerColor(2);
    gr_sigma_1060 -> SetMarkerStyle(21);
    gr_sigma_1060 -> Draw("AP"); 
    gr_sigma_1060 -> SaveAs("../");
    if(particle == "omg"){
        c8 -> SaveAs("../omg_plots/omg_sigma_1060.eps");
        c8 -> SaveAs("../omg_plots/omg_sigma_1060.jpg");
    }
    else if(particle == "antiomg"){
	c8 -> SaveAs("../antiomg_plots/antiomg_sigma_1060.eps");
	c8 -> SaveAs("../antiomg_plots/antiomg_sigma_1060.jpg");
    }

    cout<<"happy"<<endl;
    //gr_010 -> SaveAs("../test.eps");
    return 0;
}
