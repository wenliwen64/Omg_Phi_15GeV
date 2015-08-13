Int_t plot_antiomgrot_15GeV(){
    string particle("antiomg");
    ofstream scalerot_dat("./scale_rot_antiomg.dat");
    ofstream levy_dat("./levy_par_antiomg.dat");
    //ifstream eff_file("../embedding/analysis/eff_omg_exp.dat");
    ifstream fp_eff_file("../embedding/analysis/eff_omgbar_fp.dat");
    ifstream exp_eff_file("../embedding/analysis/eff_omgbar_exp.dat");
    ifstream spectra_xpos_file("./spectra_xpos_antiomg.dat");
    //gStyle -> SetOptFit(111);
    double pdgmass_omg = 1.67245;
    TFile* infile_rot;
    TFile* infile_dat;
    if(particle == "omg"){
	infile_rot = new TFile("0715_2015_omgrot.local_analysis.root", "read");
        infile_dat = new TFile("0628_2015_omg.local_analysis.root", "read");
    }
    else if(particle == "antiomg"){
        infile_rot = new TFile("0806_2015_antiomgrot.local_analysis.root", "read");
	infile_dat = new TFile("0806_2015_antiomg.local_analysis.root", "read");
    }
    else{
        std::cout<<"Wrong Particle Type!"<<std::endl;
        return -1;
    }

    const Int_t kPtBin = 6;
    const Float_t ptbd[kPtBin+1] = {0.7, 1.2, 1.6, 2.0, 2.4, 2.8, 3.6};
    double int_l = pdgmass_omg - 0.005;
    double int_u = pdgmass_omg + 0.005;
    double lb = 1.6225;
    double ub = 1.72;
    double sig_counts_010[6];
    double sig_counts_1060[6];
    //Double_t eff[2][6] = {{0.00346035, 0.0195844, 0.0393808, 0.0559498, 0.0667766, 0.0789998}, {0.00292043, 0.0164087, 0.0324728, 0.0489268, 0.0594022, 0.0699429}};
    Double_t eff[2][6] = {{0.00262769, 0.0167903, 0.0365496, 0.0529927, 0.0670094, 0.0748538}, {0.00232462, 0.0158278, 0.0324728, 0.0489268, 0.0594022, 0.0699429}};
    for(Int_t icent = 0; icent < 2; icent++){
        for(Int_t ipt = 0; ipt < 6; ipt++){
            Double_t dummy;
            if(ipt < 2){
		exp_eff_file >> dummy >> dummy >> eff[icent][ipt] >> dummy;
		fp_eff_file >> dummy >> dummy >> dummy >> dummy;
            }
            else{
		exp_eff_file >> dummy >> dummy >> dummy >> dummy;
		fp_eff_file >> dummy >> dummy >> eff[icent][ipt] >> dummy;
	    }
	}
    }
    fp_eff_file.close();
    exp_eff_file.close();

    TH1F* hsig_010[kPtBin];
    TH1F* hsig_1060[kPtBin];
    //TF1* f1 = new TF1("f1", "[0] * exp(-(x - [1])*(x - [1]) / (2 * [2] * [2])) + [3] * exp(-(x - [1])*(x - [1]) / (2 * [4] * [4])) + [5] + [6] * x + [7] * x * x", lb, ub);//lb, ub);
    //TF1* f1 = new TF1("f1", "[0] * exp(-(x - [1])*(x - [1]) / (2 * [2] * [2])) + [3] + [4] * x + [5] * x * x + [6] * x * x * x", lb, ub);//lb, ub);
    TF1* f_sig = new TF1("f_sig", "[0] * exp(-(x - [1])*(x - [1]) / (2 * [2] * [2]))", lb, ub);
    //TF1* f_bg = new TF1("f_bg", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", lb, ub);

    //f1 -> SetParName(0, "Yield");
    //f1 -> SetParName(1, "Mean");
    //f1 -> SetParName(2, "Sigma");
    //f1 -> SetParName(3, "Pol0");
    //f1 -> SetParName(4, "Pol1");
    //f1 -> SetParName(5, "Pol2");
    //f1 -> SetParName(6, "Pol3");
    f_sig->SetParName(0, "Yield");
    f_sig->SetParName(1, "Mean");
    f_sig->SetParName(2, "Sigma");
    f_sig->SetParameter(1, 1.671);
    f_sig->SetParameter(2, 0.006);
    for(int i = 0; i < kPtBin; i++){
        char hist_name_sig_010[200];
        char hist_name_sig_1060[200];
        char hist_name_rot_010[200];
        char hist_name_rot_1060[200];
        sprintf(hist_name_sig_010, "sig_xipt%dcent_010", i+1);
        sprintf(hist_name_sig_1060, "sig_xipt%dcent_1060", i+1);
        sprintf(hist_name_rot_010, "sig_xipt%dcent_010", i+1);
        sprintf(hist_name_rot_1060, "sig_xipt%dcent_1060", i+1);

        char can_name_sig_010[200];
        char can_name_sig_1060[200];
        TString csig_name_010;
        TString csig_name_1060;

        if(particle == "omg"){
	    sprintf(can_name_sig_010, "../omg_plots/rot_xipt%dcent_010.png", i+1);
	    sprintf(can_name_sig_1060, "../omg_plots/rot_xipt%dcent_1060.png", i+1);
            csig_name_010.Form("../omg_plots/pureomg_pt%dcent010.png", i+1);
            csig_name_1060.Form("../omg_plots/pureomg_pt%dcent1060.png", i+1);
	}
        else if(particle == "antiomg"){
	    sprintf(can_name_sig_010, "../antiomg_plots/rot_xipt%dcent_010.eps", i+1);
	    sprintf(can_name_sig_1060, "../antiomg_plots/rot_xipt%dcent_1060.eps", i+1);
	    csig_name_010.Form("../antiomg_plots/pureomg_pt%dcent010.eps", i+1);
            csig_name_1060.Form("../antiomg_plots/pureomg_pt%dcent1060.eps", i+1);

	}
   
        TH1F* hdat_010 = (TH1F*) infile_dat->Get(hist_name_sig_010);
        hdat_010->Sumw2();
        TH1F* hdat_1060 = (TH1F*) infile_dat->Get(hist_name_sig_1060);
        hdat_1060->Sumw2();
        TH1F* hrot_010 = (TH1F*) infile_rot->Get(hist_name_rot_010);
        hrot_010->Sumw2();
        TH1F* hrot_1060 = (TH1F*) infile_rot->Get(hist_name_rot_1060);
        hrot_1060->Sumw2();

        TCanvas* c_010 = new TCanvas("c_010");
        TCanvas* csig_010 = new TCanvas("csig_010");
        TCanvas* c_1060 = new TCanvas("c_1060");
        TCanvas* csig_1060 = new TCanvas("csig_1060");

        Float_t scale_ratio = 0.;
        Int_t ratio_l1 = hrot_010->FindBin(pdgmass_omg-0.05);
        Int_t ratio_l2 = hrot_010->FindBin(pdgmass_omg+0.01);
        Int_t ratio_u1 = hrot_010->FindBin(pdgmass_omg-0.01);
        Int_t ratio_u2 = hrot_010->FindBin(pdgmass_omg+0.05);
        scale_ratio = (hrot_010->Integral(ratio_l1, ratio_u1) + hrot_010->Integral(ratio_l2, ratio_u2)) / (hdat_010->Integral(ratio_l1, ratio_u1) + hdat_010->Integral(ratio_l2, ratio_u2));

        scalerot_dat << "010 " << i << " " << scale_ratio << endl;
        c_010->cd();
        hdat_010->SetMarkerStyle(8);
        hdat_010->Draw("Hist");
        hdat_010->GetXaxis()->SetTitle("InvMass(GeV)");
	hdat_010->GetYaxis()->SetTitle("counts");
    
        TH1F* hrot_010_copy = (TH1F*)hrot_010->Clone();
        hrot_010_copy->Sumw2();
        hrot_010_copy->Scale(1./scale_ratio);  
        hrot_010_copy->SetLineColor(2);
        hrot_010_copy->SetFillColor(2);
        hrot_010_copy->SetFillStyle(3354);
        gPad->SetTicks(1, 1);
        hrot_010_copy->Draw("Hist same");
        c_010 -> SaveAs(can_name_sig_010);

        csig_010->cd();
        hsig_010[i] = (TH1F*)hdat_010->Clone();
        hsig_010[i]->Sumw2();
        hsig_010[i]->Add(hrot_010_copy, -1);
        hsig_010[i]->SetMarkerStyle(8);
        hsig_010[i]->SetMarkerColor(1);
        hsig_010[i]->SetFillColor(kCyan);
        hsig_010[i]->SetFillStyle(3354);
        gPad->SetTicks(1, 1);
        hsig_010[i]->Draw("E");
        if(i == 5) f_sig->SetParameter(1, 1.6694); 
        //hsig_010[i]->Fit("f_sig", "REM");
        csig_010->SaveAs(csig_name_010.Data());
        Int_t lb_bin = hsig_010[i]->FindBin(int_l);
        Int_t ub_bin = hsig_010[i]->FindBin(int_u);
        sig_counts_010[i] = hsig_010[i]->Integral(lb_bin, ub_bin);
        cout << sig_counts_010[i] << "======sig_counts010" << endl;
        f_sig->SetParameter(1, 1.671);

	//================================================================
        c_1060->cd();
        scale_ratio = (hrot_1060->Integral(ratio_l1, ratio_u1) + hrot_1060->Integral(ratio_l2, ratio_u2)) / (hdat_1060->Integral(ratio_l1, ratio_u1) + hdat_1060->Integral(ratio_l2, ratio_u2));
        scalerot_dat << "1060 " << i << " " << scale_ratio << endl;
	hdat_1060->SetMarkerStyle(8);
        hdat_1060->Draw("Hist ");
	hdat_1060->GetXaxis()->SetTitle("InvMass(GeV)");
	hdat_1060->GetYaxis()->SetTitle("counts");

        TH1F* hrot_1060_copy = (TH1F*)hrot_1060->Clone();
        hrot_1060_copy->Sumw2();
        hrot_1060_copy->Scale(1./scale_ratio);
        hrot_1060_copy->SetLineColor(2);
        hrot_1060_copy->SetFillColor(2);
        hrot_1060_copy->SetFillStyle(3354);
        gPad->SetTicks(1, 1);
        hrot_1060_copy->Draw("Hist  same");
        c_1060->SaveAs(can_name_sig_1060); 
     
        csig_1060->cd(); 
	hsig_1060[i] = (TH1F*)hdat_1060->Clone();
        hsig_1060[i]->Sumw2();
        hsig_1060[i]->Add(hrot_1060_copy, -1);
        hsig_1060[i]->SetMarkerStyle(8);
        hsig_1060[i]->SetMarkerColor(1);
        hsig_1060[i]->Draw("E");
        gPad->SetTicks(1, 1);
        //hsig_1060[i]->Fit("f_sig", "REM");   
        sig_counts_1060[i] = hsig_1060[i]->Integral(lb_bin, ub_bin);//h_1060 -> Integral(lb_bin, ub_bin) - bg_counts;
        csig_1060->SaveAs(csig_name_1060.Data());
        cout << sig_counts_1060[i] << "=======>sigcounts_1060" << endl;
    }

//******************** Plot Spectrum ***************************** 
    TF1* levy = new TF1("levy","[0]*pow(1+(sqrt(x*x+1.67245*1.67245)-1.67245)/([1]*[2]),-[1])*([1]-1)*([1]-2)/(2*3.14159265*[1]*[2]*([1]*[2]+1.67245*([1]-2)))",0.,8.);
    levy->SetParName(0, "dN/dy");
    levy->SetParName(1, "a");
    levy->SetParName(2, "T");
    levy->SetParameter(0, 0.05);
    levy->SetParameter(1, 900000);
    levy->SetParameter(2, 0.26);

    //===== Initial Calculation ====
    double PI = 3.1415926;
    double x_pt_spectra[2][kPtBin];
    int nevents[9] = {1.538113e6, 2.45574e6, 2.623553e6, 2.703106e6, 2.69095e6, 2.725661e6, 2.68739e6, 1.293578e6, 1.333159e6};
    for(int i = 0; i < 2; i++){
        double dummy;
        spectra_xpos_file >> dummy >> x_pt_spectra[i][0] >> x_pt_spectra[i][1] >> x_pt_spectra[i][2] >> x_pt_spectra[i][3] >> x_pt_spectra[i][4] >> x_pt_spectra[i][5];
    }
    spectra_xpos_file.close();
    ofstream spectra_xpos_file_o("./spectra_xpos_antiomg.dat");

    double x_pterr_spectra[] = {0, 0, 0, 0, 0, 0};
    double dpt_spectra[6] = {0.5, 0.4, 0.4, 0.4, 0.4, 0.8 };
    double y_pt_spectra_010[6] = {};  
    double y_pt_spectra_1060[6] = {};  
    double y_pt_spectra_err_010[6] = {};  
    double y_pt_spectra_err_1060[6] = {};  

    for(int j = 0; j < 6; j++){
        //---- 1060 ----
	y_pt_spectra_1060[j] = 1/(2*PI) * sig_counts_1060[j] / x_pt_spectra[0][j] / dpt_spectra[j] / (nevents[6] + nevents[5] + nevents[4] + nevents[3] + nevents[2])/ eff[0][j]; 
        y_pt_spectra_err_1060[j] = 1/(2*PI) * sqrt(sig_counts_1060[j]) / x_pt_spectra[0][j] / dpt_spectra[j] / (nevents[6]+nevents[5] + nevents[4] + nevents[3] + nevents[2]) / eff[0][j];
        //---- 010 ----  
	y_pt_spectra_010[j] = 1/(2*PI) * sig_counts_010[j] / x_pt_spectra[1][j] / dpt_spectra[j]/(nevents[7] + nevents[8]) / eff[1][j];
	y_pt_spectra_err_010[j] = 1/(2*PI) * sqrt(sig_counts_010[j]) / x_pt_spectra[1][j] / dpt_spectra[j]/(nevents[7] + nevents[8]) / eff[1][j];
        cout<<"testtestestesttestest"<<y_pt_spectra_1060[j] <<  "<---->" << y_pt_spectra_010[j] << endl;
	printf("sig_count = %.10f<->%.10f\n", sig_counts_010[j]/eff[1][j], sig_counts_1060[j]/eff[0][j]);
    }

    //==== Define Fitting Function ====
    TF1* f1 = new TF1("f1", "[0] * exp(-(x - [1])*(x - [1]) / (2 * [2] * [2])) + [3] + [4] * x + [5] * x * x + [6] * x * x * x", lb, ub);//lb, ub);
    TF1* f_sig = new TF1("f_sig", "[0] * exp(-(x - [1])*(x - [1]) / (2 * [2] * [2]))", lb, ub);
    TF1* f_bg = new TF1("f_bg", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", lb, ub);
    TF1* levy_pt = new TF1("levypt","x*[0]*pow(1+(sqrt(x*x+1.67245*1.67245)-1.67245)/([1]*[2]),-[1])*([1]-1)*([1]-2)/(2*3.14159265*[1]*[2]*([1]*[2]+1.67245*([1]-2)))",0.,8.);
    TF1* levy_pt2 = new TF1("levypt2","x*x*[0]*pow(1+(sqrt(x*x+1.67245*1.67245)-1.67245)/([1]*[2]),-[1])*([1]-1)*([1]-2)/(2*3.14159265*[1]*[2]*([1]*[2]+1.67245*([1]-2)))",0.,8.);
    Double_t fitting_par_1060[3];
    Double_t fitting_par_010[3];
    
    //==== 1060 ====
    TGraphErrors* cur_g_1060 = new TGraphErrors(6, x_pt_spectra[0], y_pt_spectra_1060, x_pterr_spectra, y_pt_spectra_err_1060);
    cur_g_1060->Fit(levy, "REM0");
    spectra_xpos_file_o << "1060";
    //==== to get xpt corrected position ====
    for(int i = 0; i < 6; i++){
        levy_pt->SetParameter(0, levy->GetParameter(0));
        levy_pt->SetParameter(1, levy->GetParameter(1));
        levy_pt->SetParameter(2, levy->GetParameter(2));
        levy_pt2->SetParameter(0, levy->GetParameter(0));
        levy_pt2->SetParameter(1, levy->GetParameter(1));
        levy_pt2->SetParameter(2, levy->GetParameter(2));
	x_pt_spectra[0][i] = levy_pt2->Integral(ptbd[i], ptbd[i+1])/levy_pt->Integral(ptbd[i], ptbd[i+1]);
        spectra_xpos_file_o << " " << x_pt_spectra[0][i];
        cout << x_pt_spectra[0][i]<<"------------------"<<endl;
    }
    spectra_xpos_file_o << endl;

    Double_t y_pt_spectra_1060_scale[6];
    Double_t y_pt_spectra_err_1060_scale[6];

    for(int j = 0; j < 6; j++){
	y_pt_spectra_1060[j] = 1/(2*PI) * sig_counts_1060[j] / x_pt_spectra[0][j] / dpt_spectra[j] / (nevents[6]+nevents[5] + nevents[4] + nevents[3] + nevents[2]) / eff[0][j]; 
	y_pt_spectra_1060_scale[j] = 0.1*1/(2*PI) * sig_counts_1060[j] / x_pt_spectra[0][j] / dpt_spectra[j] / (nevents[6]+nevents[5] + nevents[4] + nevents[3] + nevents[2]) / eff[0][j]; 
        y_pt_spectra_err_1060[j] = 1/(2*PI) * sqrt(sig_counts_1060[j]) / x_pt_spectra[0][j] / dpt_spectra[j] / (nevents[6]+nevents[5] + nevents[4] + nevents[3] + nevents[2]) / eff[0][j];
        y_pt_spectra_err_1060_scale[j] = 0.1*1/(2*PI) * sqrt(sig_counts_1060[j]) / x_pt_spectra[0][j] / dpt_spectra[j] / (nevents[6]+nevents[5] + nevents[4] + nevents[3] + nevents[2]) / eff[0][j];
        cout << "1060 pt = " << x_pt_spectra[1][j] << " levy = " << y_pt_spectra_1060[j] << endl;
    }

    //TCanvas* cpt_omg_1060 = new TCanvas("cpt_omg_1060", "cpt_omg_1060", 200, 10, 600, 400);
    //cpt_omg_1060 -> SetLogy();
    TCanvas* spectra_can = new TCanvas("spectra_can", "spectra_can");
    spectra_can->SetLogy();
    spectra_can->SetTicks(1, 1);
    //TGraphErrors* cur_g_1060_new = new TGraphErrors(6, x_pt_spectra[0], y_pt_spectra_1060_scale, x_pterr_spectra, y_pt_spectra_err_1060_scale);
    TGraphErrors* cur_g_1060_new = new TGraphErrors(6, x_pt_spectra[0], y_pt_spectra_1060, x_pterr_spectra, y_pt_spectra_err_1060);
    cur_g_1060_new->SetMarkerSize(1.0);
    cur_g_1060_new->SetMarkerStyle(20);
    cur_g_1060_new->SetMarkerColor(2);
    cur_g_1060_new->SetMaximum(10E-2);
    cur_g_1060_new->SetMinimum(10E-8);
    cur_g_1060_new->GetXaxis()->SetLimits(0.5, 3.60);
    cur_g_1060_new->SetTitle("#Omega^{-} Spectra, Au+Au 14.6GeV");
    cur_g_1060_new->GetYaxis()->SetTitle("#frac{d^{2}N}{2#piNP_{T}dP_{T}dy}(GeV/c)^{-2}");
    cur_g_1060_new->GetXaxis()->SetTitle("P_{T}(GeV/c)");
    cur_g_1060_new->GetYaxis()->SetTitleOffset(1.3);
    if(particle == "antiomg")
	cur_g_1060_new -> SetTitle("#Omega^{+} Spectra, Au+Au 14.6GeV");
    cur_g_1060_new->Draw("AP");
    cur_g_1060_new->Fit(levy, "REM0");
    levy->GetParameters(fitting_par_1060);
    levy_dat << "1060 " << levy->GetParameter(0) << " " << levy->GetParameter(1) << " " << levy->GetParameter(2) << std::endl;
/*
    if(particle == "omg"){
	cpt_omg_1060->SaveAs("../omg_plots/omg_pt_spectra_1060_rot.eps");
    }
    else{
	cpt_omg_1060->SaveAs("../antiomg_plots/omg_pt_spectra_1060_rot.eps");
    }
*/
    //==== 010 ====
    TGraphErrors* cur_g = new TGraphErrors(6, x_pt_spectra[1], y_pt_spectra_010, x_pterr_spectra, y_pt_spectra_err_010);
    cur_g->Fit(levy, "REM0");
    spectra_xpos_file_o << "010";
    for(int i = 0; i < 6; i++){
        levy_pt->SetParameter(0, levy->GetParameter(0));
        levy_pt->SetParameter(1, levy->GetParameter(1));
        levy_pt->SetParameter(2, levy->GetParameter(2));
        levy_pt2->SetParameter(0, levy->GetParameter(0));
        levy_pt2->SetParameter(1, levy->GetParameter(1));
        levy_pt2->SetParameter(2, levy->GetParameter(2));
	x_pt_spectra[1][i] = levy_pt2->Integral(ptbd[i], ptbd[i+1])/levy_pt->Integral(ptbd[i], ptbd[i+1]);
        spectra_xpos_file_o << " " << x_pt_spectra[1][i];
        cout<<x_pt_spectra[1][i]<<"------------------"<<endl;
    }
    spectra_xpos_file_o << std::endl;
    spectra_xpos_file_o.close();
    for(int j = 0; j < 6; j++){
	y_pt_spectra_010[j] = 1/(2*PI) * sig_counts_010[j] / x_pt_spectra[1][j] / dpt_spectra[j]/(nevents[7] + nevents[8]) / eff[1][j];
	y_pt_spectra_err_010[j] = 1/(2*PI) * sqrt(sig_counts_010[j]) / x_pt_spectra[1][j] / dpt_spectra[j]/(nevents[7] + nevents[8]) / eff[1][j];
        cout << "010 pt = " << x_pt_spectra[1][j] << " levy = " << y_pt_spectra_010[j] << endl;
    }
    //TCanvas* cpt_omg_010 = new TCanvas("cpt_omg_010", "cpt_omg_010", 200, 10, 600, 400);
    //cpt_omg_010 -> SetLogy();
    TGraphErrors* cur_g_010_new = new TGraphErrors(6, x_pt_spectra[1], y_pt_spectra_010, x_pterr_spectra, y_pt_spectra_err_010);
    cur_g_010_new->SetMarkerSize(1.0);
    cur_g_010_new->SetMarkerStyle(20);
    cur_g_010_new->SetMarkerColor(1);
    cur_g_010_new->SetMaximum(10E-2);
    cur_g_010_new->SetMinimum(10E-8);
    //cur_g_010_new->GetXaxis()->SetLimits(0.5, 3.60);
    //cur_g_010_new->SetTitle("#Omega^{-} 10%@AuAu14.5GeV");
    if(particle == "antiomg")
	cur_g_010_new->SetTitle("#Omega^{+} 0-10%@AuAu14.5GeV");
    //cur_g_010_new->GetYaxis() -> SetTitle("#frac{d^{2}N}{2#piNP_{T}dP_{T}dy}(GeV/c)^{2}");
    //cur_g_010_new->GetXaxis() -> SetTitle("P_{T}(GeV/c)");
    //cur_g_010_new->GetYaxis() -> SetTitleOffset(1.3);
    cur_g_010_new->Draw("P sames");

    //cur_g_010_new->Fit(levy, "REM0");
    levy->GetParameters(fitting_par_010);
    levy_dat << "010 " << levy->GetParameter(0) << " " << levy->GetParameter(1) << " " << levy->GetParameter(2) << std::endl;

    TF1* levy1060 = (TF1*) levy->Clone();
    levy1060->SetParameters(fitting_par_1060);
    levy1060->SetLineStyle(2);
    levy1060->SetLineColor(4);
    levy1060->Draw("sames");

    TF1* levy010 = (TF1*) levy->Clone();
    levy010->SetParameters(fitting_par_010);
    levy010->SetLineStyle(2);
    levy010->SetLineColor(4);
    levy010->Draw("sames");

    TLegend* leg = new TLegend(0.6, 0.6, 0.75, 0.8);
    leg->SetBorderSize(0);
    leg->AddEntry(cur_g_010_new, "0-10%", "p");
    leg->AddEntry(cur_g_1060_new, "10-60%", "p");
    leg->AddEntry(levy010, "Levy Function", "l");
    leg->Draw("sames");
    gPad->SaveAs("../antiomg_plots/antiomg_spectra_corrected.eps");
    gPad->SaveAs("../antiomg_plots/antiomg_spectra_corrected.png");
/*
    if(particle == "omg"){
	cpt_omg_010->SaveAs("../omg_plots/omg_pt_spectra_010_rot.eps");
    }
    else{
	cpt_omg_010->SaveAs("../antiomg_plots/omg_pt_spectra_010_rot.eps");
    }
*/
    
//====Plotting the overview figures====
/*
    TFile* file_overview = new TFile("overview.histo.root", "read");
    TCanvas* c1 = new TCanvas("c1");
    TH1F* h_centbin9 = (TH1F*) file_overview -> Get("h_centbin9_after");
    h_centbin9 -> GetXaxis() -> SetTitle("centrality bin");
    h_centbin9 -> GetYaxis() -> SetTitle("counts");
    h_centbin9 -> Draw();
    c1 -> SaveAs("../omg_plots/centbin9_gaus.eps");

    TCanvas* c2 = new TCanvas("c2");
    c2 -> SetLogz();
    TH2F* h_vpr = (TH2F*) file_overview -> Get("h_vpr_after");
    h_vpr -> GetXaxis() -> SetTitle("x(cm)");
    h_vpr -> GetYaxis() -> SetTitle("y(cm)");
    h_vpr -> Draw("colorz"); 
    c2 -> SaveAs("../omg_plots/vpr_gaus.eps");

    TCanvas* c3 = new TCanvas("c3");
    TH1F* h_vpz = (TH1F*) file_overview -> Get("h_vpz_after");
    h_vpz -> GetXaxis() -> SetTitle("z(cm)");
    h_vpz -> GetYaxis() -> SetTitle("counts");
    h_vpz -> Draw();
    c3 -> SaveAs("../omg_plots/vpz_gaus.eps");
   
    TCanvas* c4 = new TCanvas("c4");
    TH1F* h_ximass = (TH1F*) infile -> Get("ximass");
    h_ximass -> SetTitle("#Omega^{-} Invariant Mass");
    if(particle == "antiomg")
	h_ximass -> SetTitle("#Omega^{+} Invariant Mass");
    h_ximass -> GetXaxis() -> SetTitle("Invariant Mass(GeV)");
    h_ximass -> GetYaxis() -> SetTitle("counts");
    h_ximass -> Draw();
    if(particle == "omg"){
	c4 -> SaveAs("../omg_plots/ximass_gaus.eps");
    }
    else{
	c4 -> SaveAs("../antiomg_plots/ximass_gaus.eps");
    }

    TCanvas* c5 = new TCanvas("mass_010");
    TGraphErrors* gr_mass_010 = new TGraphErrors(6, x_pt_spectra, fit_mass_010, x_pterr_spectra, fit_masserr_010);
    gr_mass_010 -> SetMarkerColor(2);
    gr_mass_010 -> SetMarkerStyle(21);
    gr_mass_010 -> Draw("AP"); 
    if(particle == "omg"){
	c5 -> SaveAs("../omg_plots/omg_mass_010_gaus.eps");
	c5 -> SaveAs("../omg_plots/omg_mass_010_gaus.jpg");
    }
    else if(particle == "antiomg"){
	c5 -> SaveAs("../antiomg_plots/antiomg_mass_010_gaus.eps");
	c5 -> SaveAs("../antiomg_plots/antiomg_mass_010_gaus.jpg");
    }

    TCanvas* c6 = new TCanvas("sigma_010");
    TGraphErrors* gr_sigma_010 = new TGraphErrors(6, x_pt_spectra, fit_sigma_010, x_pterr_spectra, fit_sigmaerr_010);
    gr_sigma_010 -> SetMarkerColor(2);
    gr_sigma_010 -> SetMarkerStyle(21);
    gr_sigma_010 -> Draw("AP"); 
    if(particle == "omg"){
	c6 -> SaveAs("../omg_plots/omg_sigma_010_gaus.eps");
	c6 -> SaveAs("../omg_plots/omg_sigma_010_gaus.jpg");
    }
    else if(particle == "antiomg"){
	c6 -> SaveAs("../antiomg_plots/antiomg_sigma_010_gaus.eps");
	c6 -> SaveAs("../antiomg_plots/antiomg_sigma_010_gaus.jpg");
    }

    TCanvas* c7 = new TCanvas("mass_1060");
    TGraphErrors* gr_mass_1060 = new TGraphErrors(6, x_pt_spectra, fit_mass_1060, x_pterr_spectra, fit_masserr_1060);
    gr_mass_1060 -> SetMarkerColor(2);
    gr_mass_1060 -> SetMarkerStyle(21);
    gr_mass_1060 -> Draw("AP"); 
    if(particle == "omg"){
	c7 -> SaveAs("../omg_plots/omg_mass_1060_gaus.eps");
	c7 -> SaveAs("../omg_plots/omg_mass_1060_gaus.jpg");
    }
    else if(particle == "antiomg"){
	c7 -> SaveAs("../antiomg_plots/antiomg_mass_1060_gaus.eps");
	c7 -> SaveAs("../antiomg_plots/antiomg_mass_1060_gaus.jpg");
    }

    TCanvas* c8 = new TCanvas("sigma_1060");
    TGraphErrors* gr_sigma_1060 = new TGraphErrors(6, x_pt_spectra, fit_sigma_1060, x_pterr_spectra, fit_sigmaerr_1060);
    gr_sigma_1060 -> SetMarkerColor(2);
    gr_sigma_1060 -> SetMarkerStyle(21);
    gr_sigma_1060 -> Draw("AP"); 
    gr_sigma_1060 -> SaveAs("../");
    if(particle == "omg"){
	c8 -> SaveAs("../omg_plots/omg_sigma_1060_gaus.eps");
	c8 -> SaveAs("../omg_plots/omg_sigma_1060_gaus.jpg");
    }
    else if(particle == "antiomg"){
	c8 -> SaveAs("../antiomg_plots/antiomg_sigma_1060_gaus.eps");
	c8 -> SaveAs("../antiomg_plots/antiomg_sigma_1060_gaus.jpg");
    }
*/
    levy_dat.close();
    scalerot_dat.close();
    return 0;

}
