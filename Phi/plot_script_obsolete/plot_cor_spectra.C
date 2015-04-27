//==================Plot Corrected Spectra====================
void plot_cor_spectra(){
    TFile* spectra_plot_file = new TFile("spectra_plot.root", "recreate");
    ofstream out_file("./cor_outfile");

    //=============Get Spectra Data Points====================================
    ifstream data_file("./data_point.dat");
    char line_from_data_file[250];
    int nEvents[9] = {};
    float sig_cal[9][11] = {};
    //int cent_no_bin[7] = {0, 5, 10, 20, 30, 40, 50};
    int cent_no_bin[9]={0, 10, 20, 30, 40, 50, 60, 70, 80};
    while(data_file.getline(line_from_data_file, 250)){
          int cent_bin = 0;
          int pt_bin = 0;
          int nevents = 0;
          float ypt = 0.0;
          float sig_cal_temp = 0.0;
          sscanf(line_from_data_file, "%d %d %f %f %d", &cent_bin, &pt_bin, &ypt, &sig_cal_temp, &nevents);
          std::cout << "cent_bin" << cent_bin << std::endl;
          nEvents[cent_bin-1] = nevents;
          std::cout << nevents << " events" <<std::endl;
          sig_cal[cent_bin-1][pt_bin-1] = sig_cal_temp;
    }
    
    char line_from_eff_file[250];
    float eff[6][11] = {};
    ifstream eff_file("./efficiency_19GeV_phi.dat");
    while(eff_file.getline(line_from_eff_file, 250)){
         int cent_bin_no = 0;
         int pt_bin_no = 0;
         float pt_eff = 0.0;
         float eff_temp = 0.0;
         sscanf(line_from_eff_file, "cent%d pt%d %f %f", &cent_bin_no, &pt_bin_no, &pt_eff, &eff_temp);
         eff[cent_bin_no-1][pt_bin_no-1] = eff_temp;
    }

    //=============Recalculate the spectra values============================
    int no_events_total[9];
    std::cout<<"so far so good"<<std::endl;
    float PI = 3.1415926;
    //float pt_bin[12] = {0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.3, 1.7, 2.0, 2.5, 3.5, 5.0}; 
    //int cent_no_bin[9]={0, 10, 20, 30, 40, 50, 60, 70, 80};
    float x_pt_spectra[11] = {0.45, 0.55, 0.65, 0.75, 0.9, 1.15, 1.5, 1.85, 2.25, 3.0, 4.25};
    float dpt_spectra[11] = {0.1, 0.1, 0.1, 0.1, 0.2, 0.3, 0.4, 0.3, 0.5, 1.0, 1.5};
    float y_pt_spectra[9][11] = {};
    float y_pt_spectra_raw[8][11] = {};
    float y_pt_spectra_cor[8][11] = {};

    for(int i = 0; i < 8; i++){
	if(i==0) continue;
	else if(i==1) no_events_total[i] = nEvents[0]+nEvents[1];
        else no_events_total[i] = nEvents[i];
	for(int j = 0; j < 11; j++){
            double scale_factor = TMath::Power(10, -(i-1));
            if(i==1){
		y_pt_spectra_raw[i-1][j] = scale_factor*(sig_cal[0][j]+sig_cal[1][j])/2/PI/x_pt_spectra[j]/dpt_spectra[j]/no_events_total[1];
		//y_pt_spectra_cor[i-1][j] = scale_factor*(sig_cal[0][j]+sig_cal[1][j])/2/PI/x_pt_spectra[j]/dpt_spectra[j]/no_events_total[1] / eff[0][j];
		out_file << i << " " << j+1 << " " << y_pt_spectra_cor[i-1][j] << " " << (sig_cal[0][j]+sig_cal[1][j]) << std::endl;
            }
            else{
                y_pt_spectra_raw[i-1][j] = scale_factor*(sig_cal[i][j])/2/PI/x_pt_spectra[j]/dpt_spectra[j]/no_events_total[1] ; 
                //y_pt_spectra_cor[i-1][j] = scale_factor*(sig_cal[i][j])/2/PI/x_pt_spectra[j]/dpt_spectra[j]/no_events_total[1] / eff[i-1][j]; 
		out_file << i << " " << j+1 << " " << y_pt_spectra_cor[i-1][j] << " " << sig_cal[i][j] << std::endl;
            }
        }
    }


    TCanvas* cpt_phi_raw = new TCanvas("cpt_phi_raw", "cpt_phi_raw", 200, 10, 600, 400);
    cpt_phi_raw -> SetLogy();
    TLegend* leg_raw = new TLegend(0.65, 0.65, 0.89, 0.89);
    for(int i = 0; i < 7; i++){
        std::cout<<"graph no."<<i<<std::endl;
        TGraph* cur_g = new TGraph(11, x_pt_spectra, y_pt_spectra_raw[i]);
	cur_g -> SetMarkerSize(1.5);
	cur_g -> SetMarkerStyle(20);
	cur_g -> SetMarkerColor(2+i);
	cur_g -> SetMaximum(10.0);
        //double scale_factor = TMath::Power(10, -i);
	cur_g -> SetMinimum(1E-20);
        cpt_phi_raw -> cd();
	if(i == 0){ 
	    cur_g -> Draw("AP");
	}
        else{
            cur_g -> Draw("P");
        }
        cout<<"sofarsogood"<<endl;

        char leg_text[100];
        sprintf(leg_text, "%d%%~%d%%*10e-%d", cent_no_bin[i], cent_no_bin[i+1], i);
        leg_raw -> AddEntry(cur_g, leg_text, "p");
    }
    leg_raw -> Draw("same");
    cpt_phi_raw -> SaveAs("phi_pt_spectra_raw.eps");
    cpt_phi_raw -> SaveAs("phi_pt_spectra_raw.png");
    spectra_plot_file -> cd();
    cpt_phi_raw -> Write();

    //=============================Corrected Spectra======================= 
    /*
    TCanvas* cpt_phi = new TCanvas("cpt_phi", "cpt_phi", 200, 10, 600, 400);
    cpt_phi -> SetLogy();
    TLegend* leg = new TLegend(0.65, 0.65, 0.89, 0.89);
    for(int i = 0; i < 7; i++){
        std::cout<<"graph no."<<i<<std::endl;
        TGraph* cur_g = new TGraph(11, x_pt_spectra, y_pt_spectra_cor[i]);
	cur_g -> SetMarkerSize(1.5);
	cur_g -> SetMarkerStyle(20);
	cur_g -> SetMarkerColor(2+i);
	cur_g -> SetMaximum(10.0);
        //double scale_factor = TMath::Power(10, -i);
	cur_g -> SetMinimum(1E-20);
        cpt_phi -> cd();
	if(i == 0){ 
	    cur_g -> Draw("AP");
	}
        else{
            cur_g -> Draw("P");
        }
        cout<<"sofarsogood"<<endl;

        char leg_text[100];
        sprintf(leg_text, "%d%%~%d%% * 10e-%d", cent_no_bin[i], cent_no_bin[i+1], i);
        leg -> AddEntry(cur_g, leg_text, "p");
    }
    leg -> Draw("same");
    cpt_phi -> SaveAs("phi_pt_spectra_cor.eps");
    cpt_phi -> SaveAs("phi_pt_spectra_cor.png");
   
    spectra_plot_file -> cd();
    cpt_phi -> Write();
    */

    data_file.close();
    out_file.close();
    eff_file.close(); 
}
