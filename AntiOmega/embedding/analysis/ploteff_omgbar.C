{
    const Int_t kPtBin = 6;
    Double_t effx[6] = {0.95, 1.4, 1.8, 2.2, 2.6, 3.2};
    //==== Flat ====
    Double_t eff_omg_flat[2][6] = {{0.00262769, 0.0167903, 0.0365496, 0.0529927, 0.0670049, 0.0748538}, {0.00232462, 0.0158278, 0.0321278, 0.0465452, 0.0582144, 0.0660179}};
    Double_t efferr_omg_flat[2][6] = {{0.00270322, 0.00637537, 0.00631142, 0.00457834, 0.0037172, 0.00525776}, {0.00211458, 0.00511403, 0.00552688, 0.00509511, 0.00343662, 0.00448935}};
    Double_t eff_omg_flat_scale[2][6] = {{0.00262769, 0.0167903, 0.0365496, 0.0529927, 0.0670049, 0.0748538}, {0.000232462, 0.00158278, 0.00321278, 0.00465452, 0.00582144, 0.00660179}};
    Double_t efferr_omg_flat_scale[2][6] = {{0.00270322, 0.00637537, 0.00631142, 0.00457834, 0.0037172, 0.00525776}, {0.000211458, 0.000511403, 0.000552688, 0.000509511, 0.000343662, 0.000448935}};

    //==== Exp. =====
    Double_t eff_omg_exp[2][6] = {{0.0033918, 0.0181645, 0.037077, 0.052368, 0.072376, 0.0813053}, {0.00288703, 0.0155606, 0.032658, 0.0447753, 0.0575065, 0.0609553}};
    Double_t efferr_omg_exp[2][6] = {{0.00249349, 0.00448822, 0.00777027, 0.00658639, 0.0125241, 0.0212422}, {0.00221113, 0.00449435, 0.00510584, 0.00375908, 0.00998374, 0.0241452}};
    Double_t eff_omg_exp_scale[2][6] = {{0.0033918, 0.0181645, 0.037077, 0.052368, 0.072376, 0.0813053}, {0.000288703, 0.00155606, 0.0032658, 0.00447753, 0.00575065, 0.00609553}};
    Double_t efferr_omg_exp_scale[2][6] = {{0.00249349, 0.00448822, 0.00777027, 0.00658639, 0.0125241, 0.0212422}, {0.000221113, 0.000449435, 0.000510584, 0.000375908, 0.000998374, 0.00241452}};

    TCanvas* ceff = new TCanvas("ceff");
    ceff->SetLogy();
    ceff->SetTicks(1, 1);

    TGraphErrors* eff_gre = new TGraphErrors(kPtBin, effx, eff_omg_flat_scale[0], 0, efferr_omg_flat_scale[0]);
    eff_gre->SetTitle("#Omega^{-} Efficiency, Au+Au 14.6GeV");
    eff_gre->GetXaxis()->SetTitle("Pt(GeV/c)");
    eff_gre->GetYaxis()->SetTitle("Efficiency");
    eff_gre->SetMarkerStyle(22);
    eff_gre->SetMarkerSize(1.35);
    eff_gre->SetMarkerColor(2);
    eff_gre->SetMaximum(0.5);
    eff_gre->SetMinimum(10e-5);
    eff_gre->Draw("AP");

    TGraphErrors* eff_gre1 = new TGraphErrors(kPtBin, effx, eff_omg_flat_scale[1], 0, efferr_omg_flat_scale[1]);
    eff_gre1->SetMarkerStyle(20);
    eff_gre1->SetMarkerSize(1.35);
    eff_gre1->SetMarkerColor(2);
    eff_gre1->SetMaximum(0.5);
    eff_gre1->SetMinimum(10e-7);
    eff_gre1->Draw("P sames");

    TGraphErrors* eff_gre2 = new TGraphErrors(kPtBin, effx, eff_omg_exp_scale[0], 0, efferr_omg_exp_scale[0]);
    eff_gre2->SetMarkerStyle(22);
    eff_gre2->SetMarkerSize(1.35);
    eff_gre2->SetMarkerColor(4);
    eff_gre2->SetMaximum(0.5);
    eff_gre2->SetMinimum(10e-7);
    eff_gre2->Draw("P sames");

    TGraphErrors* eff_gre3 = new TGraphErrors(kPtBin, effx, eff_omg_exp_scale[1], 0, efferr_omg_exp_scale[1]);
    eff_gre3->SetMarkerStyle(20);
    eff_gre3->SetMarkerSize(1.35);
    eff_gre3->SetMarkerColor(4);
    eff_gre3->SetMaximum(0.5);
    eff_gre3->SetMinimum(10e-7);
    eff_gre3->Draw("P sames");

    TLegend* leg = new TLegend(0.5, 0.15, 0.7, 0.4);
    leg->SetBorderSize(0);
    leg->AddEntry(eff_gre1, "flat 0-10%", "p");
    leg->AddEntry(eff_gre, "flat 10-60%#times 10^{-1}", "p");
    leg->AddEntry(eff_gre3, "exp. 0-10%", "p");
    leg->AddEntry(eff_gre2, "exp. 10-60%#times 10^{-1}", "p");
    leg->Draw("sames"); 
    ceff->SaveAs("../plots/eff_omgbar.png");
    ceff->SaveAs("../plots/eff_omgbar.eps");
    ceff->SaveAs("../plots/eff_omgbar.pdf");
}
