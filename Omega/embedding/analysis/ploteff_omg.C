{
    const Int_t kPtBin = 6;
    Double_t effx[6] = {0.95, 1.4, 1.8, 2.2, 2.6, 3.2};
    //==== Flat ====
    Double_t eff_omg_flat[2][6] = {{0.00300657, 0.0184336, 0.0393808, 0.0559498, 0.0667766, 0.0789998}, {0.00276196, 0.0165463, 0.0324728, 0.0489268, 0.0594022, 0.0699429}};
    Double_t efferr_omg_flat[2][6] = {{0.0028194, 0.0064705, 0.00565115, 0.0268709, 0.00392026, 0.00388809}, {0.00287008, 0.00533756, 0.00617895, 0.00425882, 0.00421918, 0.00470884}};
    Double_t eff_omg_flat_scale[2][6] = {{0.000300657, 0.00184336, 0.00393808, 0.00559498, 0.00667766, 0.00789998}, {0.00276196, 0.0165463, 0.0324728, 0.0489268, 0.0594022, 0.0699429}};
    Double_t efferr_omg_flat_scale[2][6] = {{0.00028194, 0.00064705, 0.000565115, 0.00268709, 0.000392026, 0.000388809}, {0.00287008, 0.00533756, 0.00617895, 0.00425882, 0.00421918, 0.00470884}};

    //==== Exp. =====
    Double_t eff_omg_exp[2][6] = {{0.00346035, 0.0195844, 0.0392631, 0.057325, 0.0697926, 0.0705259}, {0.00292043, 0.0164087, 0.0326063, 0.0466934, 0.0613472, 0.0628396}};
    Double_t eff_omg_exp_scale[2][6] = {{0.000346035, 0.00195844, 0.00392631, 0.0057325, 0.00697926, 0.00705259}, {0.00292043, 0.0164087, 0.0326063, 0.0466934, 0.0613472, 0.0628396}};
    Double_t efferr_omg_exp[2][6] = {{0.00266238, 0.00789188, 0.00674196, 0.00593556, 0.0135425, 0.0247056}, {0.00209784, 0.00509533, 0.00583159, 0.00729448, 0.00667792, 0.023782}};
    Double_t efferr_omg_exp_scale[2][6] = {{0.000266238, 0.000789188, 0.000674196, 0.000593556, 0.00135425, 0.00247056}, {0.00209784, 0.00509533, 0.00583159, 0.00729448, 0.00667792, 0.023782}};

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
    ceff->SaveAs("../plots/eff_omg.png");
    ceff->SaveAs("../plots/eff_omg.eps");
    ceff->SaveAs("../plots/eff_omg.pdf");
}
