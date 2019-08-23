


void annealing()
{
    
    TLegend * leg;
    
    float NLUMI  = 400; //lumi/year
//     float TOTDCR = 5.4*2.6; //DCR/ yearly lumi without annealing @2.0V
    float TOTDCR = 4.2*2.6; //DCR/ yearly lumi without annealing @1.5V

    float DCR_per_FB = 35*2.6/4000.*2.5/2;    //S12572-015C at 1.5V
//     float DCR_per_FB = 70*2.6/4000.*2.5/2;    //HDR2 at 1.5V
//     float DCR_per_FB = 40*2.6/4000.;    //FBK at 1.5V
    
    
    int NYEARS = 12;
//     int NDAYS    = 365;
    int NDAYSRUN = 243;
//     int DAYSRT   = 14;
    int DAYSRT= 14*2;
//     int DAYSRT= NDAYS-NDAYSRUN;
    
    
    TGraphErrors * gLumi_vs_Time = new TGraphErrors ("lumi_vs_time_ultimate.txt");
    TGraphErrors * gIntLumi_vs_Time = new TGraphErrors ();
    
    Double_t myInstLumi;
    Double_t myMonth;
    Double_t myIntLumi = 0.;
    float convert_coeff =  1./ 1.0e24 * (30.*24.*60.*60.)*1.0e-15/4.16;
    gIntLumi_vs_Time->SetPoint(0, 0, 0);
    
    for (int iPoint = 0; iPoint<gLumi_vs_Time->GetN(); iPoint++)
    {
        gLumi_vs_Time->GetPoint(iPoint, myMonth, myInstLumi);
        myIntLumi += myInstLumi *convert_coeff;
        std::cout << "iPoint = " << iPoint << " :: myMonth = " << myMonth << " :: myInstLumi = " << myInstLumi << " :: myIntLumi = " << myIntLumi << std::endl;
        gIntLumi_vs_Time->SetPoint(iPoint, myMonth+1, myIntLumi);
        
    }
    
    TCanvas * cLumi_vs_Time = new TCanvas ("cLumi_vs_Time", "cLumi_vs_Time", 600, 500);
    cLumi_vs_Time->cd();
    gLumi_vs_Time->GetXaxis()->SetTitle("Time [months]");
    gLumi_vs_Time->GetYaxis()->SetTitle("Luminosity [cm^{-2} s^{-1}]");
    gLumi_vs_Time->SetMarkerStyle(20);
    gLumi_vs_Time->SetMarkerColor(kRed+1);
    gLumi_vs_Time->Draw("APE");
        
    
    TCanvas * cIntLumi_vs_Time = new TCanvas ("cIntLumi_vs_Time", "cIntLumi_vs_Time", 600, 500);
    cIntLumi_vs_Time->cd();
    gIntLumi_vs_Time->GetXaxis()->SetTitle("Time [months]");
    gIntLumi_vs_Time->GetYaxis()->SetTitle("Integrated luminosity [fb^{-1}]");
//     gIntLumi_vs_Time->SetMarkerStyle(20);
//     gIntLumi_vs_Time->SetMarkerColor(kRed+1);
    gIntLumi_vs_Time->Draw("ALPE");
        
    float inst_lumi = (float) NLUMI/NDAYSRUN ;
    
    
    
    int NCENTERS = 4;
    float fraction[NCENTERS];
    float taus[NCENTERS];
    
    fraction[0] = 0.21; taus[0] = 0.4;
    fraction[1] = 0.36; taus[1] = 13;
    fraction[2] = 0.17; taus[2] = 600;
    fraction[3] = 0.26; taus[3] = 1e6;
    
    TF1 * funcAnnealingRT_tot = new TF1 ("funcAnnealingRT_tot", "[0]*exp(-x/[1]) + [2]*exp(-x/[3]) + [4]*exp(-x/[5]) + [6]*exp(-x/[7])", 0.05, 10000);
    for (int iCenter = 0; iCenter< NCENTERS; iCenter++)
    {
        funcAnnealingRT_tot->SetParameter(iCenter*2,   fraction[iCenter]); //intensity 1
        funcAnnealingRT_tot->SetParameter(iCenter*2+1, taus[iCenter]);  //time constant 1
    }
    /*
    funcAnnealingRT->SetParameter(2, 0.36); //intensity 2
    funcAnnealingRT->SetParameter(3, 13); //time constant 2
    
    funcAnnealingRT->SetParameter(4, 0.17); //intensity 3
    funcAnnealingRT->SetParameter(5, 600); //time constant 3
    
    funcAnnealingRT->SetParameter(6, 0.26); //intensity 4
    funcAnnealingRT->SetParameter(7, 1e10); //time constant 4
    */
    
    TF1 * funcAnnealingRT[NCENTERS]; 
    for (int iCenter = 0; iCenter<NCENTERS; iCenter++)
    {
        funcAnnealingRT[iCenter] = new TF1 (Form("funcAnnealingRT_%d", iCenter), "[0]*exp(-x/[1])", 0.05, 10000);
        funcAnnealingRT[iCenter] -> SetParameter(0, fraction[iCenter]);
        funcAnnealingRT[iCenter] -> SetParameter(1, taus[iCenter]);
    }
    
    
    
    TGraphErrors * gLumi = new TGraphErrors ();
    TGraphErrors * gDCR_tot = new TGraphErrors ();
    TGraphErrors * gDCR_tot_vs_Lumi = new TGraphErrors ();
    
    
    float int_DCR = 0.;
    float inst_DCR[NCENTERS];
    float DCR[NCENTERS];
    float DCR_0[NCENTERS];
    
    TGraphErrors * gDCR[NCENTERS];
    for (int iCenter = 0; iCenter < NCENTERS; iCenter++)
    {
        gDCR[iCenter] = new TGraphErrors ();
        
        DCR[iCenter] = 0;
        DCR_0[iCenter] = 0;        
        inst_DCR[iCenter] = (float) inst_lumi*TOTDCR/NLUMI*fraction[iCenter];
    }
        
    
    TGraphErrors * gDCR_naive = new TGraphErrors ();
    TGraphErrors * gDCR_naive_vs_Lumi = new TGraphErrors ();
    TGraphErrors * gDCR_anneal_naive = new TGraphErrors ();
    TGraphErrors * gDCR_hard_naive = new TGraphErrors ();
    
    float int_lumi = 0;
    
    const int NTEMP = 4;
    int temp[NTEMP] = {21, 49, 60, 80};//, 106};
    
    TF1 * funcMollAnnealingT[NTEMP];
    float alphaI[NTEMP] = {1.23, 1.28, 1.26, 1.13};//, 0.};
    float tauI[NTEMP] = {1.4e4, 260, 94, 9};//, 1.};
    float alpha0[NTEMP] = {7.07, 5.36, 4.87, 4.23};//, 3.38};
    float beta[NTEMP] = {3.29, 3.11, 3.16, 2.83};//, 2.97};
    float t0[NTEMP] = {1., 1., 1., 1.};//, 1.};
    
    
    TGraphErrors *gAlphaI = new TGraphErrors ();
    TGraphErrors *gTauI = new TGraphErrors ();
    TGraphErrors *gAlpha0 = new TGraphErrors ();
    TGraphErrors *gBeta = new TGraphErrors ();
    
    for (int iTemp = 0; iTemp < NTEMP;  iTemp++)
    {
        gAlphaI->SetPoint(iTemp , temp[iTemp], alphaI[iTemp] );
        gTauI->SetPoint(iTemp , temp[iTemp], tauI[iTemp] );
        gAlpha0->SetPoint(iTemp , temp[iTemp], alpha0[iTemp] );
        gBeta->SetPoint(iTemp , temp[iTemp], beta[iTemp] );
        
    }
    
    TCanvas * cAlphaI = new TCanvas ("cAlphaI", "cAlphaI", 600, 500);
    gAlphaI->Draw("ALPE");
    gAlphaI->SetMarkerStyle(20);
    gAlphaI->GetXaxis()->SetTitle("T [#circC]");
    gAlphaI->GetYaxis()->SetTitle("#alpha_{I}");
    gAlphaI->GetXaxis()->SetTitleSize(0.05);
    gAlphaI->GetYaxis()->SetTitleSize(0.05);
    gAlphaI->GetYaxis()->SetRangeUser(0.5, 2.5);
    
    TCanvas * cTauI = new TCanvas ("cTauI", "cTauI", 600, 500);
    gTauI->Draw("ALPE");
    gTauI->SetMarkerStyle(20);
    gTauI->GetXaxis()->SetTitle("T [#circC]");
    gTauI->GetYaxis()->SetTitle("#tau_{I}");
    gTauI->GetXaxis()->SetTitleSize(0.05);
    gTauI->GetYaxis()->SetTitleSize(0.05);
//     gTauI->GetYaxis()->SetRangeUser(0.5, .5);
    
    TCanvas * cAlpha0 = new TCanvas ("cAlpha0", "cAlpha0", 600, 500);
    gAlpha0->Draw("ALPE");
    gAlpha0->SetMarkerStyle(20);
    gAlpha0->GetXaxis()->SetTitle("T [#circC]");
    gAlpha0->GetYaxis()->SetTitle("#alpha_{0}");
    gAlpha0->GetXaxis()->SetTitleSize(0.05);
    gAlpha0->GetYaxis()->SetTitleSize(0.05);
//     gAlpha0->GetYaxis()->SetRangeUser(0.5, .5);
    
    TCanvas * cBeta = new TCanvas ("cBeta", "cBeta", 600, 500);
    gBeta->Draw("ALPE");
    gBeta->SetMarkerStyle(20);
    gBeta->GetXaxis()->SetTitle("T [#circC]");
    gBeta->GetYaxis()->SetTitle("#beta");
    gBeta->GetXaxis()->SetTitleSize(0.05);
    gBeta->GetYaxis()->SetTitleSize(0.05);
    gBeta->GetYaxis()->SetRangeUser(0.5, 5);
    
    for (int iTemp = 0; iTemp<NTEMP; iTemp++)
    {
        funcMollAnnealingT[iTemp] = new TF1 ("funcMollAnnealingT", "[0]*exp(-x/[1]) + [2] - [3]*log(x/[4])", t0[iTemp], 1e6 );
        funcMollAnnealingT[iTemp]->SetParameters(alphaI[iTemp]*1.0e-17, tauI[iTemp], alpha0[iTemp]*1.0e-17, beta[iTemp]*1.0e-18, t0[iTemp]);
    }
    
    
    //35 deg
    float myTemp = 30;
    TF1 * funcMollAnnealing35 = new TF1 ("funcMollAnnealing35", "[0]*exp(-x/[1]) + [2] - [3]*log(x/[4])", t0[0], 1e6 );
    funcMollAnnealing35->SetParameters(gAlphaI->Eval(myTemp)*1.0e-17, gTauI->Eval(myTemp), gAlpha0->Eval(myTemp)*1.0e-17, gBeta->Eval(myTemp)*1.0e-18, t0[0]);
    

    TFile * fileOutput = new TFile ("output/output_root/annealing_scenarios/annealing_scenario_30.root", "RECREATE");  
    
    int NDAYS = 150*30;    //days
    int anneal_day = 0;
    
    for (int iDay = 0; iDay < NDAYS; iDay++)
    {
        //            std::cout << "inst_lumi = " << inst_lumi << std::endl;
            if (gLumi_vs_Time->Eval(int(iDay/30)) != 0) // in cold + irradiation
            {                
                int_DCR  += gLumi_vs_Time->Eval(int(iDay/30))*convert_coeff/30.*DCR_per_FB;
                
                for (int iCenter = 0; iCenter < NCENTERS; iCenter++)
                {                                        
                    DCR[iCenter]      += gLumi_vs_Time->Eval(int(iDay/30))*convert_coeff/30.*fraction[iCenter]*DCR_per_FB; 
                    DCR_0[iCenter]     = DCR[iCenter];
                }                                
                anneal_day = 0;
            }
            else  // no irradiation + RT annealing
            {
                anneal_day ++;// iDay -NDAYS+DAYSRT;
                
                for (int iCenter = 0; iCenter < NCENTERS; iCenter++)
                {
                    float annealing_coeff; 
//                     if (anneal_day>1) annealing_coeff = funcAnnealingRT[iCenter]->Eval(anneal_day)/funcAnnealingRT[iCenter]->Eval(anneal_day);        // standard 21°C annealing
//                     else              annealing_coeff = funcAnnealingRT[iCenter]->Eval(anneal_day)/funcAnnealingRT[iCenter]->Eval(anneal_day);    // 2 days of fast annealing at 50°C                    
                    
//                     DCR[iCenter] =  DCR_0[iCenter]*annealing_coeff;                
                    
//                     DCR[iCenter] =  DCR_0[iCenter]*funcAnnealingRT[iCenter]->Eval(anneal_day)/funcAnnealingRT[iCenter]->Eval(0) * funcMollAnnealingT[3]->Eval(anneal_day*24*60)/funcMollAnnealingT[0]->Eval(anneal_day*24*60); 
                    DCR[iCenter] =  DCR_0[iCenter]*funcAnnealingRT[iCenter]->Eval(anneal_day)/funcAnnealingRT[iCenter]->Eval(0) * funcMollAnnealing35->Eval(anneal_day*24*60)/funcMollAnnealingT[0]->Eval(anneal_day*24*60); 
//                     DCR[iCenter] =  DCR_0[iCenter]*funcAnnealingRT[iCenter]->Eval(anneal_day)/funcAnnealingRT[iCenter]->Eval(0) * funcMollAnnealingT[0]->Eval(anneal_day*24*60)/funcMollAnnealingT[0]->Eval(anneal_day*24*60); 
//                      DCR[iCente0r] =  DCR_0[iCenter]*funcAnnealingRT[iCenter]->Eval(anneal_day)/1.8/funcAnnealingRT[iCenter]->Eval(0); 
                    std::cout << "annealing day# " << anneal_day << " :: annealing_coeff = " << annealing_coeff << " :: total DCR[" << iCenter << "] = "  << DCR[iCenter] << std::endl;
//                     DCR_0[iCenter] = DCR[iCenter];
                }
                
            }
            
            gLumi-> SetPoint(iDay, iDay, int_lumi/100);
            
            
            float tot_DCR = 0;
            for (int iCenter = 0; iCenter < NCENTERS; iCenter++)
            {
                tot_DCR+= DCR[iCenter];
                gDCR[iCenter] -> SetPoint(iDay, iDay, DCR[iCenter]);
            }
            
            
            
            gDCR_tot -> SetPoint(iDay, iDay, tot_DCR);            
            gDCR_tot_vs_Lumi -> SetPoint(iDay, gIntLumi_vs_Time->Eval(iDay/30.), tot_DCR);
            
            gDCR_naive -> SetPoint(iDay, iDay, int_DCR);            
            gDCR_naive_vs_Lumi -> SetPoint(iDay, gIntLumi_vs_Time->Eval(iDay/30.), int_DCR);
            
            gDCR_anneal_naive -> SetPoint(iDay, iDay, int_DCR*funcAnnealingRT_tot->Eval(100));
//             gDCR_hard_naive -> SetPoint(iDay, iDay, int_DCR*funcAnnealingRT_tot->Eval(10000));
        
        
    }
    
    TCanvas * cIntLumi = new TCanvas ("cIntLumi", "cIntLumi", 600, 600);
    gLumi->Draw("ALPE");
    gLumi->GetXaxis()->SetTitle("Time [days]");
    gPad->SetGridy();
    
    
    TCanvas *cAnnealingFunction = new TCanvas ("cAnnealingFunction", "cAnnealingFunction", 600, 600);
    gPad->SetLogy();
    gPad->SetLogx();
    for (int i= 0; i<4; i++) 
    {
        funcAnnealingRT_tot->SetParName(i*2, Form("g_{%d}", i) );
        funcAnnealingRT_tot->SetParName(i*2+1, Form("#tau_{%d}", i) );
    }
    
    
        
    funcAnnealingRT_tot->SetNpx(1000000);
    funcAnnealingRT_tot->SetLineColor(kBlack);
//     funcAnnealingRT_tot->GetYaxis()->SetRangeUser(0.0001, 1);
    funcAnnealingRT_tot->SetMinimum(0.01);
    funcAnnealingRT_tot->SetMaximum(3.);
    funcAnnealingRT_tot->GetXaxis()->SetRangeUser(0.1, 1000);
    funcAnnealingRT_tot->GetXaxis()->SetLimits(0.1, 1000);
    
    funcAnnealingRT_tot->GetXaxis()->SetTitle("Time [days]");
    funcAnnealingRT_tot->GetYaxis()->SetTitle("Relative Current");
    funcAnnealingRT_tot->GetXaxis()->SetTitleSize(0.05);
    funcAnnealingRT_tot->GetYaxis()->SetTitleSize(0.05);
    
    funcAnnealingRT_tot->Draw();
    gPad->SetLogx();
    
    
    leg = new TLegend(0.3,0.7,0.7,0.88,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
 
    leg ->AddEntry(funcAnnealingRT_tot, "total", "lp");
    
    for (int iCenter = 0; iCenter < NCENTERS; iCenter++)
    {
        funcAnnealingRT[iCenter]->SetLineColor(iCenter+2);
        funcAnnealingRT[iCenter]->Draw("same");
        
        leg ->AddEntry(funcAnnealingRT[iCenter], Form("g_{%d} = %.2f, #tau_{%d} = %.1f days", iCenter, funcAnnealingRT[iCenter]->GetParameter(0), iCenter, funcAnnealingRT[iCenter]->GetParameter(1)), "lp");
    }
    leg->Draw();
    
    
    
    
    TCanvas * cIntDCR = new TCanvas ("cIntDCR", "cIntDCR", 600, 600);
    gDCR_tot->SetLineColor(kBlack);
    gDCR_tot->SetMarkerColor(kBlack);
    gDCR_tot->SetLineWidth(2);
//     gDCR_tot->GetXaxis()->SetTitle("Integrated Luminosity [fb^{-1}]");
    gDCR_tot->GetXaxis()->SetTitle("Time [days]");
    gDCR_tot->GetYaxis()->SetTitle("DCR [GHz]");
    gDCR_tot->GetXaxis()->SetTitleSize(0.05);
    gDCR_tot->GetYaxis()->SetTitleSize(0.05);
    gDCR_tot->GetYaxis()->SetRangeUser(0,160);
    gDCR_tot->Draw("same ALPE");
    
    for (int iCenter = 0; iCenter < NCENTERS; iCenter++)
    {
        gDCR[iCenter]->SetLineColor(iCenter+2);
        gDCR[iCenter]->SetMarkerColor(iCenter+2);
//         gDCR[iCenter]->Draw("same LPE");
    }

    gDCR[3]->Draw("same LPE");
    gDCR[3]->SetLineColor(kBlue+1);
    gDCR[3]->SetMarkerColor(kBlue+1);
    
    gLumi->SetLineColor(kViolet);
    gLumi->SetMarkerColor(kViolet);
//     gLumi->Draw("same LPE");
    
    gDCR_anneal_naive->SetLineStyle(7);
    gDCR_naive->Draw("same LPE");
    
    gDCR_anneal_naive->SetMarkerColor(kGreen+1);
    gDCR_anneal_naive->SetLineColor(kGreen+1);
    gDCR_anneal_naive->SetLineStyle(7);
    gDCR_anneal_naive->Draw("same LPE");
    
//     gDCR_hard_naive->SetMarkerColor(kBlue+1);
//     gDCR_hard_naive->SetLineColor(kBlue+1);
//     gDCR_hard_naive->Draw("same LPE");
    
    gPad->SetGridy();
    
        
    
    leg = new TLegend(0.3,0.7,0.7,0.88,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    
    /*
    leg ->AddEntry(gDCR_tot, "exponential recovery model at RT from HPK APDs", "lp");
    for (int iCenter = 0; iCenter<NCENTERS; iCenter++) leg ->AddEntry(gDCR[iCenter], Form("center %d: I=%.2f, #tau = %.1f days", iCenter, fraction[iCenter], taus[iCenter]), "lp");
        
    leg ->AddEntry(gDCR_naive, "constant lumi - no annealing", "lp");
    leg ->AddEntry(gDCR_anneal_naive, "constant lumi - tot annealing of 100 days", "lp");
//     leg ->AddEntry(gLumi, "int. luminosity / 100 fb-1", "lp");
//     leg ->AddEntry(gDCR_hard_naive, "constant lumi - irreducible term", "lp");
    */
    
    leg ->AddEntry(gDCR_naive, "No annealing", "lp");
    leg ->AddEntry(gDCR_tot, "2 weeks / year", "lp");    
    leg ->AddEntry(gDCR_anneal_naive, "100 days after 2x10^{14} neq/cm^{2}", "lp");
    leg ->AddEntry(gDCR[3], "Unrecoverable component", "lp");
    leg->Draw();
    
//     gPad->SetLogy();
    
    
    
    
    TCanvas * cIntDCR_vsLumi = new TCanvas ("cIntDCR_vsLumi", "cIntDCR_vsLumi", 600, 600);
    gDCR_tot_vs_Lumi->SetLineColor(kBlack);
    gDCR_tot_vs_Lumi->SetMarkerColor(kBlack);
    gDCR_tot_vs_Lumi->SetLineWidth(2);
    gDCR_tot_vs_Lumi->GetXaxis()->SetTitle("Integrated luminosity [fb^{-1}]");
    gDCR_tot_vs_Lumi->GetYaxis()->SetTitle("DCR [GHz]");
    gDCR_tot_vs_Lumi->GetXaxis()->SetTitleSize(0.05);
    gDCR_tot_vs_Lumi->GetYaxis()->SetTitleSize(0.05);
    gDCR_tot_vs_Lumi->GetYaxis()->SetRangeUser(0,160);
    
    
    gDCR_tot_vs_Lumi->Draw("ALPE");
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

  
//     TFile * fileOutput = new TFile ("output/output_root/annealing_scenarios/annealing_scenario_20.root", "RECREATE");
//          TFile * fileOutput = new TFile ("output/output_root/annealing_scenarios/annealing_scenario_30.root", "RECREATE");
//      TFile * fileOutput = new TFile ("output/output_root/annealing_scenarios/annealing_scenario_50.root", "RECREATE");
//        TFile * fileOutput = new TFile ("output/output_root/annealing_scenarios/annealing_scenario_50_short.root", "RECREATE");
//      TFile * fileOutput = new TFile ("output/output_root/annealing_scenarios/annealing_scenario_60.root", "RECREATE");

//     TFile * fileOutput = new TFile ("output/output_root/annealing_scenarios/annealing_scenario_80.root", "RECREATE");
    
//     TFile * fileOutput = new TFile ("output_root/annealing_scenario_1.root", "RECREATE");
//     TFile * fileOutput = new TFile ("output_root/annealing_scenario_2.root", "RECREATE");
    fileOutput->cd();
    
    gDCR_naive->SetName("gDCR_naive");
    gDCR_naive->Write();
    
    gDCR_naive_vs_Lumi->SetName("gDCR_naive_vs_Lumi");
    gDCR_naive_vs_Lumi->Write();
    
    gDCR_tot->SetName("gDCR_tot");
    gDCR_tot->Write();
    gDCR_tot_vs_Lumi->SetName("gDCR_tot_vs_Lumi");
    gDCR_tot_vs_Lumi->Write();
    gDCR_anneal_naive->SetName("gDCR_anneal_naive");
    gDCR_anneal_naive->Write();
    for (int iCenter = 0; iCenter<NCENTERS; iCenter++)
    {
        
        gDCR[iCenter]->SetName(Form("gDCR_center_%d", iCenter));
        gDCR[iCenter]->Write();
    }
    
    
    fileOutput->Close();
    
    
    
    
}
