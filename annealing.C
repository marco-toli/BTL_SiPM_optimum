


void annealing()
{
    
    TLegend * leg;
    
    float NLUMI  = 400; //lumi/year
//     float TOTDCR = 5.4*2.6; //DCR/ yearly lumi without annealing @2.0V
    float TOTDCR = 4.2*2.6; //DCR/ yearly lumi without annealing @1.5V
    
    int NYEARS   = 10;
    int NDAYS    = 365;
    int NDAYSRUN = 243;
//     int DAYSRT   = 14;
    int DAYSRT= 14*2;
//     int DAYSRT= NDAYS-NDAYSRUN;
    
        
    float inst_lumi = (float) NLUMI/NDAYSRUN ;
    
    
    
    int NCENTERS = 4;
    float fraction[NCENTERS];
    float taus[NCENTERS];
    
    fraction[0] = 0.21; taus[0] = 0.4;
    fraction[1] = 0.36; taus[1] = 13;
    fraction[2] = 0.17; taus[2] = 600;
    fraction[3] = 0.26; taus[3] = 1e4;
    
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
    
    
    float int_DCR = 0;
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
    TGraphErrors * gDCR_anneal_naive = new TGraphErrors ();
    TGraphErrors * gDCR_hard_naive = new TGraphErrors ();
    
    float int_lumi = 0;
    


    for (int iYear = 0; iYear < NYEARS; iYear++)
    {
        for (int iDay = 0; iDay < NDAYS; iDay++)
        {
            if (DAYSRT>NDAYS-NDAYSRUN) 
            {
                std::cout << "IMPOSSIBLE!!! -->  DAYSRT>NDAYS-NDAYSRUN" << std::endl;
                break;
            }
    
            
//            std::cout << "inst_lumi = " << inst_lumi << std::endl;
            if (iDay <NDAYSRUN) // in cold + irradiation
            {
                int_lumi += inst_lumi;
                int_DCR  += (float) inst_lumi*TOTDCR/NLUMI;
                
                for (int iCenter = 0; iCenter < NCENTERS; iCenter++)
                {                                        
                    DCR[iCenter]      += inst_DCR[iCenter]; 
                    DCR_0[iCenter]     = DCR[iCenter];
                }                                    
            }
            else if (iDay>=NDAYS-DAYSRT) // no irradiation + RT annealing
            {
                int anneal_day = iDay -NDAYS+DAYSRT;
                
                for (int iCenter = 0; iCenter < NCENTERS; iCenter++)
                {
                    DCR[iCenter] =  DCR_0[iCenter]*funcAnnealingRT[iCenter]->Eval(anneal_day)/funcAnnealingRT[iCenter]->Eval(0); 
//                     std::cout << "annealing day# " << anneal_day << " :: anneal_coeff = " << funcAnnealingRT[iCenter]->Eval(anneal_day)/funcAnnealingRT[iCenter]->Eval(0) << " :: total DCR[" << iCenter << "] = "  << DCR[iCenter] << std::endl;
                }
                
            }
            
            
            
            int time = iDay + iYear*NDAYS;
            
            gLumi-> SetPoint(time, time, int_lumi/100);
            
            
            float tot_DCR = 0;
            for (int iCenter = 0; iCenter < NCENTERS; iCenter++)
            {
                tot_DCR+= DCR[iCenter];
                gDCR[iCenter] -> SetPoint(time, time, DCR[iCenter]);
            }
            
            
            
            gDCR_tot -> SetPoint(time, time, tot_DCR);
            gDCR_naive -> SetPoint(time, time, int_DCR);
            gDCR_anneal_naive -> SetPoint(time, time, int_DCR*funcAnnealingRT_tot->Eval(100));
            gDCR_hard_naive -> SetPoint(time, time, int_DCR*funcAnnealingRT_tot->Eval(10000));
            
        }
        
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
    
    TFile * fileOutput = new TFile ("output_root/annealing_scenario_temp.root", "RECREATE");
//     TFile * fileOutput = new TFile ("output_root/annealing_scenario_1.root", "RECREATE");
//     TFile * fileOutput = new TFile ("output_root/annealing_scenario_2.root", "RECREATE");
    fileOutput->cd();
    
    gDCR_naive->SetName("gDCR_naive");
    gDCR_naive->Write();
    gDCR_tot->SetName("gDCR_tot");
    gDCR_tot->Write();
    gDCR_anneal_naive->SetName("gDCR_anneal_naive");
    gDCR_anneal_naive->Write();
    for (int iCenter = 0; iCenter<NCENTERS; iCenter++)
    {
        
        gDCR[iCenter]->SetName(Form("gDCR_center_%d", iCenter));
        gDCR[iCenter]->Write();
    }
    
    
    fileOutput->Close();
    
    
    
    
}
