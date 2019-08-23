void SiPM_OV()
{
    
    
    TLegend * leg;
    
    const int NSIPM = 3;
    
    std::string sipm_name[NSIPM];
    sipm_name[0] = "S12572-015C";
    sipm_name[1] = "HDR2-15";
    sipm_name[2] = "FBK-thin-epi-15";
    
    std::string temp[NSIPM];
    temp[0] = "21 #circC";
    temp[1] = "21 #circC";
    temp[2] = "21 #circC";
    
    float spad_size[NSIPM];    
    spad_size[0] = 15.;
    spad_size[1] = 15.;
    spad_size[2] = 15.;
    
    float RC[NSIPM];
    RC[0] = 8.;
    RC[1] = 8.;
    RC[2] = 8.;
    
    

//     float newFLUKA_fluence = 1.;
    float newFLUKA_fluence = 1.27;  //1.9*4/3/2.1
    
    
    float QE_loss[NSIPM]; //Quantum Efficiency loss after 2e14
//     QE_loss[0] = 0.05*newFLUKA_fluence;
//     QE_loss[1] = 0.05*newFLUKA_fluence;
//     QE_loss[2] = 0.01*newFLUKA_fluence;    
    
//     QE_loss[0] = 0.1;
//     QE_loss[1] = 0.1;
//     QE_loss[2] = 0.02;
    
    QE_loss[0] = 0.;
    QE_loss[1] = 0.;
    QE_loss[2] = 0.;
    
    
    //Yuri fitted functions from measurements
    // T=-30°C, fluence = 2.1e12neq/cm², annealing RT 1 week + 80 min at 60°C, 
    // S12572 and HDR2-3015 were 9 mm² SiPMs 
    std::cout << "defining input functions ..." << std::endl;
    TF1 * funcPDE_vsWL[NSIPM];
    TF1 * funcPDE[NSIPM];
    TF1 * funcCurrent[NSIPM];
    TF1 * funcGain[NSIPM];
    TF1 * funcENC[NSIPM];
    
    
    TGraphErrors * gPDE_vs_WL[NSIPM];
    TGraphErrors * gPDE_vs_WL_norm[NSIPM];
    
    Double_t pde_x_0[45] = {360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730 , 740, 750, 760, 770, 780, 790, 800};
    Double_t pde_x_1[45] = {360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730 , 740, 750, 760, 770, 780, 790, 800};
    Double_t pde_x_2[45] = {360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730 , 740, 750, 760, 770, 780, 790, 800};
    
    Double_t pde_y_0[45] = {13.5, 15.5, 17.0, 19.0, 21.5, 23.5, 24.0, 24.5, 25.2, 26, 26.2, 26.2, 26.2, 26.1, 26.05, 26.0, 25.7, 25.5, 25, 24.5, 24.0, 23.5, 22.7, 21.5, 20.0, 18.7, 17.5, 17.3, 17.0, 16.0, 15.5, 15.0, 14.5, 14.0, 13.5, 13.0, 12.5, 12.0, 11.0, 10.0, 9.0, 8.5, 8.0, 7.5, 7.0}; //S12572-3015    
    Double_t pde_y_1[45] = {13.5, 17.0, 21.5, 25.0, 30.0, 32.5, 35.0, 37.0, 38.0, 39.0, 40.0, 41.0, 41.2, 41.3, 41.5, 41.6, 41.3, 40.0, 38.0, 37.5, 36.5, 35.0, 33.7, 32.5, 31.5, 30.0, 28.5, 27.5, 26.5, 25.2, 24.0, 23.5, 22.0, 21.0, 19.5, 18.0, 17.5, 16.5, 15.5, 14.5, 13.5, 12.5, 12.0, 11.0, 9};//HDR2-3015 
    Double_t pde_y_2[45] = {18.0, 22.0, 25.5, 27.5, 31.0, 32.0, 30.5, 29.5, 28.0, 26.0, 24.5, 23.5, 22.0, 20.0, 18.0, 17.5, 17.0, 16.0, 14.5, 14.0, 13.7, 12.0, 11.7, 11.0, 10.0, 9.8, 9.5, 9.0, 8.0, 7.7, 7.5, 7.3, 7.0, 6.5, 5.5, 5.0, 4.9, 4.8, 4.7, 4.6, 4.5, 4.3, 4.0, 3.5, 3.0}; //thin-epi
    
    gPDE_vs_WL[0] = new TGraphErrors(45, pde_x_0, pde_y_0);
    gPDE_vs_WL[1] = new TGraphErrors(45, pde_x_1, pde_y_1);
    gPDE_vs_WL[2] = new TGraphErrors(45, pde_x_2, pde_y_2);
    
    
    
    //calculated
    TF1 * funcDCR[NSIPM];
    float norm[NSIPM];
    norm[0] = 1.057;
    norm[1] = 1.082;
    norm[2] = 1.1475;
    
    for (int iSiPM = 0; iSiPM<NSIPM;iSiPM++) 
    {
        double_t x,y;
        int npoints = 45;
        float maximum = gPDE_vs_WL[iSiPM]->GetHistogram()->GetMaximum();
        gPDE_vs_WL_norm[iSiPM] = new TGraphErrors();
        std::cout << "getting PDE vs wl fit, maximum: " << maximum << std::endl;
        
        for (int iPoint = 0; iPoint < npoints; iPoint++)
        {
            gPDE_vs_WL[iSiPM]->GetPoint(iPoint, x, y);
            gPDE_vs_WL_norm[iSiPM]->SetPoint(iPoint, x, y/maximum*norm[iSiPM]);
        }
        
        //vs WL
        funcPDE_vsWL[iSiPM]= new TF1 (Form("funcPDE_vsWL_%d", iSiPM),"pol6", 350, 800);
        gPDE_vs_WL_norm[iSiPM]->Fit(funcPDE_vsWL[iSiPM], "QR");
        
//         funcPDE_vsWL[iSiPM] = funcPDE_vsWL[iSiPM]/funcPDE_vsWL[iSiPM]->GetMaximum());
        
        //vs OV
        funcPDE[iSiPM]     = new TF1 (Form("funcPDE_%d", iSiPM),"[4]*([0]+ [1]*x + [2]*x*x + [3]*x*x*x)", 0, 5);
        funcCurrent[iSiPM] = new TF1 (Form("funcCurrent_%d", iSiPM),"[4]*([0]+ [1]*x + [2]*x*x + [3]*x*x*x)", 0, 5);
//         funcDCR[iSiPM]     = new TF1 (Form("funcDCR_%d", iSiPM),"[0]+ [1]*x + [2]*x*x + [3]*x*x*x", 0, 5);        
        funcGain[iSiPM]    = new TF1 (Form("funcGain_%d", iSiPM),"[0]+[1]*x", 0, 5);        
        funcENC[iSiPM]     = new TF1 (Form("funcENC_%d", iSiPM),"[0] + [1]*x + [2]*x*x", 0, 5);
    }
    
    //LYSO input emission spectrum
    TGraphErrors* gLYSO_emission = new TGraphErrors("LYSO_emission_CPI_2777.txt");
    float max_emission = gLYSO_emission->GetHistogram()->GetMaximum();
    
    TGraphErrors* gLYSO_emission_norm  = new TGraphErrors();
    for (int iPoint = 0; iPoint < gLYSO_emission->GetN(); iPoint++)
    {
        double_t x, y;    
        gLYSO_emission->GetPoint(iPoint, x, y);
        gLYSO_emission_norm->SetPoint(iPoint, x, y/max_emission*1.1);
    }
    
    //LYSO input transmittance 1 cm
    TGraphErrors* gLYSO_transmittance = new TGraphErrors("3454_LYSO_long_transmission.txt");
    max_emission = gLYSO_transmittance->GetHistogram()->GetMaximum();
    
    TGraphErrors* gLYSO_transmittance_norm  = new TGraphErrors();
    for (int iPoint = 0; iPoint < gLYSO_transmittance->GetN(); iPoint++)
    {
        double_t x, y;    
        gLYSO_transmittance->GetPoint(iPoint, x, y);
        gLYSO_transmittance_norm->SetPoint(iPoint, x, y/max_emission*1.3);
    }
    
    
//     std::string output_folder = "./output/public_optimum_OV_sqrtN/";
//     std::string output_folder = "./output/public_optimum_OV_N/";
    std::string output_folder = "./temp/";
//     std::string output_folder = "./temp/";
    std::string output_sipm_input = "./output/public_inputs/";

//     TFile * fileOutput = new TFile ("./output/output_root/bar_performance_annealing_4monthsRT.root", "RECREATE");    
    
//     TFile * fileOutput = new TFile ("./temperature_scenarios/bar_performance_annealing_4monthsRT_T-40.root", "RECREATE");    
//     TFile * fileOutput = new TFile ("./temperature_scenarios/bar_performance_annealing_4monthsRT_T-35.root", "RECREATE");    
     TFile * fileOutput = new TFile ("./temperature_scenarios/bar_performance_annealing_4monthsRT_T-30.root", "RECREATE");    
    
//         TFile * fileOutput = new TFile ("./temperature_scenarios/bar_performance_annealing_4monthsRT_T-40_safety1.5x.root", "RECREATE");    
//     TFile * fileOutput = new TFile ("./temperature_scenarios/bar_performance_annealing_4monthsRT_T-35_safety1.5x.root", "RECREATE");    
//      TFile * fileOutput = new TFile ("./temperature_scenarios/bar_performance_annealing_4monthsRT_T-30_safety1.5x.root", "RECREATE");    
    
//     TFile * fileOutput = new TFile ("./radiation_scenarios/temp_nominal_rad_20e13_noPowerLimit.root", "RECREATE");    
//     TFile * fileOutput = new TFile ("./radiation_scenarios/temp_nominal_rad_25e13_noPowerLimit.root", "RECREATE");
//     TFile * fileOutput = new TFile ("./radiation_scenarios/temp_nominal_rad_37e13_noPowerLimit.root", "RECREATE");        
    
//     TFile * fileOutput = new TFile ("./output/output_root/bar_performance_annealing_1.root", "RECREATE");
//     TFile * fileOutput = new TFile ("./output/output_root/bar_performance_annealing_2.root", "RECREATE");
//     TFile * fileOutput = new TFile ("./output/output_root/bar_performance_annealing_3.root", "RECREATE");
//     TFile * fileOutput = new TFile ("./output/output_root/bar_performance_temp.root", "RECREATE");
    
//      TFile * fileOutput = new TFile ("./output/output_root/annealing_scenarios/bar_performance_annealing_80_fullmodel.root", "RECREATE");
//     TFile * fileOutput = new TFile ("./output/output_root/annealing_scenarios/bar_performance_annealing_30.root", "RECREATE");
//     TFile * fileOutput = new TFile ("./output/output_root/annealing_scenarios/bar_performance_annealing_50.root", "RECREATE");
//     TFile * fileOutput = new TFile ("./output/output_root/annealing_scenarios/bar_performance_annealing_60.root", "RECREATE");
//     TFile * fileOutput = new TFile ("./output/output_root/annealing_scenarios/bar_performance_annealing_80.root", "RECREATE");
    
    
    string outdaq;
    const char * daqfile;

    
    TGraphErrors * gLumi_to_DCR_annealing;
    TGraphErrors * gLumi_to_DCR_noannealing;    
    
    TFile * annealingInput = new TFile ("output/output_root/annealing_scenarios/annealing_scenario_20.root", "READ");                    
    gLumi_to_DCR_annealing = (TGraphErrors*) annealingInput->Get("gDCR_tot_vs_Lumi");        
    gLumi_to_DCR_noannealing = (TGraphErrors*) annealingInput->Get("gDCR_naive_vs_Lumi");        
    
    TCanvas * cLumi_to_DCR = new TCanvas ("cLumi_to_DCR", "cLumi_to_DCR", 600, 500);
    gLumi_to_DCR_annealing->Draw("ALPE");
    gLumi_to_DCR_annealing->GetXaxis()->SetTitle("Integrated luminosity [fb^{-1}]");
    gLumi_to_DCR_annealing->GetYaxis()->SetTitle("DCR [GHz]");
    gLumi_to_DCR_annealing->GetXaxis()->SetRangeUser(0, 4000);
    gLumi_to_DCR_annealing->GetYaxis()->SetRangeUser(0, 100);
    
    
    
    TCanvas * cComparePDE_vsWL = new TCanvas ("cComparePDE_vsWL", "cComparePDE_vsWL", 600, 500);    
    funcPDE_vsWL[0]->Draw();
    funcPDE_vsWL[0]->SetTitle("");
    funcPDE_vsWL[0]->GetXaxis()->SetTitle("Wavelength [nm]");
    funcPDE_vsWL[0]->GetYaxis()->SetTitle("PDE [%]");
    funcPDE_vsWL[0]->GetXaxis()->SetTitleSize(0.05);
    funcPDE_vsWL[0]->GetYaxis()->SetTitleSize(0.05);
    funcPDE_vsWL[0]->GetXaxis()->SetRangeUser(350, 800);
    funcPDE_vsWL[0]->GetXaxis()->SetLimits(350, 800);
    funcPDE_vsWL[0]->GetYaxis()->SetRangeUser(0, 1.6);
    
    funcPDE_vsWL[0]->SetLineColor(kBlack);        
    funcPDE_vsWL[1]->Draw("same");
    funcPDE_vsWL[1]->SetLineColor(kGreen+2);    
    funcPDE_vsWL[2]->Draw("same");
    funcPDE_vsWL[2]->SetLineColor(kRed+1);    
    
    gLYSO_emission_norm->SetLineColor(kBlue+1);
//     gLYSO_emission_norm->SetLineColor(Blue+1);
    gLYSO_emission_norm->SetLineWidth(2);
    gLYSO_emission_norm->SetLineStyle(7);
    gLYSO_emission_norm->Draw("same L");
    
    gLYSO_transmittance_norm->SetLineColor(kCyan+1);
//     gLYSO_emission_norm->SetLineColor(Blue+1);
    gLYSO_transmittance_norm->SetLineWidth(2);
    gLYSO_transmittance_norm->SetLineStyle(7);
    gLYSO_transmittance_norm->Draw("same L");
    
    TLegend * legSiPMs;
    legSiPMs = new TLegend(0.15,0.7,0.55,0.88,NULL,"brNDC");
    legSiPMs->SetBorderSize(0);
    legSiPMs->SetTextFont(42);
    legSiPMs->SetTextSize(0.03);
    legSiPMs->SetLineColor(1);
    legSiPMs->SetLineStyle(1);
    legSiPMs->SetLineWidth(1);
    legSiPMs->SetFillColor(0);
 
    for (int iSiPM = 0; iSiPM< NSIPM; iSiPM++) legSiPMs ->AddEntry(funcPDE_vsWL[iSiPM], sipm_name[iSiPM].c_str(), "lp");
    legSiPMs->Draw();
    
    legSiPMs->AddEntry(gLYSO_emission_norm, "LYSO emission spectrum", "lp");
    legSiPMs->AddEntry(gLYSO_transmittance_norm, "LYSO transmittance", "lp");
    legSiPMs->Draw();
    
    outdaq = Form("%sPDE_vs_WL.png", output_sipm_input.c_str());
    daqfile = outdaq.c_str();
    cComparePDE_vsWL->cd();
    cComparePDE_vsWL->SaveAs(daqfile);
    
    //convolute PDE correction coefficient: integral of PDE / maximum PDE
    float integral_pde[NSIPM];
    float integral_lyso;
    
    for (int iSiPM = 0; iSiPM< NSIPM; iSiPM++)
    {        
        integral_pde[iSiPM]= 0.;
        integral_lyso = 0.;
        
        for (int i = 0; i<1000; i++)
        {
            float wl = 350 + i*(800-350)/1000;
            integral_pde[iSiPM] += gLYSO_emission_norm->Eval(wl)*funcPDE_vsWL[iSiPM]->Eval(wl)/funcPDE_vsWL[iSiPM]->GetMaximum();
            integral_lyso += gLYSO_emission_norm->Eval(wl);
            
            
            
        }
        
        std::cout << "effective PDE coefficient [" << iSiPM << "] = " << integral_pde[iSiPM]/integral_lyso << std::endl;
    }
    
    
    
    
//     funcPDE_vsWL[0]->SetParameters(5.86E01, -1.88, 1.15E-02, -2.70E-05, 2.78E-08, -1.06E-11, 0);                 //S12572    @3V
//     funcPDE_vsWL[1]->SetParameters(-1.10E03, 7.3, -1.7E-02, 1.72E-05, -6.40E-09, 0, 0);                       //HDR2-3015 @2.9V
//     funcPDE_vsWL[2]->SetParameters(-1.10E04, 1.14E02, -4.84E-01, 1.08E-03, -1.34E-06, 8.75E-10, -2.36E-13);   //thin-epi  @4V
    
    
    //last coefficient is the normalization to integral lyso response --> 1/PDE(410nm)*PDE(maximum)*integral_pde
    funcPDE[0]->SetParameters(-4.1509, 14.171, -2.1475, 0.1435,   1./funcPDE_vsWL[0]->Eval(410)*funcPDE_vsWL[0]->GetMaximum()*integral_pde[0]/integral_lyso);    //S12572
    funcPDE[1]->SetParameters(3.6737, 21.277, -4.8388, 0.4198,    1./funcPDE_vsWL[1]->Eval(410)*funcPDE_vsWL[1]->GetMaximum()*integral_pde[1]/integral_lyso);     //HDR2-3015
    funcPDE[2]->SetParameters(0.9960, 11.48, -0.9055, -5.684e-14, 1./funcPDE_vsWL[2]->Eval(410)*funcPDE_vsWL[2]->GetMaximum()*integral_pde[2]/integral_lyso);  //thin-epi0
    for (int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        std::cout << "max = " <<  funcPDE_vsWL[iSiPM]->GetMaximum()<< " :: value at 410 nm = " << funcPDE_vsWL[iSiPM]->Eval(410) << " :: 1./funcPDE_vsWL[0]->Eval(410)*funcPDE_vsWL[0]->GetMaximum()*integral_pde[0] = " << 1./funcPDE_vsWL[iSiPM]->Eval(410)*funcPDE_vsWL[iSiPM]->GetMaximum()*integral_pde[iSiPM]/integral_lyso << std::endl;
    }
    
        
    
        
    //last parameter is normalization to current/mm²/fb-1 in uA
    float anneal_coeff = 1.42;  //additional annealing equivalent to 100 days at RT

    
    
    // tot DCR no annealing       --> 108 GHz @1.5V  --> 0.5
    // tot DCR with 2 weeks/year  --> 54 GHz @1.5V   --> 1
    // tot DCR with 4 weeks/year  --> 50 GHz @1.5V   --> 1.08
    // tot DCR with 4 months/year --> 43 GHz @1.5V   --> 1.26
    // tot DCR permanent only component --> 27 GHz   --> 2.0
    double insitu_recovery = 1.26;  //
    
    //      double insitu_recovery = 1.42/38*26;  //additional annealing equivalent to 100 days at 30C
//     double insitu_recovery = 1.42/38*16;  //additional annealing equivalent to 100 days at 50C
//     double insitu_recovery = 1.42/38*14;  //additional annealing equivalent to 100 days at 60C
//     double insitu_recovery = 1.42/38*13;  //additional annealing equivalent to 100 days at 80C

    
    double ref_temp = -30;
    double my_temp = -30;
    
    double T_coeff[NSIPM];
    T_coeff[0] = 1.9;
    T_coeff[1] = 1.79;
    T_coeff[2] = 1.76;
    
    double temp_coefficient[NSIPM];    
    temp_coefficient[0] = pow(T_coeff[0], (ref_temp-my_temp)/10);
    temp_coefficient[1] = pow(T_coeff[1], (ref_temp-my_temp)/10);
    temp_coefficient[2] = pow(T_coeff[2], (ref_temp-my_temp)/10);
    
    TGraphErrors * gCurrent_vs_Temperature[NSIPM];
    TGraphErrors * gDCR_vs_Temperature[NSIPM];
    
    for (int iSiPM = 0; iSiPM<NSIPM; iSiPM++)
    {
        gCurrent_vs_Temperature[iSiPM] = new TGraphErrors();
        gDCR_vs_Temperature[iSiPM] = new TGraphErrors();
    }
    
    
    
    
    
    float normal_current = 200/2.1/4000./9./anneal_coeff/insitu_recovery*newFLUKA_fluence;
//     normal_current = 1;
    
    
//     float normal_current = 1.;
    std::cout << "normalization of current = " << normal_current  << std::endl;
    funcCurrent[0]->SetParameters(1.463, 5.583, 1.737e-01, 1.253, normal_current/temp_coefficient[0]); //S12572
    funcCurrent[1]->SetParameters(1.749, 9.137, 7.317, 2.089, normal_current/temp_coefficient[0]); //HDR2-3015
    funcCurrent[2]->SetParameters(7.780e-01, 4.842, 4.502, 1.463, normal_current/temp_coefficient[0]); //HDR2-3015
    
//     funcCurrent[0]->SetParameters(1.463e-06, 5.583e-06, 1.737e-07, 1.253e-06, 1); //S12572
//     funcCurrent[1]->SetParameters(1.749e-06, 9.137e-06, 7.317e-06, 2.089e-06, 1); //HDR2-3015
//     funcCurrent[2]->SetParameters(7.780e-07, 4.842e-06, 4.502e-06, 1.463e-06, 1); //HDR2-3015
    
    funcGain[0]->SetParameters(34.847*1e3, 71.5*1e3); //S12572
    funcGain[1]->SetParameters(36.619*1e3, 97.869*1e3); //HDR2-3015
    funcGain[2]->SetParameters(43.818*1e3, 104.72*1e3); //thin-epi

    funcENC[0]->SetParameters(1.01, -1.69e-02, 0.0235); //S12572
    funcENC[1]->SetParameters(1.00644, -2.19e-03, 0.00258); //HDR2-3015
    funcENC[2]->SetParameters(1.019, -5.0e-03, 3.8e-03, 1.776e-15); //thin-epi    
    
    
    
    int TOTLUMI = 6000;
    
    float sipm_area_2x2 = 4; //mm²
    float sipm_area     = 9; //mm²
    float sipm_area_TP  = 16; //mm²

    
    float MIP_rate = 2.5e6;   //MHz
//     float MIP_rate = 0;   //MHz
    


    //define inputs
    float Vbr[NSIPM];
    Vbr[0] = 63.;        //Vb at -30°C against Vb at 21°C --> 66V
    Vbr[1] = 35.8;       //                               --> 38   
    Vbr[2] = 34.2;       //                               --> 37.7     
    
//     float maxPower      = 50000;    //mW / channel    
//     float maxPower      = 50;    //mW / channel    
    float maxPower      = 70;    //mW / channel    
    int nChannels       = 500000; //total BTL number of channels
    
    float maxPower_TP   = 40;    //mW / channel
    int nChannels_TP    = 250000; //total BTL number of channels
    
    float Edep_CMS      = 4.2; // MeV average slant thickness for muons 0.8-10 GeV
    float LY            = 40000;
    float LCE           = 0.15;  
    
    float Edep_at_TB    = 2.6;
    float PDE_at_TB     = 0.37;
    float LO            = LY*LCE*Edep_CMS;        // reference light output assumed for the time resolution "sigma_phot_ref"
    float LO_red        = LY*LCE*Edep_CMS/1.4;    // assuming 3x3 mm² SiPM on 11.5x11.5 tile
    
    float sigma_phot_TB = 42.5;
    float sigma_phot_ref= sigma_phot_TB*sqrt(Edep_at_TB/Edep_CMS);
    
    float Nphe_DCR_ref  = 9000; // re
    float DCR_ref       = 20;
    float sigma_DCR_20  = 25;
    float DCR_alpha     = 1;    
    
    float sigma_clock   = 15;
    float sigma_digi    = 9.5;
    float sigma_elect   = 11.5;
    
    
    //define steps for OV scan
    
    float min_ov = 0.5;
    float max_ov = 4.0;
    float ov_step = 0.005;
    int NOVS = (int) (max_ov-min_ov) / ov_step; //660;
    
    
    float bias_ov = 1.5;    // reference constant bias
                
                
    std::cout << "drawing input plots..." << std::endl;
    
    TCanvas * cComparePDE = new TCanvas ("cComparePDE", "cComparePDE", 600, 500);
    funcPDE[0]->Draw();
    funcPDE[0]->GetXaxis()->SetTitle("V-V_{br} [V]");
    funcPDE[0]->GetYaxis()->SetTitle("PDE for LYSO emission [%]");
    funcPDE[0]->GetXaxis()->SetTitleSize(0.05);
    funcPDE[0]->GetYaxis()->SetTitleSize(0.05);
    funcPDE[0]->GetXaxis()->SetRangeUser(0, 5);
    funcPDE[0]->GetXaxis()->SetLimits(0, 5);
    funcPDE[0]->GetYaxis()->SetRangeUser(0, 50);    
    
    funcPDE[0]->SetLineColor(kBlack);
    funcPDE[1]->Draw("same");
    funcPDE[1]->SetLineColor(kGreen+2);    
    funcPDE[2]->Draw("same");
    funcPDE[2]->SetLineColor(kRed+1);  
    gPad->SetGridx();
    
    
    legSiPMs = new TLegend(0.15,0.7,0.55,0.88,NULL,"brNDC");
    legSiPMs->SetBorderSize(0);
    legSiPMs->SetTextFont(42);
    legSiPMs->SetTextSize(0.03);
    legSiPMs->SetLineColor(1);
    legSiPMs->SetLineStyle(1);
    legSiPMs->SetLineWidth(1);
    legSiPMs->SetFillColor(0);
 
    for (int iSiPM = 0; iSiPM< NSIPM; iSiPM++) legSiPMs ->AddEntry(funcPDE[iSiPM], sipm_name[iSiPM].c_str(), "lp");
    legSiPMs->Draw();
    
    outdaq = Form("%sPDE_vs_OV.png", output_sipm_input.c_str());
    daqfile = outdaq.c_str();
    cComparePDE->cd();
    cComparePDE->SaveAs(daqfile);
    
    
    TCanvas * cCompareCurrent = new TCanvas ("cCompareCurrent", "cCompareCurrent", 600, 500);
    funcCurrent[0]->Draw();
    funcCurrent[0]->SetTitle("");
    funcCurrent[0]->GetXaxis()->SetTitle("V-V_{br} [V]");
    funcCurrent[0]->GetYaxis()->SetTitle("Current / mm^{2} / fb^{-1} [#muA]");
    funcCurrent[0]->GetXaxis()->SetTitleSize(0.05);
    funcCurrent[0]->GetYaxis()->SetTitleSize(0.05);
    funcCurrent[0]->GetXaxis()->SetRangeUser(0, 5);
    funcCurrent[0]->GetXaxis()->SetLimits(0, 5);
    funcCurrent[0]->GetYaxis()->SetRangeUser(0.003, 2);    
    
    funcCurrent[0]->SetLineColor(kBlack);        
    funcCurrent[1]->Draw("same");
    funcCurrent[1]->SetLineColor(kGreen+2);    
    funcCurrent[2]->Draw("same");
    funcCurrent[2]->SetLineColor(kRed+1);    
    legSiPMs->Draw();
    gPad->SetLogy();
    gPad->SetGridx();
    
    outdaq = Form("%sCurrent_vs_OV.png", output_sipm_input.c_str());
    daqfile = outdaq.c_str();
    cCompareCurrent->cd();
    cCompareCurrent->SaveAs(daqfile);
    
    
    TCanvas * cCompareGain = new TCanvas ("cCompareGain", "cCompareGain", 600, 500);
    funcGain[0]->Draw();
    funcGain[0]->SetTitle("");
    funcGain[0]->GetXaxis()->SetTitle("V-V_{br} [V]");
    funcGain[0]->GetYaxis()->SetTitle("Gain");
    funcGain[0]->GetXaxis()->SetTitleSize(0.05);
    funcGain[0]->GetYaxis()->SetTitleSize(0.05);
    funcGain[0]->GetXaxis()->SetRangeUser(0, 5);
    funcGain[0]->GetXaxis()->SetLimits(0, 5);
    funcGain[0]->GetYaxis()->SetRangeUser(0, 700e3);    
    
    funcGain[0]->SetLineColor(kBlack);        
    funcGain[1]->Draw("same");
    funcGain[1]->SetLineColor(kGreen+2);    
    funcGain[2]->Draw("same");
    funcGain[2]->SetLineColor(kRed+1);    
    legSiPMs->Draw();
    gPad->SetGridx();
    
    outdaq = Form("%sGain_vs_OV.png", output_sipm_input.c_str());
    daqfile = outdaq.c_str();
    cCompareGain->cd();
    cCompareGain->SaveAs(daqfile);
    
    
    TCanvas * cCompareENC = new TCanvas ("cCompareENC", "cCompareENC", 600, 500);
    funcENC[0]->Draw();
    funcENC[0]->SetTitle("");
    funcENC[0]->GetXaxis()->SetTitle("V-V_{br} [V]");
    funcENC[0]->GetYaxis()->SetTitle("ENF");
    funcENC[0]->GetXaxis()->SetTitleSize(0.05);
    funcENC[0]->GetYaxis()->SetTitleSize(0.05);
    funcENC[0]->GetXaxis()->SetRangeUser(0, 5);
    funcENC[0]->GetXaxis()->SetLimits(0, 5);
    funcENC[0]->GetYaxis()->SetRangeUser(0.9, 1.5);    
    
    funcENC[0]->SetLineColor(kBlack);        
    funcENC[1]->Draw("same");
    funcENC[1]->SetLineColor(kGreen+2);    
    funcENC[2]->Draw("same");
    funcENC[2]->SetLineColor(kRed+1);    
    legSiPMs->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    
    outdaq = Form("%sENC_vs_OV.png", output_sipm_input.c_str());
    daqfile = outdaq.c_str();
    cCompareENC->cd();
    cCompareENC->SaveAs(daqfile);
    

    
    TGraphErrors * gPDE_vs_OV[NSIPM];
    TGraphErrors * gCurrent_vs_OV[NSIPM];
    TGraphErrors * gGain_vs_OV[NSIPM];   
    TGraphErrors * gDCR_vs_OV[NSIPM];               
    
    
    for (int iSiPM = 0; iSiPM<NSIPM; iSiPM++) 
    {
        gPDE_vs_OV[iSiPM]     = new TGraphErrors ();
        gCurrent_vs_OV[iSiPM] = new TGraphErrors ();
        gGain_vs_OV[iSiPM]    = new TGraphErrors ();        
        gDCR_vs_OV[iSiPM]     = new TGraphErrors ();        
    }
    
    
    //calculation of DCR vs OV        
    for(int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        for (int i = 0; i<100; i++)
        {
            float ov = (float) i*5/100.;
//             std::cout << " ov = "  << ov << std::endl;
//             gDCR_vs_OV[iSiPM]->SetPoint(i, ov, fitCurrent[iSiPM]->Eval(ov)*1e-6/fitGain[iSiPM]->Eval(ov)/1.6e-19/1e9 * 4000 * 9);
            gDCR_vs_OV[iSiPM]->SetPoint(i, ov, funcCurrent[iSiPM]->Eval(ov)*1.e-6/funcGain[iSiPM]->Eval(ov)/1.6e-19/1e9 * 4000 * 9);
//             gDCR_vs_OV[iSiPM]->SetPoint(i, ov, funcCurrent[iSiPM]->Eval(ov)*1.e-6/funcGain[iSiPM]->Eval(ov)/1.6e-19/1e9);
        }
    }
    
    TCanvas * cDCR_vs_OV = new TCanvas ("cDCR_vs_OV", "cDCR_vs_OV", 600, 500);
    
    gDCR_vs_OV[0]->Draw("AL");    
//     gDCR_vs_OV[0]->SetTitle(sipm_name[0].c_str());
    gDCR_vs_OV[0]->GetXaxis()->SetTitle("bias OV [V]");
    gDCR_vs_OV[0]->GetYaxis()->SetTitle("DCR [GHz]");
    gDCR_vs_OV[0]->GetXaxis()->SetTitleSize(0.05);
    gDCR_vs_OV[0]->GetYaxis()->SetTitleSize(0.05);
    gDCR_vs_OV[0]->GetXaxis()->SetRangeUser(0, 8);
    gDCR_vs_OV[0]->GetXaxis()->SetLimits(0, 5);
    gDCR_vs_OV[0]->GetYaxis()->SetRangeUser(0, 500);    
    gDCR_vs_OV[0]->SetMarkerStyle(20);    

    gDCR_vs_OV[0]->SetLineColor(kBlack);    
    gDCR_vs_OV[0]->SetLineWidth(2);
    gDCR_vs_OV[1]->SetLineWidth(2);
    gDCR_vs_OV[1]->SetLineColor(kGreen+2);    
    gDCR_vs_OV[1]->Draw("same L");
    gDCR_vs_OV[2]->SetLineWidth(2);
    gDCR_vs_OV[2]->SetLineColor(kRed+1);    
    gDCR_vs_OV[2]->Draw("same L");
    
    legSiPMs->Draw();
    
    outdaq = Form("%sDCR_vs_OV.png", output_sipm_input.c_str());
    daqfile = outdaq.c_str();
    cDCR_vs_OV->cd();
    cDCR_vs_OV->SaveAs(daqfile);
    
    
    //and now extrapolations
    TGraphErrors * gPDE_vs_Lumi[NSIPM];
    TGraphErrors * gCurrent_vs_Lumi[NSIPM];
    TGraphErrors * gGain_vs_Lumi[NSIPM];    
    TGraphErrors * gDCR_vs_Lumi[NSIPM];
    TGraphErrors * gPower_vs_Lumi[NSIPM];
    TGraphErrors * gPowerDynamic_vs_Lumi[NSIPM];
    TGraphErrors * gPowerTot_vs_Lumi[NSIPM];
    TGraphErrors * gOV_vs_Lumi[NSIPM];
    TGraphErrors * gCTR_vs_Lumi[NSIPM];
    TGraphErrors * gSignal_vs_Lumi[NSIPM];
    TGraphErrors * gEleSignal_vs_Lumi[NSIPM];
    TGraphErrors * gBusyCells_vs_Lumi[NSIPM];
    
    TGraphErrors * gOptimalPDE_vs_Lumi[NSIPM];
    TGraphErrors * gOptimalDCR_vs_Lumi[NSIPM];
    TGraphErrors * gOptimalGain_vs_Lumi[NSIPM];    
    TGraphErrors * gOptimalOV_vs_Lumi[NSIPM];
    TGraphErrors * gOptimalSNR_vs_Lumi[NSIPM];
    TGraphErrors * gOptimalCTR_vs_Lumi[NSIPM];
    TGraphErrors * gOptimalCTR2x2_vs_Lumi[NSIPM];
    TGraphErrors * gOptimalCTR_Double_vs_Lumi[NSIPM];
    TGraphErrors * gOptimalCTR_Double2x2_vs_Lumi[NSIPM];
    TGraphErrors * gOptimalCTR_TP_vs_Lumi[NSIPM];
    TGraphErrors * gOptimalCTR_TP3x3_vs_Lumi[NSIPM];
    TGraphErrors * gOptimalCTR_TP2x2_vs_Lumi[NSIPM];
    TGraphErrors * gOptimalCurrent_vs_Lumi[NSIPM];
    TGraphErrors * gOptimalPower_vs_Lumi[NSIPM];
    TGraphErrors * gOptimalPowerDynamic_vs_Lumi[NSIPM];
    TGraphErrors * gOptimalPowerTot_vs_Lumi[NSIPM];
    TGraphErrors * gOptimalSignal_vs_Lumi[NSIPM];
    TGraphErrors * gOptimalEleSignal_vs_Lumi[NSIPM];
    TGraphErrors * gOptimalBusyCells_vs_Lumi[NSIPM];
    
    TGraphErrors * gOptimalPhotJitter_vs_Lumi[NSIPM];
    TGraphErrors * gOptimalDCRJitter_vs_Lumi[NSIPM];
    TGraphErrors * gOptimalPhotJitter_vs_Lumi_double[NSIPM];
    TGraphErrors * gOptimalDCRJitter_vs_Lumi_double[NSIPM];
    
    for (int iSiPM = 0; iSiPM<NSIPM; iSiPM++) 
    {
        gPDE_vs_Lumi[iSiPM]     = new TGraphErrors ();
        gCurrent_vs_Lumi[iSiPM] = new TGraphErrors ();
        gGain_vs_Lumi[iSiPM]    = new TGraphErrors ();
        gDCR_vs_Lumi[iSiPM]     = new TGraphErrors ();
        gPower_vs_Lumi[iSiPM]   = new TGraphErrors ();
        gPowerDynamic_vs_Lumi[iSiPM]   = new TGraphErrors ();
        gPowerTot_vs_Lumi[iSiPM]   = new TGraphErrors ();
        gCTR_vs_Lumi[iSiPM]   = new TGraphErrors ();
        gOV_vs_Lumi[iSiPM]    = new TGraphErrors ();
        gSignal_vs_Lumi[iSiPM]   = new TGraphErrors ();
        gEleSignal_vs_Lumi[iSiPM]   = new TGraphErrors ();
        gBusyCells_vs_Lumi[iSiPM]   = new TGraphErrors ();
        
        
        gOptimalPDE_vs_Lumi[iSiPM]   = new TGraphErrors ();
        gOptimalDCR_vs_Lumi[iSiPM]   = new TGraphErrors ();
        gOptimalGain_vs_Lumi[iSiPM]    = new TGraphErrors ();
        gOptimalOV_vs_Lumi[iSiPM]    = new TGraphErrors ();
        gOptimalSNR_vs_Lumi[iSiPM]   = new TGraphErrors ();
        gOptimalCTR_vs_Lumi[iSiPM]   = new TGraphErrors ();
        gOptimalCTR2x2_vs_Lumi[iSiPM]   = new TGraphErrors ();
        gOptimalCTR_Double_vs_Lumi[iSiPM]   = new TGraphErrors ();
        gOptimalCTR_Double2x2_vs_Lumi[iSiPM]   = new TGraphErrors ();
        gOptimalCTR_TP_vs_Lumi[iSiPM]   = new TGraphErrors ();
        gOptimalCTR_TP3x3_vs_Lumi[iSiPM]   = new TGraphErrors ();
        gOptimalCTR_TP2x2_vs_Lumi[iSiPM]   = new TGraphErrors ();
        gOptimalCurrent_vs_Lumi[iSiPM]   = new TGraphErrors ();
        gOptimalPower_vs_Lumi[iSiPM]    = new TGraphErrors ();
        gOptimalPowerDynamic_vs_Lumi[iSiPM]    = new TGraphErrors ();
        gOptimalPowerTot_vs_Lumi[iSiPM]    = new TGraphErrors ();
        gOptimalSignal_vs_Lumi[iSiPM]   = new TGraphErrors ();
        gOptimalEleSignal_vs_Lumi[iSiPM] = new TGraphErrors ();
        gOptimalBusyCells_vs_Lumi[iSiPM] = new TGraphErrors ();
        
        gOptimalPhotJitter_vs_Lumi[iSiPM] = new TGraphErrors ();
        gOptimalDCRJitter_vs_Lumi[iSiPM]  = new TGraphErrors ();
        
        gOptimalPhotJitter_vs_Lumi_double[iSiPM] = new TGraphErrors ();
        gOptimalDCRJitter_vs_Lumi_double[iSiPM]  = new TGraphErrors ();
    }
    
    TGraphAsymmErrors * gOptimalCTR_TP_vs_Lumi_Ave     = new TGraphAsymmErrors ();
    TGraphAsymmErrors * gOptimalCTR_TP3x3_vs_Lumi_Ave     = new TGraphAsymmErrors ();
    TGraphAsymmErrors * gOptimalCTR_TP2x2_vs_Lumi_Ave     = new TGraphAsymmErrors ();
    TGraphAsymmErrors * gOptimalCTR_vs_Lumi_Ave        = new TGraphAsymmErrors ();
    TGraphAsymmErrors * gOptimalCTR_2x2_vs_Lumi_Ave        = new TGraphAsymmErrors ();
    TGraphAsymmErrors * gOptimalCTR_Double_vs_Lumi_Ave = new TGraphAsymmErrors ();
    TGraphAsymmErrors * gOptimalCTR_Double2x2_vs_Lumi_Ave = new TGraphAsymmErrors ();
    
    
    
    TGraphErrors * gDigitization = new TGraphErrors ();
    TGraphErrors * gClock        = new TGraphErrors ();
    TGraphErrors * gElectronics  = new TGraphErrors ();
    
    TGraphErrors * gDigitization_double = new TGraphErrors ();
    TGraphErrors * gClock_double        = new TGraphErrors ();
    TGraphErrors * gElectronics_double  = new TGraphErrors ();
    
    
    for (int iLumi = 0; iLumi < TOTLUMI; iLumi++)
    {
    
        float min_CTR_TP4x4 = 99999;
        float max_CTR_TP4x4 = -9999;
        
        float min_CTR_TP3x3 = 99999;
        float max_CTR_TP3x3 = -9999;
        
        float min_CTR_TP2x2 = 99999;
        float max_CTR_TP2x2 = -9999;
//         float ave_CTR_TP = 0;
        
        float min_CTR_Single = 99999;
        float max_CTR_Single = -9999;
        
        float min_CTR_2x2 = 99999;
        float max_CTR_2x2 = -9999;
//         float ave_CTR_Single = 0;
        
        float min_CTR_Double = 99999;
        float max_CTR_Double = -9999;
        
        float min_CTR_Double2x2 = 99999;
        float max_CTR_Double2x2 = -9999;
//         float ave_CTR_Double = 0;
        
        gDigitization->SetPoint(iLumi, iLumi, sigma_digi);
        gClock       ->SetPoint(iLumi, iLumi, sigma_clock);
        gElectronics ->SetPoint(iLumi, iLumi, sigma_elect);
        
        gDigitization_double->SetPoint(iLumi, iLumi, sigma_digi/sqrt(2));
        gClock_double       ->SetPoint(iLumi, iLumi, sigma_clock);
        gElectronics_double ->SetPoint(iLumi, iLumi, sigma_elect/sqrt(2));
        
        for (int iSiPM = 0; iSiPM< NSIPM; iSiPM++)
        {
            
            //to include annealing model and HL-LHC scenario
//             float current = funcCurrent[iSiPM]->Eval(bias_ov) * gLumi_to_DCR_annealing->Eval(iLumi) / gLumi_to_DCR_noannealing->Eval(iLumi) * (iLumi*sipm_area *2.6);            //in uA
            
            //for fixed DCR increase
            float current = iLumi*funcCurrent[iSiPM]->Eval(bias_ov)*sipm_area;            //in uA --
            float DCR     = current/1e6 / (1.6e-19) / funcGain[iSiPM]->Eval(bias_ov)/1e9;  //current from uA to A
            
            double busy_cells = DCR/sipm_area*RC[iSiPM]*2./pow(1000./spad_size[iSiPM],2);
            if (busy_cells >1) busy_cells = 1.;
                
            double signal = funcPDE[iSiPM]->Eval(bias_ov)/100*(1.-busy_cells/1.) * LO * (1-QE_loss[iSiPM]/4000*iLumi);
            double sigma_phot = sigma_phot_ref  *sqrt(LO*PDE_at_TB/signal); 
            double sigma_DCR = sigma_DCR_20*sqrt(DCR/DCR_ref)*pow(Nphe_DCR_ref/signal, DCR_alpha);
            
            double CTR = sqrt(pow(sigma_phot,2) + pow(sigma_DCR,2)  + pow(sigma_digi,2) +  pow(sigma_elect,2) +  pow(sigma_clock,2));
            
            float power = current*(bias_ov+Vbr[iSiPM])/1e3; //converting current from uA to mA
            float tot_power = power*nChannels;
            
            double ele_signal = signal*funcGain[iSiPM]->Eval(bias_ov);
            float power_dyn   = MIP_rate*ele_signal*1.602*1e-19*(bias_ov+Vbr[iSiPM])*1e3; //converting current from Ampere to mA
            
            gCurrent_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, current);
            gGain_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, funcGain[iSiPM]->Eval(bias_ov));
            gPower_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, power);
            gPowerDynamic_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, power_dyn);
            gPowerTot_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, power_dyn+power);
            gDCR_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, DCR);
            
            gPDE_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, funcPDE[iSiPM]->Eval(bias_ov)/100*(1-QE_loss[iSiPM]/4000*iLumi));
            gOV_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, bias_ov);                        
            gCTR_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, CTR);
            
            gSignal_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, signal);
            gEleSignal_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, signal*funcGain[iSiPM]->Eval(bias_ov));
            gBusyCells_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, busy_cells);
            

            
            float best_SNR = -9999999;
            
            float best_sigma_phot = 99999;
            float best_sigma_DCR = 99999;
            
            float best_CTR  = 99999;
            float best_CTR_2x2  = 99999;
            float best_CTR_double  = 99999;
            float best_CTR_double2x2  = 99999;
            
            float best_CTR_TP4x4  = 99999;
            float best_CTR_TP3x3  = 99999;
            float best_CTR_TP2x2  = 99999;
            float best_signal = 0;
            float best_ele_signal = 0;
            float best_busy_cell = 0;
            float opt_current;
            float opt_power;
            float opt_dyn_power;
            float opt_ov;
            float opt_PDE;
            float opt_DCR;
            
            
            for (int iOV = 0; iOV<NOVS; iOV++)
            {
//                 std::cout << "testing OV = " << iOV*ov_step << std::endl;
                //calculate SNR at this OV and this lumi 
                float temp_ov      = min_ov+iOV*ov_step;
                float temp_PDE     = funcPDE[iSiPM]->Eval(temp_ov)/100 * (1-QE_loss[iSiPM]/4000*iLumi);
                
                //4x4 mm² SiPM
                float temp_current_TP4x4 = iLumi*funcCurrent[iSiPM]->Eval(temp_ov)*sipm_area_TP;                
                float temp_DCR_TP4x4     = temp_current_TP4x4/1e6 / (1.6e-19) / funcGain[iSiPM]->Eval(temp_ov)/1e9;                                
                float temp_power_TP4x4   = temp_current_TP4x4*(temp_ov+Vbr[iSiPM])/1e3;
                
                //3x3 mm² SiPM
                //to include annealing model and HL-LHC scenario
//                 float temp_current = funcCurrent[iSiPM]->Eval(temp_ov) * gLumi_to_DCR_annealing->Eval(iLumi)  / gLumi_to_DCR_noannealing->Eval(iLumi)* (iLumi*sipm_area *2.6);            //in uA
            
                //for fixed DCR increase            
                float temp_current     = iLumi*funcCurrent[iSiPM]->Eval(temp_ov)*sipm_area;                
                float temp_DCR         = temp_current/1e6 / (1.6e-19) / funcGain[iSiPM]->Eval(temp_ov)/1e9;                                
                float temp_power       = temp_current*(temp_ov+Vbr[iSiPM])/1e3; //converting current to Ampere
                
                
                
                //2x2 mm ² SiPM
                float temp_current_2x2 = iLumi*funcCurrent[iSiPM]->Eval(temp_ov)*sipm_area_2x2;                
                float temp_DCR_2x2     = temp_current_2x2/1e6 / (1.6e-19) / funcGain[iSiPM]->Eval(temp_ov)/1e9;                                
                float temp_power_2x2   = temp_current_2x2*(temp_ov+Vbr[iSiPM])/1e3;
                
                //estimating CTR

                
                double temp_busy_cells = temp_DCR/sipm_area*RC[iSiPM]*2./pow(1000./spad_size[iSiPM],2);
                if (temp_busy_cells >1) busy_cells = 1.;
                
                double temp_signal          = temp_PDE*(1.-temp_busy_cells/1.) * LO;
                double temp_signal_3x3      = temp_PDE*(1.-temp_busy_cells/1.) * LO_red;
                double temp_signal_TP2x2    = temp_PDE*(1.-temp_busy_cells/1.) * LO_red/1.5;
                double temp_signal_2x2      = temp_PDE*(1.-temp_busy_cells/1.) * LO/1.5;
                double temp_signal_2x2_bar  = temp_PDE*(1.-temp_busy_cells/1.) * LO/2.0;
                
                double temp_ele_signal = temp_signal*funcGain[iSiPM]->Eval(temp_ov);
                
                float temp_dyn_power   = MIP_rate*temp_ele_signal*1.6e-19*(bias_ov+Vbr[iSiPM])*1e3; //converting current to Ampere
                
//                 float temp_SNR     = temp_signal / sqrt(temp_DCR);
                double sigma_phot           = sigma_phot_ref*sqrt(LO*PDE_at_TB/temp_signal);                 
                double sigma_phot_3x3       = sigma_phot_ref*sqrt(LO*PDE_at_TB/temp_signal_3x3);             
                double sigma_phot_2x2       = sigma_phot_ref*sqrt(LO*PDE_at_TB/temp_signal_2x2);             
                double sigma_phot_TP2x2     = sigma_phot_ref*sqrt(LO*PDE_at_TB/temp_signal_TP2x2);           
                double sigma_phot_2x2_bar   = sigma_phot_ref*sqrt(LO*PDE_at_TB/temp_signal_2x2_bar);         

                double sigma_DCR_TP4x4   = sigma_DCR_20*sqrt(temp_DCR_TP4x4/DCR_ref)*pow(Nphe_DCR_ref/temp_signal, DCR_alpha);
                double sigma_DCR_TP3x3   = sigma_DCR_20*sqrt(temp_DCR/DCR_ref)*pow(Nphe_DCR_ref/temp_signal_3x3, DCR_alpha);
                double sigma_DCR_TP2x2   = sigma_DCR_20*sqrt(temp_DCR_2x2/DCR_ref)*pow(Nphe_DCR_ref/temp_signal_TP2x2, DCR_alpha);
                double sigma_DCR         = sigma_DCR_20*sqrt(temp_DCR/DCR_ref)*pow(Nphe_DCR_ref/temp_signal, DCR_alpha);
                double sigma_DCR_2x2     = sigma_DCR_20*sqrt(temp_DCR_2x2/DCR_ref)*pow(Nphe_DCR_ref/temp_signal_2x2, DCR_alpha);
                double sigma_DCR_2x2_bar = sigma_DCR_20*sqrt(temp_DCR_2x2/DCR_ref)*pow(Nphe_DCR_ref/temp_signal_2x2_bar, DCR_alpha);
                
                double temp_CTR_TP4x4       = sqrt(pow(sigma_phot,2) + pow(sigma_DCR_TP4x4,2)  + pow(sigma_digi,2) +  pow(sigma_elect,2) +  pow(sigma_clock,2));
                double temp_CTR_TP3x3       = sqrt(pow(sigma_phot_3x3,2) + pow(sigma_DCR,2)  + pow(sigma_digi,2) +  pow(sigma_elect,2) +  pow(sigma_clock,2));
                double temp_CTR_TP2x2       = sqrt(pow(sigma_phot_TP2x2,2) + pow(sigma_DCR_TP2x2,2)  + pow(sigma_digi,2) +  pow(sigma_elect,2) +  pow(sigma_clock,2));
                double temp_CTR             = sqrt(pow(sigma_phot,2) + pow(sigma_DCR,2)  + pow(sigma_digi,2) +  pow(sigma_elect,2) +  pow(sigma_clock,2));
                double temp_CTR_2x2         = sqrt(pow(sigma_phot_2x2,2) + pow(sigma_DCR_2x2,2)  + pow(sigma_digi,2) +  pow(sigma_elect,2) +  pow(sigma_clock,2));
                double temp_CTR_2x2_bar     = sqrt(pow(sigma_phot_2x2_bar/sqrt(2),2) + pow(sigma_DCR_2x2_bar/sqrt(2),2) + pow(sigma_digi/sqrt(2),2) +  pow(sigma_elect/sqrt(2),2) +  pow(sigma_clock,2) );
                
                
                if (true
                    && temp_CTR < best_CTR                    
//                     && temp_power<maxPower
                    && (temp_power+temp_dyn_power)<maxPower
                ) 
                {
//                     best_SNR = temp_SNR;
                    best_CTR = temp_CTR;
                    best_CTR_double = sqrt(pow(sigma_phot/sqrt(2),2) + pow(sigma_DCR/sqrt(2),2) + pow(sigma_digi/sqrt(2),2) +  pow(sigma_elect/sqrt(2),2) +  pow(sigma_clock,2) );
                    
                    
                    opt_current = temp_current;
                    opt_power = temp_power;
                    opt_dyn_power = temp_dyn_power;
                    opt_ov    = temp_ov;
                    opt_PDE    = temp_PDE;
                    opt_DCR    = temp_DCR;
                    best_signal = temp_signal;
                    best_ele_signal = temp_ele_signal;
                    best_busy_cell = temp_busy_cells;
                    
                    best_sigma_phot = sigma_phot;
                    best_sigma_DCR = sigma_DCR;                                                            
                }
                
                if (true
                    && temp_CTR_TP4x4 < best_CTR_TP4x4 
                    && temp_power_TP4x4 < maxPower_TP
                ) 
                {
                    best_CTR_TP4x4 = temp_CTR_TP4x4;
                }
                
                if (true
                    && temp_CTR_TP3x3 < best_CTR_TP3x3 
                    && temp_power < maxPower
                ) 
                {
                    best_CTR_TP2x2 = temp_CTR_TP2x2;
                }
                if (true
                    && temp_CTR_TP2x2 < best_CTR_TP2x2
                    && temp_power_2x2 < maxPower
                ) 
                {
                    best_CTR_TP3x3 = temp_CTR_TP3x3;
                }
                if (true
                    && temp_CTR_2x2 < best_CTR_2x2 
                    && temp_power_2x2 < maxPower
                ) 
                {
                    best_CTR_2x2 = temp_CTR_2x2;
                }
                if (true
                    && temp_CTR_2x2_bar < best_CTR_double2x2 
                    && temp_power_2x2 < maxPower
                ) 
                {
                    best_CTR_double2x2 = temp_CTR_2x2_bar;              
                }
                
            }
            
            if (best_CTR_TP4x4 < min_CTR_TP4x4)         min_CTR_TP4x4     = best_CTR_TP4x4;
            if (best_CTR_TP3x3 < min_CTR_TP3x3)         min_CTR_TP3x3  = best_CTR_TP3x3;
            if (best_CTR_TP2x2 < min_CTR_TP2x2)         min_CTR_TP2x2  = best_CTR_TP2x2;
            if (best_CTR < min_CTR_Single)              min_CTR_Single = best_CTR;
            if (best_CTR_2x2 < min_CTR_2x2)             min_CTR_2x2 = best_CTR_2x2;
            if (best_CTR_double < min_CTR_Double)       min_CTR_Double = best_CTR_double;
            if (best_CTR_double2x2 < min_CTR_Double2x2) min_CTR_Double2x2 = best_CTR_double2x2;
            
            if (best_CTR_TP4x4 > max_CTR_TP4x4)         max_CTR_TP4x4  = best_CTR_TP4x4;
            if (best_CTR_TP3x3 > max_CTR_TP3x3)         max_CTR_TP3x3  = best_CTR_TP3x3;
            if (best_CTR_TP2x2 > max_CTR_TP2x2)         max_CTR_TP2x2  = best_CTR_TP2x2;
            if (best_CTR > max_CTR_Single)              max_CTR_Single = best_CTR;
            if (best_CTR_2x2 > max_CTR_2x2)             max_CTR_2x2 = best_CTR_2x2;
            if (best_CTR_double > max_CTR_Double)       max_CTR_Double = best_CTR_double;
            if (best_CTR_double2x2 > max_CTR_Double2x2) max_CTR_Double2x2 = best_CTR_double2x2;
            
//             ave_CTR_TP+=best_CTR_TP
            
//             std::cout << "min_CTR_TP = " << min_CTR_TP << " :: min_CTR_Single =  " << min_CTR_Single << " :: min_CTR_Double = " << min_CTR_Double << std::endl;
//             std::cout << "max_CTR_TP = " << max_CTR_TP << " :: max_CTR_Single =  " << max_CTR_Single << " :: max_CTR_Double = " << max_CTR_Double << std::endl;
            
            
            gOptimalPDE_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, opt_PDE);
            gOptimalDCR_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, opt_DCR);
            gOptimalOV_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, opt_ov);
//             gOptimalSNR_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, best_SNR);
            gOptimalGain_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, funcGain[iSiPM]->Eval(opt_ov));
            gOptimalCTR_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, best_CTR);
            gOptimalCTR_Double_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, best_CTR_double);
            gOptimalCTR_TP_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, best_CTR_TP4x4);
            gOptimalCTR_TP3x3_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, best_CTR_TP3x3);
            gOptimalCTR_TP2x2_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, best_CTR_TP2x2);
            gOptimalCurrent_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, opt_current);
            gOptimalPower_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, opt_power);
            gOptimalPowerDynamic_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, opt_dyn_power);
            gOptimalPowerTot_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, opt_dyn_power+opt_power);
            gOptimalSignal_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, best_signal);
            gOptimalEleSignal_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, best_ele_signal);
            gOptimalBusyCells_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, best_busy_cell);
            
            gOptimalPhotJitter_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, best_sigma_phot);
            gOptimalDCRJitter_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, best_sigma_DCR);
            gOptimalPhotJitter_vs_Lumi_double[iSiPM]->SetPoint(iLumi, iLumi, best_sigma_phot/sqrt(2));
            gOptimalDCRJitter_vs_Lumi_double[iSiPM]->SetPoint(iLumi, iLumi, best_sigma_DCR/sqrt(2));  
            
            
            if (iLumi == TOTLUMI-1)
            {
                for (int iT = 0; iT < 70; iT++)
                {
                    float temperature = -50+iT;
                    gCurrent_vs_Temperature[iSiPM]->SetPoint(iT, temperature, opt_current/pow(T_coeff[iSiPM], (ref_temp-temperature)/10));
                    gDCR_vs_Temperature[iSiPM]->SetPoint(iT, temperature, opt_DCR/pow(T_coeff[iSiPM], (ref_temp-temperature)/10));
                    
                }
            }
            
            

            
        }//end of SiPM loop
        
        //to fill with average resolution
        /*
        gOptimalCTR_TP_vs_Lumi_Ave      ->SetPoint(iLumi, iLumi, (max_CTR_TP4x4+min_CTR_TP4x4)/2.);
        gOptimalCTR_TP3x3_vs_Lumi_Ave   ->SetPoint(iLumi, iLumi, (max_CTR_TP3x3+min_CTR_TP3x3)/2.);
        gOptimalCTR_TP2x2_vs_Lumi_Ave   ->SetPoint(iLumi, iLumi, (max_CTR_TP2x2+min_CTR_TP2x2)/2.);
        gOptimalCTR_vs_Lumi_Ave         ->SetPoint(iLumi, iLumi, (max_CTR_Single+min_CTR_Single)/2.);
        gOptimalCTR_2x2_vs_Lumi_Ave         ->SetPoint(iLumi, iLumi, (max_CTR_2x2+min_CTR_2x2)/2.);
        gOptimalCTR_Double_vs_Lumi_Ave  ->SetPoint(iLumi, iLumi, (max_CTR_Double+min_CTR_Double)/2.);
        gOptimalCTR_Double2x2_vs_Lumi_Ave  ->SetPoint(iLumi, iLumi, (max_CTR_Double2x2+min_CTR_Double2x2)/2.);
        
        gOptimalCTR_TP_vs_Lumi_Ave      ->SetPointError(iLumi, 0, (max_CTR_TP4x4-min_CTR_TP4x4)/2.);
        gOptimalCTR_TP3x3_vs_Lumi_Ave   ->SetPointError(iLumi, 0, (max_CTR_TP3x3-min_CTR_TP3x3)/2.);
        gOptimalCTR_TP2x2_vs_Lumi_Ave   ->SetPointError(iLumi, 0, (max_CTR_TP2x2-min_CTR_TP2x2)/2.);
        gOptimalCTR_vs_Lumi_Ave         ->SetPointError(iLumi, 0, (max_CTR_Single-min_CTR_Single)/2.);
        gOptimalCTR_2x2_vs_Lumi_Ave         ->SetPointError(iLumi, 0, (max_CTR_2x2-min_CTR_2x2)/2.);
        gOptimalCTR_Double_vs_Lumi_Ave  ->SetPointError(iLumi, 0, (max_CTR_Double-min_CTR_Double)/2.);
        gOptimalCTR_Double2x2_vs_Lumi_Ave  ->SetPointError(iLumi, 0, (max_CTR_Double2x2-min_CTR_Double2x2)/2.);
        */
        
        //to fill with best time resolution
        gOptimalCTR_TP_vs_Lumi_Ave      ->SetPoint(iLumi, iLumi, min_CTR_TP4x4);
        gOptimalCTR_TP3x3_vs_Lumi_Ave   ->SetPoint(iLumi, iLumi, min_CTR_TP3x3);
        gOptimalCTR_TP2x2_vs_Lumi_Ave   ->SetPoint(iLumi, iLumi, min_CTR_TP2x2);
        gOptimalCTR_vs_Lumi_Ave         ->SetPoint(iLumi, iLumi, min_CTR_Single);
        gOptimalCTR_2x2_vs_Lumi_Ave         ->SetPoint(iLumi, iLumi, min_CTR_2x2);
        gOptimalCTR_Double_vs_Lumi_Ave  ->SetPoint(iLumi, iLumi, min_CTR_Double);
        gOptimalCTR_Double2x2_vs_Lumi_Ave  ->SetPoint(iLumi, iLumi, min_CTR_Double2x2);
        
        gOptimalCTR_TP_vs_Lumi_Ave      ->SetPointError(iLumi, 0, 0, 0, max_CTR_TP4x4-min_CTR_TP4x4);
        gOptimalCTR_TP3x3_vs_Lumi_Ave   ->SetPointError(iLumi, 0, 0, 0, max_CTR_TP3x3-min_CTR_TP3x3);
        gOptimalCTR_TP2x2_vs_Lumi_Ave   ->SetPointError(iLumi, 0, 0, 0, max_CTR_TP2x2-min_CTR_TP2x2);
        gOptimalCTR_vs_Lumi_Ave         ->SetPointError(iLumi, 0, 0, 0, max_CTR_Single-min_CTR_Single);
        gOptimalCTR_2x2_vs_Lumi_Ave         ->SetPointError(iLumi, 0, 0, 0, max_CTR_2x2-min_CTR_2x2);
        gOptimalCTR_Double_vs_Lumi_Ave  ->SetPointError(iLumi, 0, 0, 0, max_CTR_Double-min_CTR_Double);
        gOptimalCTR_Double2x2_vs_Lumi_Ave  ->SetPointError(iLumi, 0, 0, 0, max_CTR_Double2x2-min_CTR_Double2x2);
        
    }
    
    

    for (int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        gCurrent_vs_Lumi[iSiPM]->SetLineWidth(2);
        gGain_vs_Lumi[iSiPM]->SetLineWidth(2);
        gPower_vs_Lumi[iSiPM]->SetLineWidth(2);
        gPowerDynamic_vs_Lumi[iSiPM]->SetLineWidth(2);
        gPowerTot_vs_Lumi[iSiPM]->SetLineWidth(2);
        gDCR_vs_Lumi[iSiPM]->SetLineWidth(2);
        gCTR_vs_Lumi[iSiPM]->SetLineWidth(2);
        gOV_vs_Lumi[iSiPM]->SetLineWidth(2);
        gPDE_vs_Lumi[iSiPM]->SetLineWidth(2);
        gSignal_vs_Lumi[iSiPM]->SetLineWidth(2);
        gEleSignal_vs_Lumi[iSiPM]->SetLineWidth(2);
        gBusyCells_vs_Lumi[iSiPM]->SetLineWidth(2);
        
        
        gCurrent_vs_Lumi[iSiPM]->SetLineStyle(7);
        gGain_vs_Lumi[iSiPM]->SetLineStyle(7);
        gPower_vs_Lumi[iSiPM]->SetLineStyle(7);
        gPowerDynamic_vs_Lumi[iSiPM]->SetLineStyle(7);
        gPowerTot_vs_Lumi[iSiPM]->SetLineStyle(7);       
        gDCR_vs_Lumi[iSiPM]->SetLineStyle(7);
        gCTR_vs_Lumi[iSiPM]->SetLineStyle(7);
        gOV_vs_Lumi[iSiPM]->SetLineStyle(7);
        gPDE_vs_Lumi[iSiPM]->SetLineStyle(7);
        gSignal_vs_Lumi[iSiPM]->SetLineStyle(7);
        gEleSignal_vs_Lumi[iSiPM]->SetLineStyle(7);
        gBusyCells_vs_Lumi[iSiPM]->SetLineStyle(7);
        
        gOptimalPDE_vs_Lumi[iSiPM]->SetLineWidth(2);
        gOptimalDCR_vs_Lumi[iSiPM]->SetLineWidth(2);
        gOptimalGain_vs_Lumi[iSiPM]->SetLineWidth(2);
        gOptimalOV_vs_Lumi[iSiPM]->SetLineWidth(2);
        gOptimalSNR_vs_Lumi[iSiPM]->SetLineWidth(2);
        gOptimalCTR_vs_Lumi[iSiPM]->SetLineWidth(2);
        gOptimalCurrent_vs_Lumi[iSiPM]->SetLineWidth(2);
        gOptimalPower_vs_Lumi[iSiPM]->SetLineWidth(2);
        gOptimalPowerDynamic_vs_Lumi[iSiPM]->SetLineWidth(2);
        gOptimalPowerTot_vs_Lumi[iSiPM]->SetLineWidth(2);
        gOptimalSignal_vs_Lumi[iSiPM]->SetLineWidth(2);
        gOptimalEleSignal_vs_Lumi[iSiPM]->SetLineWidth(2);
        gOptimalBusyCells_vs_Lumi[iSiPM]->SetLineWidth(2);
        gOptimalPhotJitter_vs_Lumi[iSiPM]->SetLineWidth(2);
        gOptimalDCRJitter_vs_Lumi[iSiPM]->SetLineWidth(2);
        gOptimalPhotJitter_vs_Lumi_double[iSiPM]->SetLineWidth(2);
        gOptimalDCRJitter_vs_Lumi_double[iSiPM]->SetLineWidth(2);
        
    }
    
    int color_constant[NSIPM];
    int color_optimized[NSIPM];
    color_constant[0]  = kGreen+3;
    color_optimized[0] = kGreen+1;
        
    color_constant[1]  = kBlue+2;    
    color_optimized[1] = kCyan+2;
    
    color_constant[2]  = kRed+3;    
    color_optimized[2] = kRed+1;
    
    //color codes
    for (int iSiPM = 0; iSiPM  < NSIPM; iSiPM++)
    {
        gCurrent_vs_Lumi[iSiPM]->SetLineColor(color_constant[iSiPM]);
        gGain_vs_Lumi[iSiPM]->SetLineColor(color_constant[iSiPM]);
        gPower_vs_Lumi[iSiPM]->SetLineColor(color_constant[iSiPM]);
        gPowerDynamic_vs_Lumi[iSiPM]->SetLineColor(color_constant[iSiPM]);
        gPowerTot_vs_Lumi[iSiPM]->SetLineColor(color_constant[iSiPM]);
        gDCR_vs_Lumi[iSiPM]->SetLineColor(color_constant[iSiPM]);
        gCTR_vs_Lumi[iSiPM]->SetLineColor(color_constant[iSiPM]);
        gOV_vs_Lumi[iSiPM]->SetLineColor(color_constant[iSiPM]);
        gPDE_vs_Lumi[iSiPM]->SetLineColor(color_constant[iSiPM]);
        gSignal_vs_Lumi[iSiPM]->SetLineColor(color_constant[iSiPM]);
        gEleSignal_vs_Lumi[iSiPM]->SetLineColor(color_constant[iSiPM]);
        gBusyCells_vs_Lumi[iSiPM]->SetLineColor(color_constant[iSiPM]);
        
        gOptimalPDE_vs_Lumi[iSiPM]->SetLineColor(color_optimized[iSiPM]);
        gOptimalDCR_vs_Lumi[iSiPM]->SetLineColor(color_optimized[iSiPM]);
        gOptimalGain_vs_Lumi[iSiPM]->SetLineColor(color_optimized[iSiPM]);
        gOptimalOV_vs_Lumi[iSiPM]->SetLineColor(color_optimized[iSiPM]);
        gOptimalCTR_vs_Lumi[iSiPM]->SetLineColor(color_optimized[iSiPM]);
        gOptimalCurrent_vs_Lumi[iSiPM]->SetLineColor(color_optimized[iSiPM]);
        gOptimalPower_vs_Lumi[iSiPM]->SetLineColor(color_optimized[iSiPM]);
        gOptimalPowerDynamic_vs_Lumi[iSiPM]->SetLineColor(color_optimized[iSiPM]);
        gOptimalPowerTot_vs_Lumi[iSiPM]->SetLineColor(color_optimized[iSiPM]);
        gOptimalSignal_vs_Lumi[iSiPM]->SetLineColor(color_optimized[iSiPM]);
        gOptimalEleSignal_vs_Lumi[iSiPM]->SetLineColor(color_optimized[iSiPM]);
        gOptimalBusyCells_vs_Lumi[iSiPM]->SetLineColor(color_optimized[iSiPM]);        
    
        gCurrent_vs_Lumi[iSiPM]->SetMarkerColor(color_constant[iSiPM]);
        gGain_vs_Lumi[iSiPM]->SetMarkerColor(color_constant[iSiPM]);
        gPower_vs_Lumi[iSiPM]->SetMarkerColor(color_constant[iSiPM]);
        gPowerDynamic_vs_Lumi[iSiPM]->SetMarkerColor(color_constant[iSiPM]);
        gPowerTot_vs_Lumi[iSiPM]->SetMarkerColor(color_constant[iSiPM]);
        gDCR_vs_Lumi[iSiPM]->SetMarkerColor(color_constant[iSiPM]);
        gCTR_vs_Lumi[iSiPM]->SetMarkerColor(color_constant[iSiPM]);
        gOV_vs_Lumi[iSiPM]->SetMarkerColor(color_constant[iSiPM]);
        gPDE_vs_Lumi[iSiPM]->SetMarkerColor(color_constant[iSiPM]);
        gSignal_vs_Lumi[iSiPM]->SetMarkerColor(color_constant[iSiPM]);
        gEleSignal_vs_Lumi[iSiPM]->SetMarkerColor(color_constant[iSiPM]);
        gBusyCells_vs_Lumi[iSiPM]->SetMarkerColor(color_constant[iSiPM]);
        
        gOptimalPDE_vs_Lumi[iSiPM]->SetMarkerColor(color_optimized[iSiPM]);
        gOptimalGain_vs_Lumi[iSiPM]->SetMarkerColor(color_optimized[iSiPM]);
        gOptimalDCR_vs_Lumi[iSiPM]->SetMarkerColor(color_optimized[iSiPM]);
        gOptimalOV_vs_Lumi[iSiPM]->SetMarkerColor(color_optimized[iSiPM]);
        gOptimalCTR_vs_Lumi[iSiPM]->SetMarkerColor(color_optimized[iSiPM]);
        gOptimalCurrent_vs_Lumi[iSiPM]->SetMarkerColor(color_optimized[iSiPM]);
        gOptimalPower_vs_Lumi[iSiPM]->SetMarkerColor(color_optimized[iSiPM]);
        gOptimalPowerTot_vs_Lumi[iSiPM]->SetMarkerColor(color_optimized[iSiPM]);
        gOptimalPowerDynamic_vs_Lumi[iSiPM]->SetMarkerColor(color_optimized[iSiPM]);
        gOptimalSignal_vs_Lumi[iSiPM]->SetMarkerColor(color_optimized[iSiPM]);
        gOptimalEleSignal_vs_Lumi[iSiPM]->SetMarkerColor(color_optimized[iSiPM]);
        gOptimalBusyCells_vs_Lumi[iSiPM]->SetMarkerColor(color_optimized[iSiPM]);
    }

    
    
    TCanvas * cCurrent_vs_Lumi = new TCanvas ("cCurrent_vs_Lumi", "cCurrent_vs_Lumi", 600, 500);
    cCurrent_vs_Lumi->SetLeftMargin(0.12);
    gCurrent_vs_Lumi[0]->Draw("ALE");
    gCurrent_vs_Lumi[0]->GetXaxis()->SetTitle("Integrated Luminosity [fb^{-1}]");
    gCurrent_vs_Lumi[0]->GetYaxis()->SetTitle("Current / SiPM [#muA]");
    gCurrent_vs_Lumi[0]->GetXaxis()->SetTitleSize(0.05);
    gCurrent_vs_Lumi[0]->GetYaxis()->SetTitleSize(0.05);
    gCurrent_vs_Lumi[0]->GetXaxis()->SetTitleOffset(0.85);
    gCurrent_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gCurrent_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 3500);
    gCurrent_vs_Lumi[1]->Draw("same LE");
    gCurrent_vs_Lumi[2]->Draw("same LE");
    gOptimalCurrent_vs_Lumi[0]->Draw("same LE");
    gOptimalCurrent_vs_Lumi[1]->Draw("same LE");
    gOptimalCurrent_vs_Lumi[2]->Draw("same LE");
    gPad->SetGridy();
    
    leg = new TLegend(0.16,0.64,0.65,0.88,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
 
    for (int iSiPM = 0; iSiPM<NSIPM; iSiPM++)
    {
        leg ->AddEntry(gCurrent_vs_Lumi[iSiPM], Form("constant bias = +1.5V OV (%s)", sipm_name[iSiPM].c_str()), "lp");
        leg ->AddEntry(gOptimalCurrent_vs_Lumi[iSiPM], Form("optimized bias (%s)", sipm_name[iSiPM].c_str()), "lp");
    }
    leg->Draw();
        
    cCurrent_vs_Lumi->cd();
    outdaq = Form("%sCurrent_vs_Lumi.png", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cCurrent_vs_Lumi->SaveAs(daqfile);
    outdaq = Form("%sCurrent_vs_Lumi.pdf", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cCurrent_vs_Lumi->SaveAs(daqfile);
    
    
    TCanvas * cPower_vs_Lumi = new TCanvas ("cPower_vs_Lumi", "cPower_vs_Lumi", 600, 500);
    cPower_vs_Lumi->SetLeftMargin(0.12);
    gPower_vs_Lumi[0]->Draw("ALE");
    gPower_vs_Lumi[0]->GetXaxis()->SetTitle("Integrated Luminosity [fb^{-1}]");
    gPower_vs_Lumi[0]->GetYaxis()->SetTitle("Static power / SiPM [mW]");
    gPower_vs_Lumi[0]->GetXaxis()->SetTitleSize(0.05);
    gPower_vs_Lumi[0]->GetYaxis()->SetTitleSize(0.05);    
    gPower_vs_Lumi[0]->GetXaxis()->SetTitleOffset(0.85);
    gPower_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gPower_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 160);
    gPower_vs_Lumi[1]->Draw("same LE");
    gOptimalPower_vs_Lumi[0]->Draw("same LE");
    gOptimalPower_vs_Lumi[1]->Draw("same LE");
    gPower_vs_Lumi[2]->Draw("same LE");
    gOptimalPower_vs_Lumi[2]->Draw("same LE");
    gPad->SetGridy();
    leg->Draw();
    
    cPower_vs_Lumi->cd();
    outdaq = Form("%sPower_vs_Lumi.png", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cPower_vs_Lumi->SaveAs(daqfile);
    outdaq = Form("%sPower_vs_Lumi.pdf", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cPower_vs_Lumi->SaveAs(daqfile);
    
    
    TCanvas * cPowerDyn_vs_Lumi = new TCanvas ("cPowerDyn_vs_Lumi", "cPowerDyn_vs_Lumi", 600, 500);
    cPowerDyn_vs_Lumi->SetLeftMargin(0.12);
    gPowerDynamic_vs_Lumi[0]->Draw("ALE");
    gPowerDynamic_vs_Lumi[0]->GetXaxis()->SetTitle("Integrated Luminosity [fb^{-1}]");
    gPowerDynamic_vs_Lumi[0]->GetYaxis()->SetTitle("Dynamic power / SiPM [mW]");
    gPowerDynamic_vs_Lumi[0]->GetXaxis()->SetTitleSize(0.05);
    gPowerDynamic_vs_Lumi[0]->GetYaxis()->SetTitleSize(0.05);    
    gPowerDynamic_vs_Lumi[0]->GetXaxis()->SetTitleOffset(0.85);
    gPowerDynamic_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gPowerDynamic_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 180);
    gPowerDynamic_vs_Lumi[1]->Draw("same LE");
    gOptimalPowerDynamic_vs_Lumi[0]->Draw("same LE");
    gOptimalPowerDynamic_vs_Lumi[1]->Draw("same LE");
    gPowerDynamic_vs_Lumi[2]->Draw("same LE");
    gOptimalPowerDynamic_vs_Lumi[2]->Draw("same LE");
    gPad->SetGridy();
    leg->Draw();
    
    cPowerDyn_vs_Lumi->cd();
    outdaq = Form("%sPowerDyn_vs_Lumi.png", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cPowerDyn_vs_Lumi->SaveAs(daqfile);
    outdaq = Form("%sPowerDyn_vs_Lumi.pdf", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cPowerDyn_vs_Lumi->SaveAs(daqfile);
    
    
    TCanvas * cPowerTot_vs_Lumi = new TCanvas ("cPowerTot_vs_Lumi", "cPowerTot_vs_Lumi", 600, 500);
    cPowerTot_vs_Lumi->SetLeftMargin(0.12);
    gPowerTot_vs_Lumi[0]->Draw("ALE");
    gPowerTot_vs_Lumi[0]->GetXaxis()->SetTitle("Integrated Luminosity [fb^{-1}]");
    gPowerTot_vs_Lumi[0]->GetYaxis()->SetTitle("Total power / SiPM [mW]");
    gPowerTot_vs_Lumi[0]->GetXaxis()->SetTitleSize(0.05);
    gPowerTot_vs_Lumi[0]->GetYaxis()->SetTitleSize(0.05);    
    gPowerTot_vs_Lumi[0]->GetXaxis()->SetTitleOffset(0.85);
    gPowerTot_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gPowerTot_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 180);
    gPowerTot_vs_Lumi[1]->Draw("same LE");
    gOptimalPowerTot_vs_Lumi[0]->Draw("same LE");
    gOptimalPowerTot_vs_Lumi[1]->Draw("same LE");
    gPowerTot_vs_Lumi[2]->Draw("same LE");
    gOptimalPowerTot_vs_Lumi[2]->Draw("same LE");
    gPad->SetGridy();
    leg->Draw();
    
    cPowerTot_vs_Lumi->cd();
    outdaq = Form("%sPowerTot_vs_Lumi.png", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cPowerTot_vs_Lumi->SaveAs(daqfile);
    outdaq = Form("%sPowerTot_vs_Lumi.pdf", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cPowerTot_vs_Lumi->SaveAs(daqfile);
    

    
    
    

    TCanvas * cOV_vs_Lumi = new TCanvas ("cOV_vs_Lumi", "cOV_vs_Lumi", 600, 500);    
    cOV_vs_Lumi->SetLeftMargin(0.12);
    gOptimalOV_vs_Lumi[0]->Draw("ALE");
    gOptimalOV_vs_Lumi[0]->GetXaxis()->SetTitle("Integrated Luminosity [fb^{-1}]");
    gOptimalOV_vs_Lumi[0]->GetYaxis()->SetTitle("Over Voltage [V]");
    gOptimalOV_vs_Lumi[0]->GetXaxis()->SetTitleSize(0.05);
    gOptimalOV_vs_Lumi[0]->GetYaxis()->SetTitleSize(0.05);
    gOptimalOV_vs_Lumi[0]->GetXaxis()->SetTitleOffset(0.85);
    gOptimalOV_vs_Lumi[0]->GetYaxis()->SetTitleOffset(1.1);
    gOptimalOV_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gOptimalOV_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 6);
    gOV_vs_Lumi[0]->Draw("same LE");    
    gOV_vs_Lumi[1]->Draw("same LE");    
    gOptimalOV_vs_Lumi[1]->Draw("same LE");
    gOV_vs_Lumi[2]->Draw("same LE");
    gOptimalOV_vs_Lumi[2]->Draw("same LE");
    gPad->SetGridy();
    leg->Draw();
    
    cOV_vs_Lumi->cd();
    outdaq = Form("%sOptimalOV_vs_Lumi.png", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cOV_vs_Lumi->SaveAs(daqfile);
    outdaq = Form("%sOptimalOV_vs_Lumi.pdf", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cOV_vs_Lumi->SaveAs(daqfile);

    TCanvas * cPDE_vs_Lumi = new TCanvas ("cPDE_vs_Lumi", "cPDE_vs_Lumi", 600, 500);
    cPDE_vs_Lumi->SetLeftMargin(0.12);
    gOptimalPDE_vs_Lumi[0]->Draw("ALE");
    gOptimalPDE_vs_Lumi[0]->GetXaxis()->SetTitle("Integrated Luminosity [fb^{-1}]");
    gOptimalPDE_vs_Lumi[0]->GetYaxis()->SetTitle("PDE");
    gOptimalPDE_vs_Lumi[0]->GetXaxis()->SetTitleSize(0.05);
    gOptimalPDE_vs_Lumi[0]->GetYaxis()->SetTitleSize(0.05);
    gOptimalPDE_vs_Lumi[0]->GetXaxis()->SetTitleOffset(0.85);
    gOptimalPDE_vs_Lumi[0]->GetYaxis()->SetTitleOffset(1.1);
    gOptimalPDE_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gOptimalPDE_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 0.6);
    gPDE_vs_Lumi[1]->Draw("same LE");
    gPDE_vs_Lumi[0]->Draw("same LE");
    gOptimalPDE_vs_Lumi[1]->Draw("same LE");
    gPDE_vs_Lumi[2]->Draw("same LE");
    gOptimalPDE_vs_Lumi[2]->Draw("same LE");
    gPad->SetGridy();
    leg->Draw();
    
    cPDE_vs_Lumi->cd();
    outdaq = Form("%sOptimalPDE_vs_Lumi.png", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cPDE_vs_Lumi->SaveAs(daqfile);
    outdaq = Form("%sOptimalPDE_vs_Lumi.pdf", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cPDE_vs_Lumi->SaveAs(daqfile);
    
    TCanvas * cSignal_vs_Lumi = new TCanvas ("cSignal_vs_Lumi", "cSignal_vs_Lumi", 600, 500);
    cSignal_vs_Lumi->SetLeftMargin(0.12);
    gOptimalSignal_vs_Lumi[0]->Draw("ALE");
    gOptimalSignal_vs_Lumi[0]->GetXaxis()->SetTitle("Integrated Luminosity [fb^{-1}]");
    gOptimalSignal_vs_Lumi[0]->GetYaxis()->SetTitle("Signal [phe]");
    gOptimalSignal_vs_Lumi[0]->GetXaxis()->SetTitleSize(0.05);
    gOptimalSignal_vs_Lumi[0]->GetYaxis()->SetTitleSize(0.05);
    gOptimalSignal_vs_Lumi[0]->GetXaxis()->SetTitleOffset(0.85);
    gOptimalSignal_vs_Lumi[0]->GetYaxis()->SetTitleOffset(1.15);
    gOptimalSignal_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gOptimalSignal_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 20000);
    gSignal_vs_Lumi[1]->Draw("same LE");
    gSignal_vs_Lumi[0]->Draw("same LE");
    gOptimalSignal_vs_Lumi[1]->Draw("same LE");
    gSignal_vs_Lumi[2]->Draw("same LE");
    gOptimalSignal_vs_Lumi[2]->Draw("same LE");
    gPad->SetGridy();
    leg->Draw();
    
    cSignal_vs_Lumi->cd();
    outdaq = Form("%sOptimalSignal_vs_Lumi.png", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cSignal_vs_Lumi->SaveAs(daqfile);
    outdaq = Form("%sOptimalSignal_vs_Lumi.pdf", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cSignal_vs_Lumi->SaveAs(daqfile);
    
    
    TCanvas * cEleSignal_vs_Lumi = new TCanvas ("cEleSignal_vs_Lumi", "cEleSignal_vs_Lumi", 600, 500);
    cEleSignal_vs_Lumi->SetLeftMargin(0.12);
    gOptimalEleSignal_vs_Lumi[0]->Draw("ALE");
    gOptimalEleSignal_vs_Lumi[0]->GetXaxis()->SetTitle("Integrated Luminosity [fb^{-1}]");
    gOptimalEleSignal_vs_Lumi[0]->GetYaxis()->SetTitle("Electric Signal [electrons]");
    gOptimalEleSignal_vs_Lumi[0]->GetXaxis()->SetTitleSize(0.05);
    gOptimalEleSignal_vs_Lumi[0]->GetYaxis()->SetTitleSize(0.05);
    gOptimalEleSignal_vs_Lumi[0]->GetXaxis()->SetTitleOffset(0.85);
    gOptimalEleSignal_vs_Lumi[0]->GetYaxis()->SetTitleOffset(1.1);
    gOptimalEleSignal_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gOptimalEleSignal_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 10000e6);
    gEleSignal_vs_Lumi[1]->Draw("same LE");
    gEleSignal_vs_Lumi[0]->Draw("same LE");
    gOptimalEleSignal_vs_Lumi[1]->Draw("same LE");
    gEleSignal_vs_Lumi[2]->Draw("same LE");
    gOptimalEleSignal_vs_Lumi[2]->Draw("same LE");
    gPad->SetGridy();
    leg->Draw();
    
    cEleSignal_vs_Lumi->cd();
    outdaq = Form("%sOptimalEleSignal_vs_Lumi.png", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cEleSignal_vs_Lumi->SaveAs(daqfile);
    outdaq = Form("%sOptimalEleSignal_vs_Lumi.pdf", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cEleSignal_vs_Lumi->SaveAs(daqfile);
    
    TCanvas * cBusyCells_vs_Lumi = new TCanvas ("cBusyCells_vs_Lumi", "cBusyCells_vs_Lumi", 600, 500);
    cBusyCells_vs_Lumi->SetLeftMargin(0.12);
    gOptimalBusyCells_vs_Lumi[0]->Draw("ALE");
    gOptimalBusyCells_vs_Lumi[0]->GetXaxis()->SetTitle("Integrated Luminosity [fb^{-1}]");
    gOptimalBusyCells_vs_Lumi[0]->GetYaxis()->SetTitle("DCR cell occupancy");
    gOptimalBusyCells_vs_Lumi[0]->GetXaxis()->SetTitleSize(0.05);
    gOptimalBusyCells_vs_Lumi[0]->GetYaxis()->SetTitleSize(0.05);
    gOptimalBusyCells_vs_Lumi[0]->GetXaxis()->SetTitleOffset(0.85);
    gOptimalBusyCells_vs_Lumi[0]->GetYaxis()->SetTitleOffset(1.1);
    gOptimalBusyCells_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gOptimalBusyCells_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 0.06);
    gBusyCells_vs_Lumi[1]->Draw("same LE");
    gBusyCells_vs_Lumi[0]->Draw("same LE");
    gOptimalBusyCells_vs_Lumi[1]->Draw("same LE");
    gBusyCells_vs_Lumi[2]->Draw("same LE");
    gOptimalBusyCells_vs_Lumi[2]->Draw("same LE");
    gPad->SetGridy();
    leg->Draw();
    
    cBusyCells_vs_Lumi->cd();
    outdaq = Form("%sBusyCells_vs_Lumi.png", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cBusyCells_vs_Lumi->SaveAs(daqfile);
    outdaq = Form("%sBusyCells_vs_Lumi.pdf", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cBusyCells_vs_Lumi->SaveAs(daqfile);

    TCanvas * cDCR_vs_Lumi = new TCanvas ("cDCR_vs_Lumi", "cDCR_vs_Lumi", 600, 500);
    cDCR_vs_Lumi->SetLeftMargin(0.12);
    gDCR_vs_Lumi[0]->Draw("ALE");
    gDCR_vs_Lumi[0]->GetXaxis()->SetTitle("Integrated Luminosity [fb^{-1}]");
    gDCR_vs_Lumi[0]->GetYaxis()->SetTitle("DCR / SiPM [GHz]");
    gDCR_vs_Lumi[0]->GetXaxis()->SetTitleSize(0.05);
    gDCR_vs_Lumi[0]->GetYaxis()->SetTitleSize(0.05);  
    gDCR_vs_Lumi[0]->GetXaxis()->SetTitleOffset(0.85);
    gDCR_vs_Lumi[0]->GetYaxis()->SetTitleOffset(1.1);
    gDCR_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gDCR_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 120);
    gDCR_vs_Lumi[1]->Draw("same LE");
    gOptimalDCR_vs_Lumi[0]->Draw("same LE");
    gOptimalDCR_vs_Lumi[1]->Draw("same LE");
    gDCR_vs_Lumi[2]->Draw("same LE");
    gOptimalDCR_vs_Lumi[2]->Draw("same LE");
    gPad->SetGridy();
    leg->Draw();
    
    cDCR_vs_Lumi->cd();
    outdaq = Form("%sDCR_vs_Lumi.png", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cDCR_vs_Lumi->SaveAs(daqfile);
    outdaq = Form("%sDCR_vs_Lumi.pdf", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cDCR_vs_Lumi->SaveAs(daqfile);
    
    
    TCanvas * cGain_vs_Lumi = new TCanvas ("cGain_vs_Lumi", "cGain_vs_Lumi", 600, 500);
    cGain_vs_Lumi->SetLeftMargin(0.12);
    gGain_vs_Lumi[0]->Draw("ALE");
    gGain_vs_Lumi[0]->GetXaxis()->SetTitle("Integrated Luminosity [fb^{-1}]");
    gGain_vs_Lumi[0]->GetYaxis()->SetTitle("Gain");
    gGain_vs_Lumi[0]->GetXaxis()->SetTitleSize(0.05);
    gGain_vs_Lumi[0]->GetYaxis()->SetTitleSize(0.05);  
    gGain_vs_Lumi[0]->GetXaxis()->SetTitleOffset(0.85);
    gGain_vs_Lumi[0]->GetYaxis()->SetTitleOffset(1.1);
    gGain_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gGain_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 8e5);
    gGain_vs_Lumi[1]->Draw("same LE");
    gOptimalGain_vs_Lumi[0]->Draw("same LE");
    gOptimalGain_vs_Lumi[1]->Draw("same LE");
    gGain_vs_Lumi[2]->Draw("same LE");
    gOptimalGain_vs_Lumi[2]->Draw("same LE");
    gPad->SetGridy();
    leg->Draw();
    
    cGain_vs_Lumi->cd();
    outdaq = Form("%sGain_vs_Lumi.png", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cGain_vs_Lumi->SaveAs(daqfile);
    outdaq = Form("%sGain_vs_Lumi.pdf", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cGain_vs_Lumi->SaveAs(daqfile);
    
    
    
    TCanvas * cPowerCompare[NSIPM];
    
    for (int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        cPowerCompare[iSiPM] = new TCanvas (Form("cPowerCompare_%d", iSiPM), Form("cPowerCompare_%d", iSiPM), 600, 500);
        cPowerCompare[iSiPM]->cd();
        gOptimalPowerTot_vs_Lumi[iSiPM]->Draw("ALE");
        gOptimalPowerTot_vs_Lumi[iSiPM]->SetTitle(sipm_name[iSiPM].c_str());
        gOptimalPowerTot_vs_Lumi[iSiPM]->GetXaxis()->SetTitle("Integrated Luminosity [fb^{-1}]");
        gOptimalPowerTot_vs_Lumi[iSiPM]->GetYaxis()->SetTitle("Power / SiPM [mW]");
        gOptimalPowerTot_vs_Lumi[iSiPM]->GetXaxis()->SetTitleSize(0.05);
        gOptimalPowerTot_vs_Lumi[iSiPM]->GetYaxis()->SetTitleSize(0.05);    
        gOptimalPowerTot_vs_Lumi[iSiPM]->GetXaxis()->SetTitleOffset(0.85);
        gOptimalPowerTot_vs_Lumi[iSiPM]->GetXaxis()->SetRangeUser(0, 4000);
        gOptimalPowerTot_vs_Lumi[iSiPM]->GetYaxis()->SetRangeUser(0, 180);
        
        gOptimalPowerTot_vs_Lumi[iSiPM]->SetLineColor(kBlack);
        gOptimalPowerDynamic_vs_Lumi[iSiPM]->SetLineColor(kYellow+2);
        gOptimalPower_vs_Lumi[iSiPM]->SetLineColor(kCyan+1);
        
        gOptimalPowerDynamic_vs_Lumi[iSiPM]->Draw("same LE");
        gOptimalPower_vs_Lumi[iSiPM]->Draw("same LE");                        
        
        TLegend *leg_power = new TLegend(0.15,0.7,0.55,0.88,NULL,"brNDC");
        leg_power->SetBorderSize(0);
        leg_power->SetTextFont(42);
        leg_power->SetTextSize(0.03);
        leg_power->SetLineColor(1);
        leg_power->SetLineStyle(1);
        leg_power->SetLineWidth(1);
        leg_power->SetFillColor(0);
        
        leg_power->AddEntry(gOptimalPowerTot_vs_Lumi[iSiPM], "Total power consumtpion", "l");
        leg_power->AddEntry(gOptimalPower_vs_Lumi[iSiPM], "Static power consumtpion", "l");
        leg_power->AddEntry(gOptimalPowerDynamic_vs_Lumi[iSiPM], "Dynamic power consumtpion", "l");
        leg_power->Draw();
        
        outdaq = Form("%sPowerCompare_vs_Lumi_SiPM_bound50mW%d.png", output_folder.c_str(), iSiPM);
        daqfile = outdaq.c_str();    
        cPowerCompare[iSiPM]->SaveAs(daqfile);
        outdaq = Form("%sPowerCompare_vs_Lumi_SiPM_bound50mW%d.pdf", output_folder.c_str(), iSiPM);
        daqfile = outdaq.c_str();    
        cPowerCompare[iSiPM]->SaveAs(daqfile);
    }
    
    
    TCanvas * cCTR_vs_Lumi = new TCanvas ("cCTR_vs_Lumi", "cCTR_vs_Lumi", 600, 500);    
    cCTR_vs_Lumi->SetLeftMargin(0.12);
    gOptimalCTR_vs_Lumi[0]->GetXaxis()->SetTitle("Integrated Luminosity [fb^{-1}]");
    gOptimalCTR_vs_Lumi[0]->GetYaxis()->SetTitle("Time resolution [ps]");
    gOptimalCTR_vs_Lumi[0]->GetXaxis()->SetTitleSize(0.05);
    gOptimalCTR_vs_Lumi[0]->GetYaxis()->SetTitleSize(0.05);
    gOptimalCTR_vs_Lumi[0]->GetXaxis()->SetTitleOffset(0.85);
    gOptimalCTR_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 140);
    gOptimalCTR_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gOptimalCTR_vs_Lumi[0]->SetMaximum(140);
    gOptimalCTR_vs_Lumi[0]->SetLineColor(kGreen+1);
    gOptimalCTR_vs_Lumi[0]->SetMarkerColor(kGreen+1);
    gOptimalCTR_vs_Lumi[0]->Draw("ALE");
    
//     gCTR_vs_Lumi[0]->Draw("same LPE");
//     gCTR_vs_Lumi[1]->Draw("same LPE");
    gOptimalCTR_vs_Lumi[1]->Draw("same LE");
    gOptimalCTR_vs_Lumi[2]->Draw("same LE");
    gOptimalCTR_Double_vs_Lumi[2]->Draw("same LE"); 
    
    gOptimalCTR_Double_vs_Lumi[0]->SetLineStyle(7);
    gOptimalCTR_Double_vs_Lumi[0]->SetLineColor(kYellow+2);
    gOptimalCTR_Double_vs_Lumi[0]->SetMarkerColor(kYellow+2);
    gOptimalCTR_Double_vs_Lumi[0]->Draw("same LE");
    
    
    gOptimalCTR_Double_vs_Lumi[1]->SetLineStyle(5);
    gOptimalCTR_Double_vs_Lumi[1]->SetLineColor(kOrange+2);
    gOptimalCTR_Double_vs_Lumi[1]->SetMarkerColor(kOrange+2);
    gOptimalCTR_Double_vs_Lumi[1]->Draw("same LE");
    
    gOptimalCTR_TP_vs_Lumi[0]->SetLineStyle(7);
    gOptimalCTR_TP_vs_Lumi[0]->SetLineColor(kViolet+2);
    gOptimalCTR_TP_vs_Lumi[0]->SetMarkerColor(kViolet+2);
    gOptimalCTR_TP_vs_Lumi[0]->Draw("same LE");
    
    
    gOptimalCTR_TP_vs_Lumi[1]->SetLineStyle(5);
    gOptimalCTR_TP_vs_Lumi[1]->SetLineColor(kViolet+2);
    gOptimalCTR_TP_vs_Lumi[1]->SetMarkerColor(kViolet+2);
    gOptimalCTR_TP_vs_Lumi[1]->Draw("same LE");
    
    gOptimalCTR_TP_vs_Lumi[2]->SetLineStyle(5);
    gOptimalCTR_TP_vs_Lumi[2]->SetLineColor(kViolet+2);
    gOptimalCTR_TP_vs_Lumi[2]->SetMarkerColor(kViolet+2);
    gOptimalCTR_TP_vs_Lumi[2]->Draw("same LE");
    
    gPad->SetGridy();
    
    TLegend * leg2;
    
    leg2 = new TLegend(0.15,0.7,0.55,0.88,NULL,"brNDC");
    leg2->SetBorderSize(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.03);
    leg2->SetLineColor(1);
    leg2->SetLineStyle(1);
    leg2->SetLineWidth(1);
    leg2->SetFillColor(0);
 
    for (int iSiPM = 0; iSiPM< NSIPM; iSiPM++)
    {
//         leg2->AddEntry(gCTR_vs_Lumi[iSiPM], Form("constant bias = +1.5V OV (%s)", sipm_name[iSiPM].c_str()), "lp");    
//         leg2->AddEntry(gOptimalCTR_vs_Lumi[iSiPM], Form("optimized bias (%s)", sipm_name[iSiPM].c_str()), "lp");
//         leg2->AddEntry(gOptimalCTR_Double_vs_Lumi[iSiPM], Form("optimized bias / sqrt(2) [double read-out] (%s)", sipm_name[iSiPM].c_str()), "lp");
        
        
        leg2->AddEntry(gOptimalCTR_Double_vs_Lumi[iSiPM], Form("SiPM type: (%s)", sipm_name[iSiPM].c_str()), "lp");
        leg2->AddEntry(gOptimalCTR_vs_Lumi[iSiPM], Form("SiPM type: (%s) [single read-out]", sipm_name[iSiPM].c_str()), "lp");
        //     leg2->AddEntry(gOptimalCTR_Double_vs_Lumi[1], "optimized bias / sqrt(2) [double read-out]", "lp");
    }
    leg2->Draw();
    
    cCTR_vs_Lumi->cd();
    outdaq = Form("%sTimeRes_vs_Lumi.png", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cCTR_vs_Lumi->SaveAs(daqfile);
    outdaq = Form("%sTimeRes_vs_Lumi.pdf", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cCTR_vs_Lumi->SaveAs(daqfile);
    
    
    int myFillStyle = 3004;
    TCanvas * cCTR_vs_Lumi_Ave = new TCanvas ("cCTR_vs_Lumi_Ave", "cCTR_vs_Lumi_Ave", 600, 500);
    cCTR_vs_Lumi_Ave->SetLeftMargin(0.12);
    
    gOptimalCTR_TP_vs_Lumi_Ave->GetXaxis()->SetTitle("Integrated Luminosity [fb^{-1}]");
    gOptimalCTR_TP_vs_Lumi_Ave->GetYaxis()->SetTitle("Time resolution [ps]");
    gOptimalCTR_TP_vs_Lumi_Ave->GetXaxis()->SetTitleSize(0.05);
    gOptimalCTR_TP_vs_Lumi_Ave->GetYaxis()->SetTitleSize(0.05);
    gOptimalCTR_TP_vs_Lumi_Ave->GetYaxis()->SetRangeUser(0, 180);
    gOptimalCTR_TP_vs_Lumi_Ave->GetXaxis()->SetRangeUser(0, 4000);
    gOptimalCTR_TP_vs_Lumi_Ave->SetMinimum(0);
    gOptimalCTR_TP_vs_Lumi_Ave->SetMaximum(140);
    gOptimalCTR_TP_vs_Lumi_Ave->SetFillColor(kOrange+1);
    gOptimalCTR_TP_vs_Lumi_Ave->SetFillStyle(myFillStyle);
    gOptimalCTR_TP_vs_Lumi_Ave->SetLineColor(kOrange+1);
    gOptimalCTR_TP_vs_Lumi_Ave->SetMarkerColor(kOrange+1);
    gOptimalCTR_TP_vs_Lumi_Ave->SetLineWidth(2);
    gOptimalCTR_TP_vs_Lumi_Ave->Draw("3ALE");
    
    gOptimalCTR_TP3x3_vs_Lumi_Ave->SetFillColor(kYellow+1);
    gOptimalCTR_TP3x3_vs_Lumi_Ave->SetFillStyle(myFillStyle);
    gOptimalCTR_TP3x3_vs_Lumi_Ave->SetLineColor(kYellow+1);
    gOptimalCTR_TP3x3_vs_Lumi_Ave->SetMarkerColor(kYellow+1);
    gOptimalCTR_TP3x3_vs_Lumi_Ave->SetLineWidth(2);
    gOptimalCTR_TP3x3_vs_Lumi_Ave->Draw("same 3LE");
    
    gOptimalCTR_TP2x2_vs_Lumi_Ave->SetFillColor(kRed+1);
    gOptimalCTR_TP2x2_vs_Lumi_Ave->SetFillStyle(myFillStyle);
    gOptimalCTR_TP2x2_vs_Lumi_Ave->SetLineColor(kRed+1);
    gOptimalCTR_TP2x2_vs_Lumi_Ave->SetMarkerColor(kRed+1);
    gOptimalCTR_TP2x2_vs_Lumi_Ave->SetLineWidth(2);
    gOptimalCTR_TP2x2_vs_Lumi_Ave->Draw("same 3LE");
    
    gOptimalCTR_vs_Lumi_Ave->SetLineWidth(2);
    gOptimalCTR_vs_Lumi_Ave->SetFillColor(kBlue+1);
    gOptimalCTR_vs_Lumi_Ave->SetLineColor(kBlue+1);
    gOptimalCTR_vs_Lumi_Ave->SetMarkerColor(kBlue+1);
    gOptimalCTR_vs_Lumi_Ave->SetFillStyle(myFillStyle);
    gOptimalCTR_vs_Lumi_Ave->Draw("same 3LE");
    
    gOptimalCTR_2x2_vs_Lumi_Ave->SetLineWidth(2);
    gOptimalCTR_2x2_vs_Lumi_Ave->SetFillColor(kViolet+1);
    gOptimalCTR_2x2_vs_Lumi_Ave->SetLineColor(kViolet+1);
    gOptimalCTR_2x2_vs_Lumi_Ave->SetMarkerColor(kViolet+1);
    gOptimalCTR_2x2_vs_Lumi_Ave->SetFillStyle(myFillStyle);
    gOptimalCTR_2x2_vs_Lumi_Ave->Draw("same 3LE");
    
    gOptimalCTR_Double_vs_Lumi_Ave->SetLineWidth(2);
    gOptimalCTR_Double_vs_Lumi_Ave->SetLineColor(kGreen+1);
    gOptimalCTR_Double_vs_Lumi_Ave->SetFillColor(kGreen+1);
    gOptimalCTR_Double_vs_Lumi_Ave->SetMarkerColor(kGreen+1);
    gOptimalCTR_Double_vs_Lumi_Ave->SetFillStyle(myFillStyle);
    gOptimalCTR_Double_vs_Lumi_Ave->Draw("same 3LE");
    
    gOptimalCTR_Double2x2_vs_Lumi_Ave->SetLineWidth(2);
    gOptimalCTR_Double2x2_vs_Lumi_Ave->SetLineColor(kCyan+1);
    gOptimalCTR_Double2x2_vs_Lumi_Ave->SetFillColor(kCyan+1);
    gOptimalCTR_Double2x2_vs_Lumi_Ave->SetMarkerColor(kCyan+1);
    gOptimalCTR_Double2x2_vs_Lumi_Ave->SetFillStyle(myFillStyle);
    gOptimalCTR_Double2x2_vs_Lumi_Ave->Draw("same 3LE");
    
    leg2 = new TLegend(0.15,0.7,0.55,0.88,NULL,"brNDC");
    leg2->SetBorderSize(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.03);
    leg2->SetLineColor(1);
    leg2->SetLineStyle(1);
    leg2->SetLineWidth(1);
    leg2->SetFillColor(0);
    
    leg2->AddEntry(gOptimalCTR_TP_vs_Lumi_Ave, "TP reference (4x4 mm^{2} SiPM)", "flpe");
    leg2->AddEntry(gOptimalCTR_TP3x3_vs_Lumi_Ave, "TP reference (3x3 mm^{2} SiPM)", "flpe");
    leg2->AddEntry(gOptimalCTR_TP2x2_vs_Lumi_Ave, "TP reference (2x2 mm^{2} SiPM)", "flpe");
    leg2->AddEntry(gOptimalCTR_vs_Lumi_Ave, "TDR Option T (3x3 mm^{2} SiPM)", "flpe");
    leg2->AddEntry(gOptimalCTR_2x2_vs_Lumi_Ave, "TDR Option T (2x2 mm^{2} SiPM)", "flpe");
    leg2->AddEntry(gOptimalCTR_Double_vs_Lumi_Ave, "TDR Option B (3x3 mm^{2} SiPM)", "flpe");
    leg2->AddEntry(gOptimalCTR_Double2x2_vs_Lumi_Ave, "TDR Option B (2x2 mm^{2} SiPM)", "flpe");
    leg2->Draw();
    
//     outdaq = Form("%sTimeResAve_vs_Lumi.png", output_folder.c_str());
//     daqfile = outdaq.c_str();
//     cCTR_vs_Lumi_Ave->cd();
//     cCTR_vs_Lumi_Ave->SaveAs(daqfile);
    
    //PLOTS for CDR/TDR
        
    TCanvas * cCTR_vs_Lumi_brokedown[NSIPM];
    for (int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        cCTR_vs_Lumi_brokedown[iSiPM] = new TCanvas (Form("cCTR_vs_Lumi_brokedown_%s", sipm_name[iSiPM].c_str()), Form("cCTR_vs_Lumi_brokedown_%s", sipm_name[iSiPM].c_str()), 600, 500);    
        cCTR_vs_Lumi_brokedown[iSiPM]->SetLeftMargin(0.12);
        gOptimalCTR_vs_Lumi[iSiPM]->Draw("ALPE"); 
        gOptimalCTR_vs_Lumi[iSiPM]->SetTitle(sipm_name[iSiPM].c_str());
        gOptimalCTR_vs_Lumi[iSiPM]->GetXaxis()->SetTitle("Integrated Luminosity [fb^{-1}]");
        gOptimalCTR_vs_Lumi[iSiPM]->GetYaxis()->SetTitle("Time resolution [ps]");
        gOptimalCTR_vs_Lumi[iSiPM]->GetXaxis()->SetTitleSize(0.05);
        gOptimalCTR_vs_Lumi[iSiPM]->GetYaxis()->SetTitleSize(0.05);
        gOptimalCTR_vs_Lumi[iSiPM]->GetXaxis()->SetTitleOffset(0.85);
        gOptimalCTR_vs_Lumi[iSiPM]->GetYaxis()->SetRangeUser(0, 120);
        gOptimalCTR_vs_Lumi[iSiPM]->SetLineColor(kGreen+1);
        gOptimalCTR_vs_Lumi[iSiPM]->SetMarkerColor(kGreen+1);
        gOptimalCTR_vs_Lumi[iSiPM]->GetXaxis()->SetRangeUser(0, 4000);    
    
        gOptimalCTR_Double_vs_Lumi[iSiPM]->SetLineWidth(2);
        gOptimalCTR_Double_vs_Lumi[iSiPM]->SetLineStyle(7);
        gOptimalCTR_Double_vs_Lumi[iSiPM]->SetLineColor(kOrange+2);
        gOptimalCTR_Double_vs_Lumi[iSiPM]->SetMarkerColor(kOrange+2);
//         gOptimalCTR_Double_vs_Lumi[iSiPM]->Draw("same LE");
        
        gCTR_vs_Lumi[iSiPM]->SetLineWidth(2);        
        gCTR_vs_Lumi[iSiPM]->SetLineColor(kBlack);
        gCTR_vs_Lumi[iSiPM]->SetMarkerColor(kBlack);
        gCTR_vs_Lumi[iSiPM]->Draw("same LE");
        
        
        gOptimalPhotJitter_vs_Lumi[iSiPM]->SetLineWidth(2);        
        gOptimalPhotJitter_vs_Lumi[iSiPM]->SetLineColor(kYellow+1);
        gOptimalPhotJitter_vs_Lumi[iSiPM]->SetMarkerColor(kYellow+1);
        gOptimalPhotJitter_vs_Lumi[iSiPM]->Draw("same LE");
        
        
        gOptimalDCRJitter_vs_Lumi[iSiPM]->SetLineWidth(2);    
        gOptimalDCRJitter_vs_Lumi[iSiPM]->SetLineColor(kRed+1);
        gOptimalDCRJitter_vs_Lumi[iSiPM]->SetMarkerColor(kRed+1);
        gOptimalDCRJitter_vs_Lumi[iSiPM]->Draw("same LE");
        
        
        gDigitization->SetLineWidth(2);
        gDigitization->SetLineColor(kBlue+1);
        gDigitization->SetMarkerColor(kBlue+1);
        gDigitization->Draw("same LE");
        
        gClock->SetLineWidth(2);
        gClock->SetLineColor(kViolet+2);
        gClock->SetMarkerColor(kViolet+2);
        gClock->Draw("same LE");
        
        gElectronics->SetLineWidth(2);
        gElectronics->SetLineColor(kAzure+1);
        gElectronics->SetMarkerColor(kAzure+1);
        gElectronics->Draw("same LE");
        
        
        
        leg2 = new TLegend(0.15,0.6,0.55,0.88,NULL,"brNDC");
        leg2->SetBorderSize(0);
        leg2->SetTextFont(42);
        leg2->SetTextSize(0.03);
        leg2->SetLineColor(1);
        leg2->SetLineStyle(1);
        leg2->SetLineWidth(1);
        leg2->SetFillColor(0);
 
//         leg2->AddEntry(gOptimalCTR_vs_Lumi[iSiPM], Form("Total time resolution (%s)", sipm_name[iSiPM].c_str()), "lp");        
        leg2->AddEntry(gCTR_vs_Lumi[iSiPM], "Total time resolution (constant bias - no power limit)", "lp");        
        leg2->AddEntry(gOptimalCTR_vs_Lumi[iSiPM], "Total time resolution", "lp");        
//         leg2->AddEntry(gOptimalCTR_Double_vs_Lumi[iSiPM], "Total / sqrt(2) [double read-out]", "lp");
        leg2->AddEntry(gOptimalPhotJitter_vs_Lumi[iSiPM], "Photostatistics", "lp");
        leg2->AddEntry(gOptimalDCRJitter_vs_Lumi[iSiPM], "DCR noise", "lp");
        leg2->AddEntry(gElectronics, "Electronics", "lp");
        leg2->AddEntry(gDigitization, "Digitization", "lp");
        leg2->AddEntry(gClock, "Clock", "lp");
        leg2->Draw();
        
        
        gPad->SetGridy();
        cCTR_vs_Lumi_brokedown[iSiPM]->cd();
        outdaq = Form("%sTimeResSingleBrokeDown_vs_Lumi_%s.png", output_folder.c_str(), sipm_name[iSiPM].c_str());
        daqfile = outdaq.c_str();        
        cCTR_vs_Lumi_brokedown[iSiPM]->SaveAs(daqfile);
        outdaq = Form("%sTimeResSingleBrokeDown_vs_Lumi_%s.pdf", output_folder.c_str(), sipm_name[iSiPM].c_str());
        daqfile = outdaq.c_str();        
        cCTR_vs_Lumi_brokedown[iSiPM]->SaveAs(daqfile);
        
        
    }
    
    
    
    TCanvas * cCTR_vs_Lumi_brokedown_sqrt2[NSIPM];
    for (int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        cCTR_vs_Lumi_brokedown_sqrt2[iSiPM] = new TCanvas (Form("cCTR_vs_Lumi_brokedown_sqrt2_%s", sipm_name[iSiPM].c_str()), Form("cCTR_vs_Lumi_brokedown_sqrt2_%s", sipm_name[iSiPM].c_str()), 600, 500);    
        cCTR_vs_Lumi_brokedown_sqrt2[iSiPM]->SetLeftMargin(0.12);
//         gOptimalCTR_vs_Lumi[iSiPM]->Draw("ALE"); 
        gOptimalCTR_vs_Lumi[iSiPM]->SetTitle(sipm_name[iSiPM].c_str());
        gOptimalCTR_vs_Lumi[iSiPM]->GetXaxis()->SetTitle("Integrated Luminosity [fb^{-1}]");
        gOptimalCTR_vs_Lumi[iSiPM]->GetYaxis()->SetTitle("Time resolution [ps]");
        gOptimalCTR_vs_Lumi[iSiPM]->GetXaxis()->SetTitleSize(0.05);
        gOptimalCTR_vs_Lumi[iSiPM]->GetYaxis()->SetTitleSize(0.05);
//         gOptimalCTR_vs_Lumi[iSiPM]->GetXaxis()->SetTitleOffset(0.85);
        gOptimalCTR_vs_Lumi[iSiPM]->GetYaxis()->SetRangeUser(0, 120);
        
//         gOptimalCTR_vs_Lumi[iSiPM]->GetXaxis()->SetRangeUser(0, 4000);    
//         gOptimalCTR_vs_Lumi[iSiPM]->SetLineStyle(7);
//         gOptimalCTR_vs_Lumi[iSiPM]->SetLineWidth(1);
//         gOptimalCTR_vs_Lumi[iSiPM]->SetLineColor(kOrange+2);
//         gOptimalCTR_vs_Lumi[iSiPM]->SetMarkerColor(kOrange+2);
    
        gOptimalCTR_Double_vs_Lumi[iSiPM]->SetLineWidth(2);        
        gOptimalCTR_Double_vs_Lumi[iSiPM]->SetLineStyle(1);
        gOptimalCTR_Double_vs_Lumi[iSiPM]->SetLineColor(kGreen+1);
        gOptimalCTR_Double_vs_Lumi[iSiPM]->SetMarkerColor(kGreen+1);
//         gOptimalCTR_Double_vs_Lumi[iSiPM]->Draw("same LE");
        gOptimalCTR_Double_vs_Lumi[iSiPM]->Draw("ALE");
        gOptimalCTR_Double_vs_Lumi[iSiPM]->SetTitle(sipm_name[iSiPM].c_str());
        gOptimalCTR_Double_vs_Lumi[iSiPM]->GetXaxis()->SetTitle("Integrated Luminosity [fb^{-1}]");
        gOptimalCTR_Double_vs_Lumi[iSiPM]->GetYaxis()->SetTitle("Time resolution [ps]");
        gOptimalCTR_Double_vs_Lumi[iSiPM]->GetXaxis()->SetTitleSize(0.05);
        gOptimalCTR_Double_vs_Lumi[iSiPM]->GetXaxis()->SetTitleOffset(0.85);
        gOptimalCTR_Double_vs_Lumi[iSiPM]->GetYaxis()->SetTitleSize(0.05);
        gOptimalCTR_Double_vs_Lumi[iSiPM]->GetYaxis()->SetRangeUser(0, 120);        
        gOptimalCTR_Double_vs_Lumi[iSiPM]->GetXaxis()->SetRangeUser(0, 4000);    
        
        
        gOptimalPhotJitter_vs_Lumi_double[iSiPM]->SetLineWidth(2);        
        gOptimalPhotJitter_vs_Lumi_double[iSiPM]->SetLineColor(kYellow+1);
        gOptimalPhotJitter_vs_Lumi_double[iSiPM]->SetMarkerColor(kYellow+1);
        gOptimalPhotJitter_vs_Lumi_double[iSiPM]->Draw("same LE");
        
        
        gOptimalDCRJitter_vs_Lumi_double[iSiPM]->SetLineWidth(2);    
        gOptimalDCRJitter_vs_Lumi_double[iSiPM]->SetLineColor(kRed+1);
        gOptimalDCRJitter_vs_Lumi_double[iSiPM]->SetMarkerColor(kRed+1);
        gOptimalDCRJitter_vs_Lumi_double[iSiPM]->Draw("same LE");
        
        
        gDigitization_double->SetLineWidth(2);
        gDigitization_double->SetLineColor(kBlue+1);
        gDigitization_double->SetMarkerColor(kBlue+1);
        gDigitization_double->Draw("same LE");
        
        gClock_double->SetLineWidth(2);
        gClock_double->SetLineColor(kViolet+2);
        gClock_double->SetMarkerColor(kViolet+2);
        gClock_double->Draw("same LE");
        
        gElectronics_double->SetLineWidth(2);
        gElectronics_double->SetLineColor(kAzure+1);
        gElectronics_double->SetMarkerColor(kAzure+1);
        gElectronics_double->Draw("same LE");
        
        
        leg2 = new TLegend(0.15,0.6,0.55,0.88,NULL,"brNDC");
        leg2->SetBorderSize(0);
        leg2->SetTextFont(42);
        leg2->SetTextSize(0.03);
        leg2->SetLineColor(1);
        leg2->SetLineStyle(1);
        leg2->SetLineWidth(1);
        leg2->SetFillColor(0);
 
//         leg2->AddEntry(gOptimalCTR_vs_Lumi[iSiPM], Form("Total time resolution (%s)", sipm_name[iSiPM].c_str()), "lp");        
            
        leg2->AddEntry(gOptimalCTR_Double_vs_Lumi[iSiPM], "Total time resolution", "lp");
//         leg2->AddEntry(gOptimalCTR_vs_Lumi[iSiPM], "Total * sqrt(2) [single read-out]", "lp");    
        leg2->AddEntry(gOptimalPhotJitter_vs_Lumi_double[iSiPM], "Photostatistics", "lp");
        leg2->AddEntry(gOptimalDCRJitter_vs_Lumi_double[iSiPM], "DCR noise", "lp");
        leg2->AddEntry(gElectronics_double, "Electronics", "lp");
        leg2->AddEntry(gDigitization_double, "Digitization", "lp");
        leg2->AddEntry(gClock_double, "Clock", "lp");
        leg2->Draw();
        
        TLatex latex;
        latex.SetTextSize(0.05);
        latex.SetTextAlign(13);  //align at top
        latex.SetTextFont(42);  //helvetica, 82 for courier, 62 for helvetica bold
        latex.DrawLatex(3000,113,Form("T = %.0f#circC", my_temp));
        
        
        gPad->SetGridy();
        cCTR_vs_Lumi_brokedown_sqrt2[iSiPM]->cd();
        outdaq = Form("%sTimeResDoubleBrokeDown_vs_Lumi_%s.png", output_folder.c_str(), sipm_name[iSiPM].c_str());
        daqfile = outdaq.c_str();        
        cCTR_vs_Lumi_brokedown_sqrt2[iSiPM]->SaveAs(daqfile);
        outdaq = Form("%sTimeResDoubleBrokeDown_vs_Lumi_%s.pdf", output_folder.c_str(), sipm_name[iSiPM].c_str());
        daqfile = outdaq.c_str();        
        cCTR_vs_Lumi_brokedown_sqrt2[iSiPM]->SaveAs(daqfile);
        
        
           
        
    }
        
    leg2 = new TLegend(0.15,0.7,0.55,0.88,NULL,"brNDC");
    leg2->SetBorderSize(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.03);
    leg2->SetLineColor(1);
    leg2->SetLineStyle(1);
    leg2->SetLineWidth(1);
    leg2->SetFillColor(0);
        
    TCanvas * cTemperatureCoeff = new TCanvas ("cTemperatureCoeff", "cTemperatureCoeff", 600, 500);
    cTemperatureCoeff->cd();
    gCurrent_vs_Temperature[0]->Draw("ALPE");
    gCurrent_vs_Temperature[0]->GetXaxis()->SetTitle("Temperature [#circC]");
    gCurrent_vs_Temperature[0]->GetYaxis()->SetTitle("Current [#muA]");
    gCurrent_vs_Temperature[0]->GetYaxis()->SetRangeUser(1e2, 1e5);
    for (int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        gCurrent_vs_Temperature[iSiPM]->Draw("same LPE");
        gCurrent_vs_Temperature[iSiPM]->SetLineColor(iSiPM+1);
        gCurrent_vs_Temperature[iSiPM]->SetMarkerColor(iSiPM+1);
        gCurrent_vs_Temperature[iSiPM]->SetLineWidth(2);
        
        leg2->AddEntry(gCurrent_vs_Temperature[iSiPM], Form("%s, T_{coeff} = %.2f", sipm_name[iSiPM].c_str(), T_coeff[iSiPM]), "lp");
    }
    gCurrent_vs_Temperature[1]->SetLineColor(kGreen+2);
    gCurrent_vs_Temperature[1]->SetMarkerColor(kGreen+2);
    gCurrent_vs_Temperature[2]->SetLineColor(kRed+1);
    gCurrent_vs_Temperature[2]->SetMarkerColor(kRed+1);
        
    leg2->Draw();
    gPad->SetLogy();
    
    
    
    TCanvas * cDCR_TemperatureCoeff = new TCanvas ("cDCR_TemperatureCoeff", "cDCR_TemperatureCoeff", 600, 500);
    cDCR_TemperatureCoeff->cd();
    gDCR_vs_Temperature[0]->Draw("ALPE");
    gDCR_vs_Temperature[0]->GetXaxis()->SetTitle("Temperature [#circC]");
    gDCR_vs_Temperature[0]->GetYaxis()->SetTitle("DCR [GHz]");
    gDCR_vs_Temperature[0]->GetYaxis()->SetRangeUser(9, 1000);
    gDCR_vs_Temperature[0]->SetMinimum(9);
    for (int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        gDCR_vs_Temperature[iSiPM]->Draw("same LPE");
        gDCR_vs_Temperature[iSiPM]->SetLineColor(iSiPM+1);
        gDCR_vs_Temperature[iSiPM]->SetMarkerColor(iSiPM+1);
        gDCR_vs_Temperature[iSiPM]->SetLineWidth(2);        
    }
    gDCR_vs_Temperature[1]->SetLineColor(kGreen+2);
    gDCR_vs_Temperature[1]->SetMarkerColor(kGreen+2);
    gDCR_vs_Temperature[2]->SetLineColor(kRed+1);
    gDCR_vs_Temperature[2]->SetMarkerColor(kRed+1);
        
    leg2->Draw();
    gPad->SetLogy();
    
    
    
    
    
    
    
    fileOutput->cd();
    for (int iSiPM = 0; iSiPM<NSIPM; iSiPM++)
    {
        funcCurrent[iSiPM]->Write();
//         funcDCR[iSiPM]->Write();
        funcENC[iSiPM]->Write();
        funcGain[iSiPM]->Write();
        funcPDE[iSiPM]->Write();
        funcPDE_vsWL[iSiPM]->Write();
        
        gDCR_vs_OV[iSiPM]->SetName(Form("gDCR_vs_OV_%d", iSiPM));
//         gGain_vs_OV[iSiPM]->SetName(Form("gGain_vs_OV_%d", iSiPM));
//         gPDE_vs_OV[iSiPM]->SetName(Form("gPDE_vs_OV_%d", iSiPM));
//         gCurrent_vs_OV[iSiPM]->SetName(Form("gCurrent_vs_OV_%d", iSiPM));
    
    
        gPDE_vs_Lumi[iSiPM]     ->SetName(Form("gPDE_vs_Lumi_%d", iSiPM));
        gCurrent_vs_Lumi[iSiPM] ->SetName(Form("gCurrent_vs_Lumi_%d", iSiPM));
        gGain_vs_Lumi[iSiPM]    ->SetName(Form("gGain_vs_Lumi_%d", iSiPM));
        gDCR_vs_Lumi[iSiPM]     ->SetName(Form("gDCR_vs_Lumi_%d", iSiPM));
        gPower_vs_Lumi[iSiPM]   ->SetName(Form("gPower_vs_Lumi_%d", iSiPM));
        gCTR_vs_Lumi[iSiPM]     ->SetName(Form("gCTR_vs_Lumi_%d", iSiPM));
        gOV_vs_Lumi[iSiPM]      ->SetName(Form("gOV_vs_Lumi_%d", iSiPM));
        gSignal_vs_Lumi[iSiPM]  ->SetName(Form("gSignal_vs_Lumi_%d", iSiPM));
        gEleSignal_vs_Lumi[iSiPM]   ->SetName(Form("gEleSignal_vs_Lumi_%d", iSiPM));
        gBusyCells_vs_Lumi[iSiPM]   ->SetName(Form("gBusyCells_vs_Lumi_%d", iSiPM));
        
        
        gOptimalPDE_vs_Lumi[iSiPM]   ->SetName(Form("gOptimalPDE_vs_Lumi_%d", iSiPM));
        gOptimalDCR_vs_Lumi[iSiPM]   ->SetName(Form("gOptimalDCR_vs_Lumi_%d", iSiPM));
        gOptimalOV_vs_Lumi[iSiPM]    ->SetName(Form("gOptimalOV_vs_Lumi_%d", iSiPM));
//         gOptimalSNR_vs_Lumi[iSiPM]   ->SetName(Form("gOptimalSNR_vs_Lumi_%d", iSiPM));
        gOptimalCTR_vs_Lumi[iSiPM]   ->SetName(Form("gOptimalCTR_vs_Lumi_%d", iSiPM));
//         gOptimalCTR2x2_vs_Lumi[iSiPM]->SetName(Form("gOptimalCTR2x2_vs_Lumi_%d", iSiPM));
        gOptimalCTR_Double_vs_Lumi[iSiPM]   ->SetName(Form("gOptimalCTR_Double_vs_Lumi_%d", iSiPM));
//         gOptimalCTR_Double2x2_vs_Lumi[iSiPM]->SetName(Form("gOptimalCTR_Double2x2_vs_Lumi_%d", iSiPM));
//         gOptimalCTR_TP_vs_Lumi[iSiPM]       ->SetName(Form("gOptimalCTR_TP_vs_Lumi_%d", iSiPM));
//         gOptimalCTR_TP3x3_vs_Lumi[iSiPM]   ->SetName(Form("gOptimalCTR_TP3x3_vs_Lumi_%d", iSiPM));
//         gOptimalCTR_TP2x2_vs_Lumi[iSiPM]   ->SetName(Form("gOptimalCTR_TP2x2_vs_Lumi_%d", iSiPM));
        gOptimalCurrent_vs_Lumi[iSiPM]   ->SetName(Form("gOptimalCurrent_vs_Lumi_%d", iSiPM));
        gOptimalPower_vs_Lumi[iSiPM]    ->SetName(Form("gOptimalPower_vs_Lumi_%d", iSiPM));
        gOptimalSignal_vs_Lumi[iSiPM]   ->SetName(Form("gOptimalSignal_vs_Lumi_%d", iSiPM));
        gOptimalEleSignal_vs_Lumi[iSiPM] ->SetName(Form("gOptimalEleSignal_vs_Lumi_%d", iSiPM));
        gOptimalBusyCells_vs_Lumi[iSiPM] ->SetName(Form("gOptimalBusyCells_vs_Lumi_%d", iSiPM));
        
        gOptimalPhotJitter_vs_Lumi[iSiPM] ->SetName(Form("gOptimalPhotJitter_vs_Lumi_%d", iSiPM));
        gOptimalDCRJitter_vs_Lumi[iSiPM]  ->SetName(Form("gOptimalDCRJitter_vs_Lumi_%d", iSiPM));
        
        gOptimalPhotJitter_vs_Lumi_double[iSiPM] ->SetName(Form("gOptimalPhotJitter_vs_Lumi_double_%d", iSiPM));
        gOptimalDCRJitter_vs_Lumi_double[iSiPM]  ->SetName(Form("gOptimalDCRJitter_vs_Lumi_double_%d", iSiPM));
        
    
        gDCR_vs_OV[iSiPM]->Write();
//         gGain_vs_OV[iSiPM]->Write();
//         gPDE_vs_OV[iSiPM]->Write();
//         gCurrent_vs_OV[iSiPM]->Write();
    
    
        gPDE_vs_Lumi[iSiPM]     ->Write();
        gCurrent_vs_Lumi[iSiPM] ->Write();
        gGain_vs_Lumi[iSiPM]    ->Write();
        gDCR_vs_Lumi[iSiPM]     ->Write();
        gPower_vs_Lumi[iSiPM]   ->Write();
        gCTR_vs_Lumi[iSiPM]     ->Write();
        gOV_vs_Lumi[iSiPM]      ->Write();
        gSignal_vs_Lumi[iSiPM]  ->Write();
        gEleSignal_vs_Lumi[iSiPM]   ->Write();
        gBusyCells_vs_Lumi[iSiPM]   ->Write();
        
        
        gOptimalPDE_vs_Lumi[iSiPM]   ->Write();
        gOptimalDCR_vs_Lumi[iSiPM]   ->Write();
        gOptimalOV_vs_Lumi[iSiPM]    ->Write();
//         gOptimalSNR_vs_Lumi[iSiPM]   ->Write();
        gOptimalCTR_vs_Lumi[iSiPM]   ->Write();
//         gOptimalCTR2x2_vs_Lumi[iSiPM]->Write();
        gOptimalCTR_Double_vs_Lumi[iSiPM]   ->Write();
//         gOptimalCTR_Double2x2_vs_Lumi[iSiPM]->Write();
//         gOptimalCTR_TP_vs_Lumi[iSiPM]       ->Write();
//         gOptimalCTR_TP3x3_vs_Lumi[iSiPM]   ->Write();
//         gOptimalCTR_TP2x2_vs_Lumi[iSiPM]   ->Write();
        gOptimalCurrent_vs_Lumi[iSiPM]   ->Write();
        gOptimalPower_vs_Lumi[iSiPM]    ->Write();
        gOptimalSignal_vs_Lumi[iSiPM]   ->Write();
        gOptimalEleSignal_vs_Lumi[iSiPM] ->Write();
        gOptimalBusyCells_vs_Lumi[iSiPM] ->Write();
        
        gOptimalPhotJitter_vs_Lumi[iSiPM] ->Write();
        gOptimalDCRJitter_vs_Lumi[iSiPM]  ->Write();
        
        gOptimalPhotJitter_vs_Lumi_double[iSiPM] ->Write();
        gOptimalDCRJitter_vs_Lumi_double[iSiPM]  ->Write();
    }
    
    
    gOptimalCTR_Double_vs_Lumi_Ave->SetName("gOptimalCTR_Double_vs_Lumi_Ave");
    gOptimalCTR_Double_vs_Lumi_Ave->Write();
    
//     gOptimalCTR_TP_vs_Lumi_Ave->SetName("gOptimalCTR_TP_vs_Lumi_Ave");
//     gOptimalCTR_TP_vs_Lumi_Ave->Write();
    
    gOptimalCTR_vs_Lumi_Ave->SetName("gOptimalCTR_vs_Lumi_Ave");
    gOptimalCTR_vs_Lumi_Ave->Write();
    
    fileOutput->Close();

    
    
}
