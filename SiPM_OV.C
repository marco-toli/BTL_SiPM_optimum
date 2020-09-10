void SiPM_OV()
{
    
    
    gStyle->SetTitleXOffset (1.00) ;                                                                                        
    gStyle->SetTitleYOffset (1.2) ;                                                                                                                                                                                                                 
    gStyle->SetPadLeftMargin (0.13) ;                                                                                       
    gStyle->SetPadBottomMargin (0.13) ;                                                                                                                                                                                                              
    gStyle->SetTitleSize (0.05, "xyz") ;                                                                                    
    gStyle->SetLabelSize (0.035,"xyz") ;  
        
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(42);
    gStyle->SetLegendTextSize(0.035);
    TLegend * leg;
  
  
    const int NSIPM = 8;
    
    std::string sipm_name[NSIPM];
    //TDR technologies
    sipm_name[0] = "S12572-015";
    sipm_name[1] = "HDR2-015-TDR";
    sipm_name[2] = "FBK-TE-TDR";
    //post-TDR technologies
    sipm_name[3] = "HDR2-015-v2";    
    sipm_name[4] = "FBK-W2C";
    sipm_name[5] = "FBK-W4C";
    sipm_name[6] = "FBK-W7C";
    sipm_name[7] = "FBK-W9C";
        
    //pre-series technologies


    std::string output_folder = "./temp/";
//     TFile * fileOutput = new TFile ("./temp/temp.root", "RECREATE");    
    TFile * fileOutput = new TFile ("./output/new_data_full_fluence.root", "RECREATE");    
//     TFile * fileOutput = new TFile ("./output/default_TDR.root", "RECREATE");    

        
    string outdaq;
    const char * daqfile;       
    
    //inputs from command line
    float my_temp   = -40;    
    float sipm_area = 9; //mm²   
    

    //define inputs
    int TOTLUMI = 6000;
        
//     float MIP_rate = 2.5e6;   //MHz  
    float MIP_rate = 0;   //MHz  
    float maxPower      = 50;    //mW / channel    
    int nChannels       = 500000; //total BTL number of channels
    
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
    
    float ref_temp = -30;   //temperature to which the data in the root files refer to
    float fluence_effects_norm = 3160;// luminosity corresponding to 2e14 neq/cm²

    float temp_coefficient[NSIPM];

    // import SiPM data
        
    TFile * input_file_SiPM[NSIPM];
    
    float data_Vbr[NSIPM];
    float data_Vbr_rind_drift[NSIPM];
    float data_Vbr_t_coeff[NSIPM];
    float data_dcr_t_coeff[NSIPM];
    float data_temperature[NSIPM];
    float data_spad_size[NSIPM];   
    float data_spad_RC[NSIPM];
    float data_QE_loss[NSIPM]; //Quantum Efficiency loss after 2e14
    float data_gain_loss[NSIPM]; //residual gain after 2e14        
    
    TF1 * fPDE_vs_WL[NSIPM];
    TF1 * fPDE_vs_OV[NSIPM];
    TF1 * fCurrent_vs_OV[NSIPM];
    TF1 * fGain_vs_OV[NSIPM];
    TF1 * fENF_vs_OV[NSIPM];
    
    TGraphErrors * gDCR_vs_OV[NSIPM];
        
    for (int iSiPM = 0; iSiPM<NSIPM; iSiPM++) 
    {
        std::cout << "********************************************************** " << std::endl;
        std::cout << "reading input specifications for SiPM: " << sipm_name[iSiPM].c_str() << std::endl;
        
        
        if (iSiPM <3)      input_file_SiPM[iSiPM] = new TFile ( Form("sipm_spec_input_%s.root", sipm_name[iSiPM].c_str()), "READ");                     
        else if (iSiPM == 4) input_file_SiPM[iSiPM] = new TFile ( Form("sipm_spec_input_%s-1e13.root", sipm_name[iSiPM].c_str()), "READ");             
//         else if (iSiPM <8) input_file_SiPM[iSiPM] = new TFile ( Form("sipm_spec_input_%s-1e13.root", sipm_name[iSiPM].c_str()), "READ");             
        else               input_file_SiPM[iSiPM] = new TFile ( Form("sipm_spec_input_%s-2e14.root", sipm_name[iSiPM].c_str()), "READ");             
        
        fPDE_vs_OV[iSiPM]     = (TF1*) input_file_SiPM[iSiPM]->Get("fPDE_LYSO_vs_OV");        
        fPDE_vs_WL[iSiPM]     = (TF1*) input_file_SiPM[iSiPM]->Get("fPDE_vs_WL");
        fCurrent_vs_OV[iSiPM] = (TF1*) input_file_SiPM[iSiPM]->Get("fCurrent_vs_OV");
        fGain_vs_OV[iSiPM]    = (TF1*) input_file_SiPM[iSiPM]->Get("fGain_vs_OV");
        fENF_vs_OV[iSiPM]     = (TF1*) input_file_SiPM[iSiPM]->Get("fENF_vs_OV");
        gDCR_vs_OV[iSiPM]     = (TGraphErrors*) input_file_SiPM[iSiPM]->Get("gDCR_vs_OV");

                        
        //floats
        TGraphErrors * temp;
        
        temp = (TGraphErrors*) input_file_SiPM[iSiPM]->Get("data_dcr_t_coeff");        
        data_dcr_t_coeff[iSiPM] = temp->Eval(1);
        std::cout << "DCR temp. coefficient is: " << data_dcr_t_coeff[iSiPM] << std::endl;
        
        temp = (TGraphErrors*) input_file_SiPM[iSiPM]->Get("data_Vbr_t_coeff");        
        data_Vbr_t_coeff[iSiPM] = temp->Eval(1);
        std::cout << "Vbr temp. coefficient is: " << data_Vbr_t_coeff[iSiPM] << std::endl;
        
        temp = (TGraphErrors*) input_file_SiPM[iSiPM]->Get("data_QE_loss");        
        data_QE_loss[iSiPM] = temp->Eval(1);
        std::cout << "Quantum eff. loss after 2e14 is: " << data_QE_loss[iSiPM] << std::endl;
        
                
        temp = (TGraphErrors*) input_file_SiPM[iSiPM]->Get("data_gain_loss");        
        data_gain_loss[iSiPM] = temp->Eval(1);
        std::cout << "Residual gain after 2e14 is: " << data_gain_loss[iSiPM] << std::endl;
        
        temp = (TGraphErrors*) input_file_SiPM[iSiPM]->Get("data_Vbr_rind_drift");        
        data_Vbr_rind_drift[iSiPM] = temp->Eval(1);
        std::cout << "Vbr drift after 2e14: " << data_Vbr_rind_drift[iSiPM] << std::endl;
        
        
        temp = (TGraphErrors*) input_file_SiPM[iSiPM]->Get("data_spad_size");        
        data_spad_size[iSiPM] = temp->Eval(1);
        std::cout << "spad size is: " << data_spad_size[iSiPM] << std::endl;
        
        temp = (TGraphErrors*) input_file_SiPM[iSiPM]->Get("data_Vbr");        
        data_Vbr[iSiPM] = temp->Eval(1);
        std::cout << "Vbr at -30 is: " << data_Vbr[iSiPM] << std::endl;
        
        temp = (TGraphErrors*) input_file_SiPM[iSiPM]->Get("data_spad_RC");        
        data_spad_RC[iSiPM] = temp->Eval(1);   
        std::cout << "spad RC is: " << data_spad_RC[iSiPM] << std::endl;
        
        
                
        //adjust some values to input temperature: e.g. Vbr
        std::cout << " breakdown voltage (" << ref_temp << "deg C)[" << iSiPM << "] = "  << data_Vbr[iSiPM] << std::endl;
        
        data_Vbr[iSiPM] += data_Vbr_t_coeff[iSiPM] * (my_temp-ref_temp);        
        std::cout << " breakdown voltage (" << my_temp << "deg C)[" << iSiPM << "] = "  << data_Vbr[iSiPM] << std::endl;
        
        //hack some values
//         data_QE_loss[iSiPM] = 0;    // e.g. ignore QE losses
        
        //normalize current if needed
//         float normal_current = 200/2.1/4000./9./anneal_coeff/insitu_recovery*newFLUKA_fluence;
        temp_coefficient[iSiPM] = 1./pow(data_dcr_t_coeff[iSiPM], (ref_temp-my_temp)/10.);
        
        std::cout << " temp_coefficient[" << iSiPM << "] = "  << temp_coefficient[iSiPM] << std::endl;
        
                //saving graphs and functions   
        
    }
    
    
    //get DCR annealing scenarios
    
    std::cout << "********************************************************** " << std::endl;
    std::cout << "reading input annealing scenario... "<< std::endl;
        
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
    
    

    
    //defining graphs for extrapolations
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
    TGraphErrors * gOptimalCTR_Double_vs_Lumi[NSIPM];

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
    
        
    TGraphErrors * gCurrent_vs_Temperature[NSIPM];
    TGraphErrors * gDCR_vs_Temperature[NSIPM];
    
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
        gOptimalCTR_Double_vs_Lumi[iSiPM]   = new TGraphErrors ();
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
        
        gCurrent_vs_Temperature[iSiPM] = new TGraphErrors();
        gDCR_vs_Temperature[iSiPM] = new TGraphErrors();
    }
    
    TGraphAsymmErrors * gOptimalCTR_vs_Lumi_Ave        = new TGraphAsymmErrors ();
    TGraphAsymmErrors * gOptimalCTR_Double_vs_Lumi_Ave = new TGraphAsymmErrors ();
    
    TGraphErrors * gDigitization = new TGraphErrors ();
    TGraphErrors * gClock        = new TGraphErrors ();
    TGraphErrors * gElectronics  = new TGraphErrors ();
    
    TGraphErrors * gDigitization_double = new TGraphErrors ();
    TGraphErrors * gClock_double        = new TGraphErrors ();
    TGraphErrors * gElectronics_double  = new TGraphErrors ();
    
    
    
    std::cout << "********************************************************** " << std::endl;
    std::cout << "starting lumi scan... "<< std::endl;
    
    for (int iLumi = 0; iLumi < TOTLUMI; iLumi++)
    {
    

//         float ave_CTR_TP = 0;
        
        float min_CTR_Single = 99999;
        float max_CTR_Single = -9999;
                
        float min_CTR_Double = 99999;
        float max_CTR_Double = -9999;
        
        gDigitization->SetPoint(iLumi, iLumi, sigma_digi);
        gClock       ->SetPoint(iLumi, iLumi, sigma_clock);
        gElectronics ->SetPoint(iLumi, iLumi, sigma_elect);
        
        gDigitization_double->SetPoint(iLumi, iLumi, sigma_digi/sqrt(2));
        gClock_double       ->SetPoint(iLumi, iLumi, sigma_clock);
        gElectronics_double ->SetPoint(iLumi, iLumi, sigma_elect/sqrt(2));
        
        for (int iSiPM = 0; iSiPM< NSIPM; iSiPM++)
        {
            
            //to include annealing model and HL-LHC scenario
//             float current = fCurrent_vs_OV[iSiPM]->Eval(bias_ov) * gLumi_to_DCR_annealing->Eval(iLumi) / gLumi_to_DCR_noannealing->Eval(iLumi) * (iLumi*sipm_area *2.6);            //in uA
            
            //for fixed DCR increase
            float current = iLumi*fCurrent_vs_OV[iSiPM]->Eval(bias_ov)*sipm_area*temp_coefficient[iSiPM];            //in uA --
            float DCR     = current/1e6 / (1.6e-19) / fGain_vs_OV[iSiPM]->Eval(bias_ov)/1e9 / (1-iLumi/fluence_effects_norm*(1-data_gain_loss[iSiPM]));  //current from uA to A
            
            double busy_cells = DCR/sipm_area*data_spad_RC[iSiPM]*2./pow(1000./data_spad_size[iSiPM],2);
            if (busy_cells >1) busy_cells = 1.;
                
            double signal = fPDE_vs_OV[iSiPM]->Eval(bias_ov)/100*(1.-busy_cells/1.) * LO * (1-data_QE_loss[iSiPM]/fluence_effects_norm*iLumi);
            double sigma_phot = sigma_phot_ref  *sqrt(LO*PDE_at_TB/signal); 
            double sigma_DCR = sigma_DCR_20*sqrt(DCR/DCR_ref)*pow(Nphe_DCR_ref/signal, DCR_alpha);
            
            double CTR = sqrt(pow(sigma_phot,2) + pow(sigma_DCR,2)  + pow(sigma_digi,2) +  pow(sigma_elect,2) +  pow(sigma_clock,2));
            
            float power = current*(bias_ov+data_Vbr[iSiPM])/1e3; //converting current from uA to mA
            float tot_power = power*nChannels;
            
            double ele_signal = signal*fGain_vs_OV[iSiPM]->Eval(bias_ov)*(1-iLumi/fluence_effects_norm*(1-data_gain_loss[iSiPM]));
            float power_dyn   = MIP_rate*ele_signal*1.602*1e-19*(bias_ov+data_Vbr[iSiPM])*1e3; //converting current from Ampere to mA
            
            gCurrent_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, current);
            gGain_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, fGain_vs_OV[iSiPM]->Eval(bias_ov)*(1-iLumi/fluence_effects_norm*(1-data_gain_loss[iSiPM])));
            gPower_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, power);
            gPowerDynamic_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, power_dyn);
            gPowerTot_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, power_dyn+power);
            gDCR_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, DCR);
            
            gPDE_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, fPDE_vs_OV[iSiPM]->Eval(bias_ov)/100*(1-data_QE_loss[iSiPM]/fluence_effects_norm*iLumi));
            gOV_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, bias_ov);                        
            gCTR_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, CTR);
            
            gSignal_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, signal);
            gEleSignal_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, signal*fGain_vs_OV[iSiPM]->Eval(bias_ov)*(1-iLumi/fluence_effects_norm*(1-data_gain_loss[iSiPM])));
            gBusyCells_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, busy_cells);
            

            
            float best_SNR = -9999999;
            
            float best_sigma_phot = 99999;
            float best_sigma_DCR = 99999;
            
            float best_CTR  = 99999;
            float best_CTR_double  = 99999;
            
            float best_signal = 0;
            float best_ele_signal = 0;
            float best_busy_cell = 0;
            float opt_current;
            float opt_power;
            float opt_dyn_power;
            float opt_ov;
            float opt_PDE;
            float opt_DCR;
            
//             std::cout << "starting OV scan for lumi: " << iLumi << " ..." << std::endl;
            
            for (int iOV = 0; iOV<NOVS; iOV++)
            {
//                 std::cout << "testing OV = " << iOV*ov_step << std::endl;
                //calculate SNR at this OV and this lumi 
                float temp_ov      = min_ov+iOV*ov_step;
                float temp_PDE     = fPDE_vs_OV[iSiPM]->Eval(temp_ov)/100 * (1-data_QE_loss[iSiPM]/fluence_effects_norm*iLumi);
                
            
                //for fixed DCR increase            
                float temp_current     = iLumi*fCurrent_vs_OV[iSiPM]->Eval(temp_ov)*sipm_area*temp_coefficient[iSiPM];                
                float temp_DCR         = temp_current/1e6 / (1.6e-19) / fGain_vs_OV[iSiPM]->Eval(temp_ov)/1e9/(1-iLumi/fluence_effects_norm*(1-data_gain_loss[iSiPM])) ;                                
                float temp_power       = temp_current*(temp_ov+data_Vbr[iSiPM])/1e3; //converting current to Ampere
                
            
                //estimating CTR
                double temp_busy_cells = temp_DCR/sipm_area*data_spad_RC[iSiPM]*2./pow(1000./data_spad_size[iSiPM],2);
                if (temp_busy_cells >1) busy_cells = 1.;
                
                double temp_signal     = temp_PDE*(1.-temp_busy_cells/1.) * LO;
                double temp_ele_signal = temp_signal*fGain_vs_OV[iSiPM]->Eval(temp_ov)*(1-iLumi/fluence_effects_norm*(1-data_gain_loss[iSiPM]));
                double temp_dyn_power  = MIP_rate*temp_ele_signal*1.6e-19*(bias_ov+data_Vbr[iSiPM])*1e3; //converting current to Ampere
                
                double sigma_phot      = sigma_phot_ref*sqrt(LO*PDE_at_TB/temp_signal);       
                double sigma_DCR       = sigma_DCR_20*sqrt(temp_DCR/DCR_ref)*pow(Nphe_DCR_ref/temp_signal, DCR_alpha);                             
                double temp_CTR        = sqrt(pow(sigma_phot,2) + pow(sigma_DCR,2)  + pow(sigma_digi,2) +  pow(sigma_elect,2) +  pow(sigma_clock,2));
             
                
                if (true
                    && temp_CTR < best_CTR                    
                    && temp_power<maxPower
//                     && (temp_power+temp_dyn_power)<maxPower
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
                
            }
           
            if (best_CTR < min_CTR_Single)              min_CTR_Single = best_CTR;
            if (best_CTR_double < min_CTR_Double)       min_CTR_Double = best_CTR_double;

            if (best_CTR > max_CTR_Single)              max_CTR_Single = best_CTR;
            if (best_CTR_double > max_CTR_Double)       max_CTR_Double = best_CTR_double;
                        
//             std::cout << "min_CTR_TP = " << min_CTR_TP << " :: min_CTR_Single =  " << min_CTR_Single << " :: min_CTR_Double = " << min_CTR_Double << std::endl;
//             std::cout << "max_CTR_TP = " << max_CTR_TP << " :: max_CTR_Single =  " << max_CTR_Single << " :: max_CTR_Double = " << max_CTR_Double << std::endl;
            
            
            gOptimalPDE_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, opt_PDE);
            gOptimalDCR_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, opt_DCR);
            gOptimalOV_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, opt_ov);
            gOptimalGain_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, fGain_vs_OV[iSiPM]->Eval(opt_ov)*(1-iLumi/fluence_effects_norm*(1-data_gain_loss[iSiPM])));
            gOptimalCTR_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, best_CTR);
            gOptimalCTR_Double_vs_Lumi[iSiPM]->SetPoint(iLumi, iLumi, best_CTR_double);
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
                    
//                     std::cout << " temperature scan for T = " << temperature << " :: giving a coeff. in DCR reduction of = " << pow(data_dcr_t_coeff[iSiPM], (ref_temp-temperature)/10) << std::endl;
                    
                    gCurrent_vs_Temperature[iSiPM]->SetPoint(iT, temperature, opt_current/pow(data_dcr_t_coeff[iSiPM], (ref_temp-temperature)/10));
                    gDCR_vs_Temperature[iSiPM]->SetPoint(iT, temperature, opt_DCR/pow(data_dcr_t_coeff[iSiPM], (ref_temp-temperature)/10));
                    
                }
            }
        }//end of SiPM loop
        
        //to fill with best time resolution
        gOptimalCTR_vs_Lumi_Ave         ->SetPoint(iLumi, iLumi, min_CTR_Single);
        gOptimalCTR_Double_vs_Lumi_Ave  ->SetPoint(iLumi, iLumi, min_CTR_Double);
        
        gOptimalCTR_vs_Lumi_Ave         ->SetPointError(iLumi, 0, 0, 0, max_CTR_Single-min_CTR_Single);
        gOptimalCTR_Double_vs_Lumi_Ave  ->SetPointError(iLumi, 0, 0, 0, max_CTR_Double-min_CTR_Double);
        
    }
    
    
    
    std::cout << "********************************************************** " << std::endl;
    std::cout << "setting graphical cosmetics... "<< std::endl;

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
    
    for (int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        color_constant[iSiPM]  = (int) kRainBow+iSiPM*6;    
        color_optimized[iSiPM] = (int) kRainBow+iSiPM*6;
    }
    
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
    gCurrent_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gCurrent_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 3500);
    
    leg = new TLegend(0.16,0.64,0.65,0.88,NULL,"brNDC");
    
    for (int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        gCurrent_vs_Lumi[iSiPM]->Draw("same LE");    
        gOptimalCurrent_vs_Lumi[iSiPM]->Draw("same LE");
//         leg ->AddEntry(gCurrent_vs_Lumi[iSiPM], Form("constant bias = +1.5V OV (%s)", sipm_name[iSiPM].c_str()), "lp");
        leg ->AddEntry(gOptimalCurrent_vs_Lumi[iSiPM], Form("optimized bias (%s)", sipm_name[iSiPM].c_str()), "lp");
    }
    gPad->SetGridy();
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
    gPower_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gPower_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 160);
    for (int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        gPower_vs_Lumi[iSiPM]->Draw("same LE");    
        gOptimalPower_vs_Lumi[iSiPM]->Draw("same LE");
    }
    
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
    gPowerDynamic_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gPowerDynamic_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 180);
    for (int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        gPowerDynamic_vs_Lumi[iSiPM]->Draw("same LE");    
        gOptimalPowerDynamic_vs_Lumi[iSiPM]->Draw("same LE");
    }
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
    gPowerTot_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gPowerTot_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 180);
    for (int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        gPowerTot_vs_Lumi[iSiPM]->Draw("same LE");    
        gOptimalPowerTot_vs_Lumi[iSiPM]->Draw("same LE");
    }
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
    gOptimalOV_vs_Lumi[0]->GetYaxis()->SetTitleOffset(1.1);
    gOptimalOV_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gOptimalOV_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 6);
    for (int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        gOV_vs_Lumi[iSiPM]->Draw("same LE");
        gOptimalOV_vs_Lumi[iSiPM]->Draw("same LE");            
    }
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
    gOptimalPDE_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gOptimalPDE_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 0.6);
    for (int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        gPDE_vs_Lumi[iSiPM]->Draw("same LE");
        gOptimalPDE_vs_Lumi[iSiPM]->Draw("same LE");            
    }
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
    gOptimalSignal_vs_Lumi[0]->GetXaxis()->SetTitleOffset(0.85);
    gOptimalSignal_vs_Lumi[0]->GetYaxis()->SetTitleOffset(1.15);
    gOptimalSignal_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gOptimalSignal_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 20000);
    for (int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        gSignal_vs_Lumi[iSiPM]->Draw("same LE");
        gOptimalSignal_vs_Lumi[iSiPM]->Draw("same LE");            
    }
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
    gOptimalEleSignal_vs_Lumi[0]->GetXaxis()->SetTitleOffset(0.85);
    gOptimalEleSignal_vs_Lumi[0]->GetYaxis()->SetTitleOffset(1.1);
    gOptimalEleSignal_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gOptimalEleSignal_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 10000e6);
    for (int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        gEleSignal_vs_Lumi[iSiPM]->Draw("same LE");
        gOptimalEleSignal_vs_Lumi[iSiPM]->Draw("same LE");            
    }
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
    gOptimalBusyCells_vs_Lumi[0]->GetXaxis()->SetTitleOffset(0.85);
    gOptimalBusyCells_vs_Lumi[0]->GetYaxis()->SetTitleOffset(1.1);
    gOptimalBusyCells_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gOptimalBusyCells_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 0.06);
    for (int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        gBusyCells_vs_Lumi[iSiPM]->Draw("same LE");
        gOptimalBusyCells_vs_Lumi[iSiPM]->Draw("same LE");            
    }
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
    gDCR_vs_Lumi[0]->GetXaxis()->SetTitleOffset(0.85);
    gDCR_vs_Lumi[0]->GetYaxis()->SetTitleOffset(1.1);
    gDCR_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gDCR_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 120);
    for (int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        gDCR_vs_Lumi[iSiPM]->Draw("same LE");
        gOptimalDCR_vs_Lumi[iSiPM]->Draw("same LE");            
    }
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
    gGain_vs_Lumi[0]->GetXaxis()->SetTitleOffset(0.85);
    gGain_vs_Lumi[0]->GetYaxis()->SetTitleOffset(1.1);
    gGain_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gGain_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 8e5);
    gGain_vs_Lumi[1]->Draw("same LE");
    for (int iSiPM = 0; iSiPM < NSIPM; iSiPM++)
    {
        gGain_vs_Lumi[iSiPM]->Draw("same LE");
        gOptimalGain_vs_Lumi[iSiPM]->Draw("same LE");            
    }
    gPad->SetGridy();
    leg->Draw();
    
    cGain_vs_Lumi->cd();
    outdaq = Form("%sGain_vs_Lumi.png", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cGain_vs_Lumi->SaveAs(daqfile);
    outdaq = Form("%sGain_vs_Lumi.pdf", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cGain_vs_Lumi->SaveAs(daqfile);
    
    /*
    
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
    
    */

    
    TLegend * leg2;
    
    
    //PLOTS for CDR/TDR
    /*
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
    */
    
    
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
        gCurrent_vs_Temperature[iSiPM]->SetLineColor(color_optimized[iSiPM]);
        gCurrent_vs_Temperature[iSiPM]->SetMarkerColor(color_optimized[iSiPM]);
        gCurrent_vs_Temperature[iSiPM]->SetLineWidth(2);
        
        leg2->AddEntry(gCurrent_vs_Temperature[iSiPM], Form("%s, T_{coeff} = %.2f", sipm_name[iSiPM].c_str(), data_dcr_t_coeff[iSiPM]), "lp");
    }
//     gCurrent_vs_Temperature[1]->SetLineColor(kGreen+2);
//     gCurrent_vs_Temperature[1]->SetMarkerColor(kGreen+2);
//     gCurrent_vs_Temperature[2]->SetLineColor(kRed+1);
//     gCurrent_vs_Temperature[2]->SetMarkerColor(kRed+1);
        
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
        gDCR_vs_Temperature[iSiPM]->SetLineColor(color_optimized[iSiPM]);
        gDCR_vs_Temperature[iSiPM]->SetMarkerColor(color_optimized[iSiPM]);
        gDCR_vs_Temperature[iSiPM]->SetLineWidth(2);        
    }
//     gDCR_vs_Temperature[1]->SetLineColor(kGreen+2);
//     gDCR_vs_Temperature[1]->SetMarkerColor(kGreen+2);
//     gDCR_vs_Temperature[2]->SetLineColor(kRed+1);
//     gDCR_vs_Temperature[2]->SetMarkerColor(kRed+1);
        
    leg2->Draw();
    gPad->SetLogy();
    
    
    
        
    TCanvas * cCTR_vs_Lumi = new TCanvas ("cCTR_vs_Lumi", "cCTR_vs_Lumi", 600, 500);    
    cCTR_vs_Lumi->SetLeftMargin(0.12);
    gOptimalCTR_Double_vs_Lumi[0]->GetXaxis()->SetTitle("Integrated Luminosity [fb^{-1}]");
    gOptimalCTR_Double_vs_Lumi[0]->GetYaxis()->SetTitle("Time resolution [ps]");
    gOptimalCTR_Double_vs_Lumi[0]->GetXaxis()->SetTitleSize(0.05);
    gOptimalCTR_Double_vs_Lumi[0]->GetYaxis()->SetTitleSize(0.05);
    gOptimalCTR_Double_vs_Lumi[0]->GetXaxis()->SetTitleOffset(0.85);
    gOptimalCTR_Double_vs_Lumi[0]->GetYaxis()->SetRangeUser(0, 140);
    gOptimalCTR_Double_vs_Lumi[0]->GetXaxis()->SetRangeUser(0, 4000);
    gOptimalCTR_Double_vs_Lumi[0]->SetMaximum(140);
//     gOptimalCTR_Double_vs_Lumi[0]->SetLineColor(kGreen+1);
//     gOptimalCTR_Double_vs_Lumi[0]->SetMarkerColor(kGreen+1);
    gOptimalCTR_Double_vs_Lumi[0]->Draw("ALE");
    
// //     gCTR_vs_Lumi[0]->Draw("same LPE");
// //     gCTR_vs_Lumi[1]->Draw("same LPE");
// //     gOptimalCTR_vs_Lumi[1]->Draw("same LE");
// //     gOptimalCTR_vs_Lumi[2]->Draw("same LE");
//     
//     gOptimalCTR_Double_vs_Lumi[2]->Draw("same LE"); 
//     
//     gOptimalCTR_Double_vs_Lumi[0]->SetLineStyle(7);
//     gOptimalCTR_Double_vs_Lumi[0]->SetLineColor(kYellow+2);
//     gOptimalCTR_Double_vs_Lumi[0]->SetMarkerColor(kYellow+2);
//     gOptimalCTR_Double_vs_Lumi[0]->Draw("same LE");
//     
//     
//     gOptimalCTR_Double_vs_Lumi[1]->SetLineStyle(5);
//     gOptimalCTR_Double_vs_Lumi[1]->SetLineColor(kOrange+2);
//     gOptimalCTR_Double_vs_Lumi[1]->SetMarkerColor(kOrange+2);
//     gOptimalCTR_Double_vs_Lumi[1]->Draw("same LE");
    gPad->SetGridy();
    
    
    
    leg2 = new TLegend(0.15,0.7,0.55,0.88,NULL,"brNDC");
 
    for (int iSiPM = 0; iSiPM< NSIPM; iSiPM++)
    {
//         gOptimalCTR_Double_vs_Lumi[0]->SetLineStyle(7);
        gOptimalCTR_Double_vs_Lumi[iSiPM]->Draw("same LE");
        gOptimalCTR_Double_vs_Lumi[iSiPM]->SetLineColor(color_optimized[iSiPM]);
        gOptimalCTR_Double_vs_Lumi[iSiPM]->SetMarkerColor(color_optimized[iSiPM]);
        leg2->AddEntry(gOptimalCTR_Double_vs_Lumi[iSiPM], Form("SiPM type: (%s)", sipm_name[iSiPM].c_str()), "lp");
//         leg2->AddEntry(gOptimalCTR_vs_Lumi[iSiPM], Form("SiPM type: (%s) [single read-out]", sipm_name[iSiPM].c_str()), "lp");
    }
    leg2->Draw();
    
    cCTR_vs_Lumi->cd();
    outdaq = Form("%sTimeRes_vs_Lumi.png", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cCTR_vs_Lumi->SaveAs(daqfile);
    outdaq = Form("%sTimeRes_vs_Lumi.pdf", output_folder.c_str());
    daqfile = outdaq.c_str();    
    cCTR_vs_Lumi->SaveAs(daqfile);
    
    
    
    
    
    fileOutput->cd();
    for (int iSiPM = 0; iSiPM<NSIPM; iSiPM++)
    {
        fCurrent_vs_OV[iSiPM]->Write();
        fENF_vs_OV[iSiPM]->Write();
        fGain_vs_OV[iSiPM]->Write();
        fPDE_vs_OV[iSiPM]->Write();
        fPDE_vs_WL[iSiPM]->Write();
        
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
        gOptimalCTR_vs_Lumi[iSiPM]   ->SetName(Form("gOptimalCTR_vs_Lumi_%d", iSiPM));
        gOptimalCTR_Double_vs_Lumi[iSiPM]   ->SetName(Form("gOptimalCTR_Double_vs_Lumi_%d", iSiPM));
        gOptimalCurrent_vs_Lumi[iSiPM]   ->SetName(Form("gOptimalCurrent_vs_Lumi_%d", iSiPM));
        gOptimalPower_vs_Lumi[iSiPM]    ->SetName(Form("gOptimalPower_vs_Lumi_%d", iSiPM));
        gOptimalSignal_vs_Lumi[iSiPM]   ->SetName(Form("gOptimalSignal_vs_Lumi_%d", iSiPM));
        gOptimalEleSignal_vs_Lumi[iSiPM] ->SetName(Form("gOptimalEleSignal_vs_Lumi_%d", iSiPM));
        gOptimalBusyCells_vs_Lumi[iSiPM] ->SetName(Form("gOptimalBusyCells_vs_Lumi_%d", iSiPM));
        
        gOptimalPhotJitter_vs_Lumi[iSiPM] ->SetName(Form("gOptimalPhotJitter_vs_Lumi_%d", iSiPM));
        gOptimalDCRJitter_vs_Lumi[iSiPM]  ->SetName(Form("gOptimalDCRJitter_vs_Lumi_%d", iSiPM));
        
        gOptimalPhotJitter_vs_Lumi_double[iSiPM] ->SetName(Form("gOptimalPhotJitter_vs_Lumi_double_%d", iSiPM));
        gOptimalDCRJitter_vs_Lumi_double[iSiPM]  ->SetName(Form("gOptimalDCRJitter_vs_Lumi_double_%d", iSiPM));
        
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
        gOptimalCTR_vs_Lumi[iSiPM]   ->Write();
        gOptimalCTR_Double_vs_Lumi[iSiPM]   ->Write();
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
    
    gOptimalCTR_vs_Lumi_Ave->SetName("gOptimalCTR_vs_Lumi_Ave");
    gOptimalCTR_vs_Lumi_Ave->Write();
    
    fileOutput->Close();

    
    
}
