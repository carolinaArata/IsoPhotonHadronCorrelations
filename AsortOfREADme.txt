
Before running analysis and plot, ctrl+f and check the directories of these files in IsoGammaHadron.C:
  - TFile *fPurity
  - TFile *fileData 
  - TFile *fileDataShStd 
  - TFile *fileDataShBkg 
  - TFile *fileDataMix 
  - TFile *fileMC[nShSh]

  In the scripts and in the codes CTRL+F for ~/work/histogram/ or ~/wor/histogram/FromScratch/, set yours correct folders.

  The analysis is structured in three main blocks:
    - analysis: RunAnalysis.sh
    - combine the 0-10% and the 10-30% centralities: RunCombine.sh
    - computation of all the systematics: systematic_sw/RunSystematics.sh

General macro: IsoGammaHadron.C; PlotIsoGammaHadron.C

1) Run Analysis and plots with: 

   RunAnalysis.sh

   This script contains the analysis and it also produces all the root files for the systematic studies 

    * Purity systematic is obtained from data, using the systematic contained in the purity.root. Setting systPur = 0.9 or 1.1, we define if we add the systematic (systPur == 1.1) or we subtract (systPur == 0.9) the systematic from the purity value. To do this there is an if condition to add or subtract the systematic based on the systPur value, then the **relative systematic error** is taken and included in the purity value. The purity value is obtained from the fit. (lines 367-393)

    * ShSh Background systematic is obtained using different files for the ShSh background: 
      - ShShSystMultBkg contains only Iso1 with multi Background: 0.40-2.00, 0.50-2.00, 0.60-2.00, 0.40-1.50 
      - ShShSyst03510 contains only Iso 1 and background 0.35-1.00
      
      Since these trains were run for ISOLATION ONLY, when the systematics on the ShShBkg is selected (systShSh = true), there is a condition to use standard background (0.40-1.00) from normal data if (iSh == 1 && iso == 0 && systShSh == true) and to access to root files with different background if (iSh == 1 && iso == 1 && systShSh == true). (lines 133-146)
      
      To run these studies on different shShBkg it is suffisant to change the shshBkg (TString shshBkg = "xx.xx-xx.xx") to the corresponding range and set systShSh = true

    * Tracking Efficiency 
      Set systTrackIneff = true : this will allow to run the analysis using the MC file contained in MC_GJShSh150/MC_GJTrackInEff.root;

    * Number of centrality bins used for the mixed event
      Set systNMix18 = true and systNMix45 = false; then do the opposite systNMix18 = false and systNMix45 = true
      There are two different booleans for selecting the number of centrality bins used 18 and 45.
      
      The analysis run with old files NCentrBinMix18 and NCentrBinMix45; there is single file.rot containing all centralities: EMCAL_MB0_90.root; the ShShBkg = 0.40-2.00 and the pThadr > 200 GeV;
      The code contains conditions in the way that when it is reading fileMix and the ishsh=1 (iSh == 1 && systNMix), the ShShBkg is 0.40-2.00. (line 156 and following) 
      
      The shshbkg need to be defined by hand in the code this can be set where the root files are declared before calling the Exec(...) : sShShNCentMix = "_ShSh0.40-2.00" for old files, sShShNCentMix = "_ShSh0.40-1.00" for new file.

      The last train with N Cent Mix generate strange value in the systematics, not expected 

      The pT trig distribution for the old files were saved with ShShBkg: "_ShSh0.40-2.00" for mixed only, the name is updated in line 173-174 with 
      if (systNMix)
        hTriggerMix[iso][iSh]->SetName(sHistName + sIso + sShSh + sCent + "_hPtTriggerMixed");
      The output histograms are all saved with _ShSh0.40-1.00 in the name, to be able to run all the macros and plotting with the others files results.

    * ZYAM, set bZYAM = true

2) Combine 0-10% and 10-30% centrality bins
   
   Combine0_30.sh, this script combines all the previous root files

3) Systematics

   For running systematics, use the script:
  
   RunSystematics.sh

   Systematics:
     - Purity:
       The code for estimating the systematic on the purity is in systematic_sw. 
       The reference files are saved in checkCode, while the systematics on purity are saved in checkCodeSystPur09, checkCodeSystPur11
       The directory of the files used for estimating the purity is passed to the function as  TString sDirRefFiles = "~/work/histogram/FromScratch/checkCode" and the last part of the name of the directory is defined inside the macro with: TString PurUsed[] = {"", "SystPur11", "SystPur09"};

     - ShSh Background: 
       The code for estimating the systematic on the ShSh background is in systematic_sw.

     - N Cent Mix: SystematicsNCentrBin (with NCent45)
       The systematic is calculated with a macro in systematic_sw. The files are defined before the macro in TString CentMixUsed[] = {"~/work/histogram/FromScratch/checkCode" /*Reference*/, "~/work/histogram/FromScratch/checkCodeSystNCentrMix" /*DiffNCentrMix*/};
  

     - Tracking Efficiency: 
       The macro is in systematic_sw. The files root used are:
       dirRefData = "~/work/histogram/FromScratch/checkCode" for the reference,   dirFiles = "~/work/histogram/FromScratch/checkCodeTrackEff" for Tracking Efficiency variation

      - ZYAM systematics:
        ZYAM is used for a comparison with UE residual so it has to be run BEFORE UEresidual
        The ZYAM systematic is estimated with systematic_sw. The reference files are contained in dirFileResults = "~/work/histogram/FromScratch/checkCode" passed to the function. The ZYAM results are in dirFileResults = "~/work/histogram/FromScratch/checkCodeZYAM; only the string ZYAM is added to access ZYAM zT functions.

     - UE Residual: 
       It has to be run always after ZYAM.
       It saves all the results in Systematic_checkCode. It does not distinguish between 0-30 and 0-10, 10-30.
       Already set the good directory in the macro use to combine all together all the systematics.

      

  PlotAllSystem.C
  It takes all the results from all the systFiles.root and combine all systematics together
