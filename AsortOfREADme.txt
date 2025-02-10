
Before running analysis and plot, ctrl+f and check the directories of these files in IsoGammaHadron.C:
  - TFile *fPurity
  - TFile *fileData 
  - TFile *fileDataShStd 
  - TFile *fileDataShBkg 
  - TFile *fileDataMix 
  - TFile *fileMC[nShSh]

General macro: IsoGammaHadron.C; PlotIsoGammaHadron.C

1) Run Analysis and plots with: 

   RunAnalysis.sh

   This script contains the analysis and it also produces all the root files for the systematic studies 
    * Purity systematic is obtained from data, setting systPur = 0.9 or 1.1

    * ShSh Background systematic is obtained using different files for the ShSh background: 
      ShShSystMultBkg contains only Iso1 with multi Background: 0.40-2.00, 0.50-2.00, 0.60-2.00, 0.40-1.50 
      ShShSyst03510 contains Iso 1 and background 0.35-1.00
      REMEMBER TO SET nIso = 1, these trains were runned for Isolated ONLY.
      Change the shshBkg to the corresponding range and set systShSh = true

    * Tracking Efficiency, set systTrackIneff = true and this change the MC file to MC_GJShSh150/MC_GJTrackInEff.root;

    * number of centrality bins used, set systNMix = true
      NOTE: define two different bool for selecting number of centrality bins used
      Analysis run with old files NCentrBinMix18 and NCentrBinMix45: sigle root file containing all centralities: EMCAL_MB0_90.root; ShShBkg = 0.40-2.00; pThadr > 200 GeV;
      Update the code in the way that when it is reading fileMix the ShShBkg is 0.40-2.00; this can be set where the root files are declared before calling the Exec(...) : sShShNCentMix = "_ShSh0.40-2.00" for old files, sShShNCentMix = "_ShSh0.40-1.00" for new file.
      The shshbkg need to be defined by hand in the code

      The last train with N Cent Mix generate strange value in the systematics, not expected 

      The pT trig distribution for the old files were saved with ShShBkg: "_ShSh0.40-2.00" for mixed only, the name is updated in line 173-174 with 
      if (systNMix)
        hTriggerMix[iso][iSh]->SetName(sHistName + sIso + sShSh + sCent + "_hPtTriggerMixed");

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
