######################################################################
############ Run the analysis and the plotting          ##############
############ This procedure is repeated for systematics ##############
######################################################################

#IsoGammaHadron(float ptTrMin, float ptTrMax, TString sFileDirShSig = directory where the trains for the analysis are, bMirror = Turn Mirror ON/OFF, TString shshBkg = "0.40-1.00", TString dirFiles = directory to save output files, double systPur = for Systematic on Purity: 1, 0.9, 1.1, bool bZYAM = Turn ZYAM analysis ON/OFF , bool bPlot = true, double phiMin, double phiMax, bool systShSh = ON/ OFF, bool systTrackIneff = ON / OFF)


#PlotIsoGammaHadron(float ptMin , float ptMax , TString dirPlot = "Where to save output plot", TString shshBkg = "0.40-1.00", TString dirFiles = "Analysis root files obtained from IsoGammaHadron.C")

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
IsoGammaHadron(18, 40, "~/work/histogram/DataSh100_AssocPt500", true, "0.40-1.00", "~/work/histogram/FromScratch/checkCode", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), false, false)
.q
EOF
##

root -b -l << EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/PlotIsoGammaHadron.C
PlotIsoGammaHadron(18, 40,  "~/work/histogram/FromScratch/FigcheckCode", "0.40-1.00", "~/work/histogram/FromScratch/checkCode")
.q
EOF

#----------------------------------------------------------------------------------------------
###############################################################################################
########################## Run the analysis for systematics on purity #########################
###############################################################################################
#----------------------------------------------------------------------------------------------

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
IsoGammaHadron(18, 40, "~/work/histogram/DataSh100_AssocPt500", true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeSystPur09", 0.9, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), false ,false)
.q
EOF
#

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
IsoGammaHadron(18, 40, "~/work/histogram/DataSh100_AssocPt500", true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeSystPur11", 1.1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), false ,false)
.q
EOF

#----------------------------------------------------------------------------------------------
################################################################################################
############ Run the analysis for systematics on N centrality bins used for mixed ##############
################################################################################################
#----------------------------------------------------------------------------------------------
#root -b -l <<EOF
#.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
#IsoGammaHadron(18, 40, "OLDDatasetConfigPtAssoc200_ptIso2GeV_diffShSh/AllTogether", true, "0.40-2.00", "SystematicsNMix45_ZtMergedMore", 1, false, TMath::Pi() * 3 / 5, TMath::Pi())
#.q
#EOF

#root -b -l <<EOF
#.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
#IsoGammaHadron(18, 40, "OLDDatasetConfigPtAssoc200_ptIso2GeV_diffShSh/AllTogether", true, "0.40-2.00", "SystematicsNMix9_ZtMergedMore", 1, false, TMath::Pi() * 3 / 5, TMath::Pi())
#.q
#EOF

#
#root -b -l <<EOF
#.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
#IsoGammaHadron(16, 40, "~/work/histogram/DataSh100_AssocPt500", true, "0.40-1.00", "SystematicsNMix45_ZtMerged", 1, false, TMath::Pi() * 3 / 5, TMath::Pi())
#.q
#EOF

#-----------------------------------------------------------------------------------------------
################################################################################################
########################### Run the analysis for systematics on ShSh ###########################
################################################################################################
#-----------------------------------------------------------------------------------------------

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
IsoGammaHadron(18, 40, "~/work/histogram/DataSh100_AssocPt500", true, "0.35-1.00", "~/work/histogram/FromScratch/checkCodeSystShSh", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), true, false)
.q
EOF
##
##
root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
IsoGammaHadron(18, 40, "~/work/histogram/DataSh100_AssocPt500", true, "0.40-1.50", "~/work/histogram/FromScratch/checkCode/checkCodeSystShSh", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), true, false)
.q
EOF
##
##
root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
IsoGammaHadron(18, 40, "~/work/histogram/DataSh100_AssocPt500", true, "0.40-2.00", "~/work/histogram/FromScratch/checkCode/checkCodeSystShSh", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), true, false)
.q
EOF

#-----------------------------------------------------------------------------------------------
################################################################################################
################## Run the analysis for systematics on Tracking Inefficiency ###################
################################################################################################
#-----------------------------------------------------------------------------------------------

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
IsoGammaHadron(18, 40, "~/work/histogram/DataSh100_AssocPt500", true, "0.40-1.00", "~/work/histogram/FromScratch/checkCode/checkCodeTrackEff", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), false, true)
.q
EOF

#----------------------------------------------------------------------------------------------
################################################################################################
############################## Run the analysis for ZYAM #######################################
################################################################################################
#----------------------------------------------------------------------------------------------

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
IsoGammaHadron(18, 40, "~/work/histogram/DataSh100_AssocPt500", true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeZYAM", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), false, false)
.q
EOF

##-------------------------------
#
#######################




