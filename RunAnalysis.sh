######################################################################
############ Run the analysis and the plotting          ##############
############ This procedure is repeated for systematics ##############
######################################################################

#IsoGammaHadron(ptTrMin, ptTrMax, sFileDirShSig = directory/where/the/trains/for/analysis/are, bMirror = Turn Mirror ON/OFF, TString shshBkg = "0.40-1.00", TString dirFiles = directory/to/save/output/files, double systPur = for Systematic on Purity: 1, 0.9, 1.1, bool bZYAM = Turn ZYAM analysis ON/OFF , bool bPlot = plotting for checking ON/OFF, phiMin, phiMax, bool systShSh = ShSh systematic ON/OFF, bool systTrackIneff = track effic systematic ON/OFF, bool systNMix18 = systematic for N Centr Bin 18 for Mix ON/OFF, bool systNMix45 = systematic for N Centr Bin 45 for Mix ON/OFF)

#PlotIsoGammaHadron(float ptMin , float ptMax , TString dirPlot = "Where to save output plot", TString shshBkg = "0.40-1.00", TString dirFiles = "Analysis root files obtained from IsoGammaHadron.C")

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
IsoGammaHadron(18, 40, "~/work/histogram/DataSh100_AssocPt500", true, "0.40-1.00", "~/work/histogram/FromScratch/checkCode", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), false, false, false, false)
.q
EOF
#

root -b -l <<EOF
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
IsoGammaHadron(18, 40, "~/work/histogram/DataSh100_AssocPt500", true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeSystPur09", 0.9, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), false ,false, false, false)
.q
EOF
#

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
IsoGammaHadron(18, 40, "~/work/histogram/DataSh100_AssocPt500", true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeSystPur11", 1.1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), false ,false, false, false)
.q
EOF

#----------------------------------------------------------------------------------------------
################################################################################################
############ Run the analysis for systematics on N centrality bins used for mixed ##############
################################################################################################
#----------------------------------------------------------------------------------------------
root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
IsoGammaHadron(18, 40, "~/work/histogram/DataSh100_AssocPt500", true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeSystNMix18", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), false, false, true, false)
.q
EOF

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
IsoGammaHadron(18, 40, "~/work/histogram/DataSh100_AssocPt500", true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeSystNMix45", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), false, false, false, true)
.q
EOF

#-----------------------------------------------------------------------------------------------
################################################################################################
########################### Run the analysis for systematics on ShSh ###########################
################################################################################################
#-----------------------------------------------------------------------------------------------

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
IsoGammaHadron(18, 40, "~/work/histogram/DataSh100_AssocPt500", true, "0.35-1.00", "~/work/histogram/FromScratch/checkCodeSystShSh", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), true, false, false, false)
.q
EOF
##
##
root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
IsoGammaHadron(18, 40, "~/work/histogram/DataSh100_AssocPt500", true, "0.40-1.50", "~/work/histogram/FromScratch/checkCodeSystShSh", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), true, false, false, false)
.q
EOF
##
##
root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
IsoGammaHadron(18, 40, "~/work/histogram/DataSh100_AssocPt500", true, "0.40-2.00", "~/work/histogram/FromScratch/checkCodeSystShSh", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), true, false, false, false)
.q
EOF

#-----------------------------------------------------------------------------------------------
################################################################################################
################## Run the analysis for systematics on Tracking Inefficiency ###################
################################################################################################
#-----------------------------------------------------------------------------------------------

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
IsoGammaHadron(18, 40, "~/work/histogram/DataSh100_AssocPt500", true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeTrackEff", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), false, true, false, false)
.q
EOF

#----------------------------------------------------------------------------------------------
################################################################################################
############################## Run the analysis for ZYAM #######################################
################################################################################################
#----------------------------------------------------------------------------------------------

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/IsoGammaHadron.C
IsoGammaHadron(18, 40, "~/work/histogram/DataSh100_AssocPt500", true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeZYAM", 1, true, true, TMath::Pi() * 3 / 5, TMath::Pi(), false, false, false, false)
.q
EOF

##-------------------------------
#
#######################
