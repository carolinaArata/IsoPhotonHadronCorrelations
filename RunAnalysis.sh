######################################################################
############ Run the analysis and the plotting          ##############
############ This procedure is repeated for systematics ##############
######################################################################

#IsoGammaHadron(ptTrMin, ptTrMax, sFileDirShSig = directory/where/the/trains/for/analysis/are, bMirror = Turn Mirror ON/OFF, TString shshBkg = "0.40-1.00", TString dirFiles = directory/to/save/output/files, double systPur = for Systematic on Purity: 1, 0.9, 1.1, bool bZYAM = Turn ZYAM analysis ON/OFF , bool bPlot = plotting for checking ON/OFF, phiMin, phiMax, bool systShSh = ShSh systematic ON/OFF, bool systTrackIneff = track effic systematic ON/OFF, bool systNMix18 = systematic for N Centr Bin 18 for Mix ON/OFF, bool systNMix45 = systematic for N Centr Bin 45 for Mix ON/OFF)

#PlotIsoGammaHadron(float ptMin , float ptMax , TString dirPlot = "Where to save output plot", TString shshBkg = "0.40-1.00", TString dirFiles = "Analysis root files obtained from IsoGammaHadron.C")

root -b -l <<EOF
.L IsoGammaHadron.C
IsoGammaHadron(18, 40, "RootFiles/DataSh100_AssocPt500", true, "0.40-1.00", "Output_checkCode", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), false, false, false, false)
.q
EOF
#

root -b -l <<EOF
.L PlotIsoGammaHadron.C
PlotIsoGammaHadron(18, 40,  "Output_FigcheckCode", "0.40-1.00", "Output_checkCode")
.q
EOF

#----------------------------------------------------------------------------------------------
###############################################################################################
########################## Run the analysis for systematics on purity #########################
###############################################################################################
#----------------------------------------------------------------------------------------------

root -b -l <<EOF
.L IsoGammaHadron.C
IsoGammaHadron(18, 40, "RootFiles/DataSh100_AssocPt500", true, "0.40-1.00", "Output_checkCodeSystPur09", 0.9, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), false ,false, false, false)
.q
EOF
#

root -b -l <<EOF
.L IsoGammaHadron.C
IsoGammaHadron(18, 40, "RootFiles/DataSh100_AssocPt500", true, "0.40-1.00", "Output_checkCodeSystPur11", 1.1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), false ,false, false, false)
.q
EOF

#----------------------------------------------------------------------------------------------
################################################################################################
############ Run the analysis for systematics on N centrality bins used for mixed ##############
################################################################################################
#----------------------------------------------------------------------------------------------
root -b -l <<EOF
.L IsoGammaHadron.C
IsoGammaHadron(18, 40, "RootFiles/DataSh100_AssocPt500", true, "0.40-1.00", "Output_checkCodeSystNMix18", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), false, false, true, false)
.q
EOF

root -b -l <<EOF
.L PlotIsoGammaHadron.C
PlotIsoGammaHadron(18, 40,  "Output_FigcheckCodeSystNCentrMix18", "0.40-1.00", "Output_checkCodeSystNMix18")
.q
EOF

root -b -l <<EOF
.L IsoGammaHadron.C
IsoGammaHadron(18, 40, "RootFiles/DataSh100_AssocPt500", true, "0.40-1.00", "Output_checkCodeSystNMix45", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), false, false, false, true)
.q
EOF

root -b -l <<EOF
.L PlotIsoGammaHadron.C
PlotIsoGammaHadron(18, 40,  "Output_FigcheckCodeSystNCentrMix45", "0.40-1.00", "Output_checkCodeSystNMix45")
.q
EOF

#-----------------------------------------------------------------------------------------------
################################################################################################
########################### Run the analysis for systematics on ShSh ###########################
################################################################################################
#-----------------------------------------------------------------------------------------------

root -b -l <<EOF
.L IsoGammaHadron.C
IsoGammaHadron(18, 40, "RootFiles/DataSh100_AssocPt500", true, "0.35-1.00", "Output_checkCodeSystShSh", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), true, false, false, false)
.q
EOF
##
##
root -b -l <<EOF
.L IsoGammaHadron.C
IsoGammaHadron(18, 40, "RootFiles/DataSh100_AssocPt500", true, "0.40-0.80", "Output_checkCodeSystShSh", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), true, false, false, false)
.q
EOF
##
##
root -b -l <<EOF
.L IsoGammaHadron.C
IsoGammaHadron(18, 40, "RootFiles/DataSh100_AssocPt500", true, "0.50-1.00", "Output_checkCodeSystShSh", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), true, false, false, false)
.q
EOF

#-----------------------------------------------------------------------------------------------
################################################################################################
################## Run the analysis for systematics on Tracking Inefficiency ###################
################################################################################################
#-----------------------------------------------------------------------------------------------

root -b -l <<EOF
.L IsoGammaHadron.C
IsoGammaHadron(18, 40, "RootFiles/DataSh100_AssocPt500", true, "0.40-1.00", "Output_checkCodeTrackEff", 1, false, true, TMath::Pi() * 3 / 5, TMath::Pi(), false, true, false, false)
.q
EOF

#----------------------------------------------------------------------------------------------
################################################################################################
############################## Run the analysis for ZYAM #######################################
################################################################################################
#----------------------------------------------------------------------------------------------

root -b -l <<EOF
.L IsoGammaHadron.C
IsoGammaHadron(18, 40, "RootFiles/DataSh100_AssocPt500", true, "0.40-1.00", "Output_checkCodeZYAM", 1, true, true, TMath::Pi() * 3 / 5, TMath::Pi(), false, false, false, false)
.q
EOF

##-------------------------------
#
#######################
