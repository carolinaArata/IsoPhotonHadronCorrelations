###########################################################################################
############ Script for combining 0-10% and 10-30% centralities for analysis ##############
############ This procedure is repeated for systematics    ################################
###########################################################################################
#void Combine0_30(float ptMin, float ptMax, int iCen = 0, bool bMirror = Mirror ON/OFF, TString shshBkg = "0.40-1.00", TString dirFiles = "Directory/wih/the/files/from aanalysis", double systPur = 1, 0.9, 1.1, bool bZYAM = ZYAM ON/OFF, bool bPlot = true, TString dirPlot = "where/to/save/the/plot")

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/Combine0_30.C
Combine0_30( 18, 40, 0, true, "0.40-1.00", "~/work/histogram/FromScratch/checkCode", 1, false, true, "~/work/histogram/FromScratch/FigcheckCode")
.q
EOF
#
#----------------------------------------------------------------------------------------------
###############################################################################################
################# Combine 0-10% and 10-30% for systematics on purity ##########################
###############################################################################################
#----------------------------------------------------------------------------------------------
root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/Combine0_30.C
Combine0_30(18, 40, 0, true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeSystPur09", 0.9, false, true, "~/work/histogram/FromScratch/FigcheckCodeSystPur09")
.q
EOF
###
root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/Combine0_30.C
Combine0_30(18, 40, 0, true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeSystPur11", 1.1, false, true, "~/work/histogram/FromScratch/FigcheckCodeSystPur11")
.q
EOF
#
#----------------------------------------------------------------------------------------------
################################################################################################
########Combine 0-10% and 10-30% for systematics on N centrality bins used for mixed ###########
################################################################################################
#----------------------------------------------------------------------------------------------
root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/Combine0_30.C
Combine0_30( 18, 40, 0, true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeSystNCentrMix", 1, false, true, "~/work/histogram/FromScratch/FigcheckCodeSystNCentrMix")
.q
EOF
#-----------------------------------------------------------------------------------------------
################################################################################################
########################### Run the analysis for systematics on ShSh ###########################
################################################################################################
#-----------------------------------------------------------------------------------------------

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/Combine0_30.C
Combine0_30(18, 40, 0, true, "0.35-1.00", "~/work/histogram/FromScratch/checkCodeSystShSh", 1, false, true, "~/work/histogram/FromScratch/FigcheckCodeSystShSh")
.q
EOF

##
##

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/Combine0_30.C
Combine0_30(18, 40, 0, true, "0.40-1.50", "~/work/histogram/FromScratch/checkCodeSystShSh", 1, false, true, "~/work/histogram/FromScratch/FigcheckCodeSystShSh")
.q
EOF

##
##

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/Combine0_30.C
Combine0_30(18, 40, 0, true, "0.40-2.00", "~/work/histogram/FromScratch/checkCodeSystShSh", 1, false, true, "~/work/histogram/FromScratch/FigcheckCodeSystShSh")
.q
EOF

#-----------------------------------------------------------------------------------------------
################################################################################################
################## Run the analysis for systematics on Tracking Inefficiency ###################
################################################################################################
#-----------------------------------------------------------------------------------------------
root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/Combine0_30.C
Combine0_30(18, 40, 0, true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeTrackEff", 1, false, true, "~/work/histogram/FromScratch/FigcheckCodeTrackEff")
.q
EOF
#

#----------------------------------------------------------------------------------------------
################################################################################################
############################## Run the analysis for ZYAM #######################################
################################################################################################
#----------------------------------------------------------------------------------------------
root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/Combine0_30.C
Combine0_30(18, 40, 0, true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeZYAM", 1, true, true, "~/work/histogram/FromScratch/FigcheckCodeZYAM")
.q
EOF

##-------------------------------
#
#######################




