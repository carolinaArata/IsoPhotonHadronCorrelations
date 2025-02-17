###########################################################################################
############ Script for combining 0-10% and 10-30% centralities for analysis ##############
############ This procedure is repeated for systematics    ################################
###########################################################################################
#void Combine0_30(float ptMin, float ptMax, int iCen = 0, bool bMirror = Mirror ON/OFF, TString shshBkg = "0.40-1.00", TString dirFiles = "Directory/wih/the/files/from aanalysis", bool bPlot = true, TString dirPlot = "where/to/save/the/plot")

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/Combine0_30.C
Combine0_30( 18, 40, 0, true, "0.40-1.00", "~/work/histogram/FromScratch/checkCode", true, "~/work/histogram/FromScratch/FigcheckCode")
.q
EOF

#----------------------------------------------------------------------------------------------
###############################################################################################
################# Combine 0-10% and 10-30% for systematics on purity ##########################
###############################################################################################
#----------------------------------------------------------------------------------------------

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/Combine0_30.C
Combine0_30(18, 40, 0, true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeSystPur09", true, "~/work/histogram/FromScratch/FigcheckCodeSystPur09")
.q
EOF

##
##

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/Combine0_30.C
Combine0_30(18, 40, 0, true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeSystPur11", true, "~/work/histogram/FromScratch/FigcheckCodeSystPur11")
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
Combine0_30( 18, 40, 0, true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeSystNMix18", true, "~/work/histogram/FromScratch/FigcheckCodeSystNCentrMix18")
.q
EOF

##
##

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/Combine0_30.C
Combine0_30( 18, 40, 0, true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeSystNMix45", true, "~/work/histogram/FromScratch/FigcheckCodeSystNCentrMix45")
.q
EOF

#-----------------------------------------------------------------------------------------------
################################################################################################
########################### Run the analysis for systematics on ShSh ###########################
################################################################################################
#-----------------------------------------------------------------------------------------------

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/Combine0_30.C
Combine0_30(18, 40, 0, true, "0.35-1.00", "~/work/histogram/FromScratch/checkCodeSystShSh", true, "~/work/histogram/FromScratch/FigcheckCodeSystShSh")
.q
EOF

##
##

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/Combine0_30.C
Combine0_30(18, 40, 0, true, "0.40-1.50", "~/work/histogram/FromScratch/checkCodeSystShSh", true, "~/work/histogram/FromScratch/FigcheckCodeSystShSh")
.q
EOF

##
##

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/Combine0_30.C
Combine0_30(18, 40, 0, true, "0.40-2.00", "~/work/histogram/FromScratch/checkCodeSystShSh", true, "~/work/histogram/FromScratch/FigcheckCodeSystShSh")
.q
EOF

#-----------------------------------------------------------------------------------------------
################################################################################################
################## Run the analysis for systematics on Tracking Inefficiency ###################
################################################################################################
#-----------------------------------------------------------------------------------------------
root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/Combine0_30.C
Combine0_30(18, 40, 0, true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeTrackEff", true, "~/work/histogram/FromScratch/FigcheckCodeTrackEff")
.q
EOF

#----------------------------------------------------------------------------------------------
################################################################################################
############################## Run the analysis for ZYAM #######################################
################################################################################################
#----------------------------------------------------------------------------------------------
root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/Combine0_30.C
Combine0_30(18, 40, 0, true, "0.40-1.00", "~/work/histogram/FromScratch/checkCodeZYAM", true, "~/work/histogram/FromScratch/FigcheckCodeZYAM")
.q
EOF






