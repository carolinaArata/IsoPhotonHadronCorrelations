###########################################################################################
############ Script for combining 0-10% and 10-30% centralities for analysis ##############
############ This procedure is repeated for systematics    ################################
###########################################################################################
#void Combine0_30(float ptMin, float ptMax, int iCen = 0, bool bMirror = Mirror ON/OFF, TString shshBkg = "0.40-1.00", TString dirFiles = "Directory/wih/the/files/from aanalysis", bool bPlot = true, TString dirPlot = "where/to/save/the/plot")

root -b -l <<EOF
.L Combine0_30.C
Combine0_30( 18, 40, 0, true, "0.40-1.00", "Output_checkCode", true, "Output_FigcheckCode")
.q
EOF

#----------------------------------------------------------------------------------------------
###############################################################################################
################# Combine 0-10% and 10-30% for systematics on purity ##########################
###############################################################################################
#----------------------------------------------------------------------------------------------

root -b -l <<EOF
.L Combine0_30.C
Combine0_30(18, 40, 0, true, "0.40-1.00", "Output_checkCodeSystPur09", true, "Output_FigcheckCodeSystPur09")
.q
EOF

##
##

root -b -l <<EOF
.L Combine0_30.C
Combine0_30(18, 40, 0, true, "0.40-1.00", "Output_checkCodeSystPur11", true, "Output_FigcheckCodeSystPur11")
.q
EOF

#
#----------------------------------------------------------------------------------------------
################################################################################################
########Combine 0-10% and 10-30% for systematics on N centrality bins used for mixed ###########
################################################################################################
#----------------------------------------------------------------------------------------------

root -b -l <<EOF
.L Combine0_30.C
Combine0_30( 18, 40, 0, true, "0.40-1.00", "Output_checkCodeSystNMix18", true, "Output_FigcheckCodeSystNCentrMix18")
.q
EOF

##
##

root -b -l <<EOF
.L Combine0_30.C
Combine0_30( 18, 40, 0, true, "0.40-1.00", "Output_checkCodeSystNMix45", true, "Output_FigcheckCodeSystNCentrMix45")
.q
EOF

#-----------------------------------------------------------------------------------------------
################################################################################################
########################### Run the analysis for systematics on ShSh ###########################
################################################################################################
#-----------------------------------------------------------------------------------------------

root -b -l <<EOF
.L Combine0_30.C
Combine0_30(18, 40, 0, true, "0.35-1.00", "Output_checkCodeSystShSh", true, "Output_FigcheckCodeSystShSh")
.q
EOF

##
##

root -b -l <<EOF
.L Combine0_30.C
Combine0_30(18, 40, 0, true, "0.40-0.80", "Output_checkCodeSystShSh", true, "Output_FigcheckCodeSystShSh")
.q
EOF

##
##

root -b -l <<EOF
.L Combine0_30.C
Combine0_30(18, 40, 0, true, "0.50-1.00", "Output_checkCodeSystShSh", true, "Output_FigcheckCodeSystShSh")
.q
EOF

#-----------------------------------------------------------------------------------------------
################################################################################################
################## Run the analysis for systematics on Tracking Inefficiency ###################
################################################################################################
#-----------------------------------------------------------------------------------------------
root -b -l <<EOF
.L Combine0_30.C
Combine0_30(18, 40, 0, true, "0.40-1.00", "Output_checkCodeTrackEff", true, "Output_FigcheckCodeTrackEff")
.q
EOF

#----------------------------------------------------------------------------------------------
################################################################################################
############################## Run the analysis for ZYAM #######################################
################################################################################################
#----------------------------------------------------------------------------------------------
root -b -l <<EOF
.L Combine0_30.C
Combine0_30(18, 40, 0, true, "0.40-1.00", "Output_checkCodeZYAM", true, "Output_FigcheckCodeZYAM")
.q
EOF






