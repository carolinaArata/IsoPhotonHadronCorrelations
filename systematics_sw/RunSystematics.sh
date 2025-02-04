######################################
######### Purity systematics #########
######################################
#
# PlotPuritySyst(pTMin, pTMax, "Mixed", ReferenceDirectoryWheretoAdd,  Set Mirror=On/Off, "0.40-1.00", Set 0_30 ON/OFF)
# Check the directory of the various purities in the code
root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/systematics_sw/PlotPuritySyst.C
PlotPuritySyst(18, 40, "Mixed", "~/work/histogram/FromScratch/checkCode", true, "0.40-1.00", true)
.q
EOF
#
#
#root -b -l <<EOF
#.L ~/work/histogram//IsoPhotonHadronCorrelations/systematics_sw/PlotUESyst.C
#PlotUESyst(18, 40, "ResultsNewMixMC_ZtMergedMore", "0.40-1.00",  "Mixed", true)
#.q
#EOF
#
######################################
######### ShSh systematics ###########
######################################
#
# ShShSyst(pTMin, pTMax, "Mixed",  mirrorOn, systematicDirectory, "0.40-1.00",  ReferenceDirectory, Set 0_30 ON)

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/systematics_sw/ShShSyst.C
ShShSyst(18, 40, "Mixed",  true,  "~/work/histogram/FromScratch/checkCodeSystShSh", "0.40-1.00", "~/work/histogram/FromScratch/checkCode", true)
.q
EOF
##
##

#####################################################
######### Tracking efficiency systematics ###########
#####################################################
#TrackIneffSyst(pTMin, pTMax,  mirrorOn/Off, "fPlot", "0.40-1.00", "directory/reference/data", "directory/to/systematics/data", 0_30 Centrality ON/OFF)

root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/systematics_sw/TrackIneffSyst.C
TrackIneffSyst(18, 40, true, "fPlot", "0.40-1.00", "~/work/histogram/FromScratch/checkCode", "~/work/histogram/FromScratch/checkCodeTrackEff", true)
.q
EOF
#
#
#####################################################
######### Number of centr bin systematics ###########
#####################################################
root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/systematics_sw/PlotNCentMix.C
PlotNCentMix(18, 40, true, "0.40-1.00", true)
.q
EOF
#
#
#####################################################
######### Number of centr bin systematics ###########
#####################################################
root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/systematics_sw/SystematicsZYAM.c
SystematicsZYAM(18, 40, true, "Mixed", "0.40-1.00", "~/work/histogram/FromScratch/checkCode", true)
.q
EOF
#
#
###########################################
######### Residual UE systematics #########
###########################################
root -b -l <<EOF
.L ~/work/histogram/IsoPhotonHadronCorrelations/systematics_sw/UEforClusterCheck.C
UEforClusterCheck(18, 40, "0.40-1.00",  "~/work/histogram/FromScratch/checkCode", true)
.q
EOF