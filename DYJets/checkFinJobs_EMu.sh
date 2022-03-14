echo Data_SMu_30
ll EMu_13TeV_Data_SMu_Syst_0*.root | wc
hadd EMu_13TeV_Data_SMu_Syst_0.root EMu_13TeV_Data_SMu_Syst_0_*.root
echo DYJETS_34
ll EMu_13TeV_DYJets_UNFOLDING*.root | wc
hadd EMu_13TeV_DYJets_UNFOLDING_Syst_0.root EMu_13TeV_DYJets_UNFOLDING_Syst_0_*.root
echo TT_20
ll EMu_13TeV_TT_Syst_0*.root | wc
hadd EMu_13TeV_TT_Syst_0.root EMu_13TeV_TT_Syst_0*.root
echo Top_20
ll EMu_13TeV_Top_Syst_0*.root | wc
hadd EMu_13TeV_Top_Syst_0.root EMu_13TeV_Top_Syst_0*.root
echo ZZ_20
ll EMu_13TeV_ZZ_Syst_0_*.root | wc
hadd EMu_13TeV_ZZ_Syst_0.root EMu_13TeV_ZZ_Syst_0_*.root
echo WZ_20
ll  EMu_13TeV_WZ_Syst_0*.root | wc
hadd EMu_13TeV_WZ_Syst_0.root EMu_13TeV_WZ_Syst_0_*.root
echo WW_20
ll  EMu_13TeV_WWTo2L2Nu_Syst_0*.root | wc
hadd EMu_13TeV_WWTo2L2Nu_Syst_0.root EMu_13TeV_WWTo2L2Nu_Syst_0_*.root
echo W_20
ll EMu_13TeV_WToLNu_Syst*.root  | wc
hadd EMu_13TeV_WToLNu_Syst_0.root EMu_13TeV_WToLNu_Syst_0_*.root
echo VV_20
ll EMu_13TeV_VV_Syst*.root  | wc
hadd EMu_13TeV_VV_Syst_0.root EMu_13TeV_VV_Syst_0_*.root
mv EMu_13TeV_Data_SMu_Syst_0.root EMu_13TeV_Data_Syst_0.root



