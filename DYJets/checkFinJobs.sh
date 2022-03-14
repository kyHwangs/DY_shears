echo Data_30
ll DMu_13TeV_Data_Syst_0*.root | wc
hadd DMu_13TeV_Data_DMu_Syst_0.root DMu_13TeV_Data_Syst_0_*.root
echo Data_SMu_30
ll DMu_13TeV_Data_SMu_Syst_0*.root | wc
hadd DMu_13TeV_Data_SMu_Syst_0.root DMu_13TeV_Data_SMu_Syst_0_*.root
echo DYJETS_34
ll DMu_13TeV_DYJets_UNFOLDING*.root | wc
hadd DMu_13TeV_DYJets_UNFOLDING_Syst_0.root DMu_13TeV_DYJets_UNFOLDING_Syst_0_*.root
echo TT_20
ll DMu_13TeV_TT_Syst_0*.root | wc
hadd DMu_13TeV_TT_Syst_0.root DMu_13TeV_TT_Syst_0*.root
echo Top_20
ll DMu_13TeV_Top_Syst_0*.root | wc
hadd DMu_13TeV_Top_Syst_0.root DMu_13TeV_Top_Syst_0*.root
echo ZZ_20
ll DMu_13TeV_ZZ_Syst_0_*.root | wc
hadd DMu_13TeV_ZZ_Syst_0.root DMu_13TeV_ZZ_Syst_0_*.root
echo WZ_20
ll  DMu_13TeV_WZ_Syst_0*.root | wc
hadd DMu_13TeV_WZ_Syst_0.root DMu_13TeV_WZ_Syst_0_*.root
echo WW_20
ll  DMu_13TeV_WWTo2L2Nu_Syst_0*.root | wc
hadd DMu_13TeV_WWTo2L2Nu_Syst_0.root DMu_13TeV_WWTo2L2Nu_Syst_0_*.root
echo W_20
ll DMu_13TeV_WToLNu_Syst*.root  | wc
hadd DMu_13TeV_WToLNu_Syst_0.root DMu_13TeV_WToLNu_Syst_0_*.root
echo VV_20
ll DMu_13TeV_VV_Syst*.root  | wc
hadd DMu_13TeV_VV_Syst_0.root DMu_13TeV_VV_Syst_0_*.root
echo merging DMu and SMu
hadd DMu_13TeV_Data_Syst_0.root DMu_13TeV_Data_DMu_Syst_0.root DMu_13TeV_Data_SMu_Syst_0.root



