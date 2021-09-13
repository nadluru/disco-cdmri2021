# initialization
export dataroot=DiSCo_datasets
export coderoot=disco-cdmri2021
export trainroot=$dataroot/Training_DiSCo1
export validationroot=$dataroot/Validation_DiSCo3
export testroot=$dataroot/Test_DiSCo2

export dwi=($trainroot/DiSCo1_DWI_RicianNoise-snr30.nii.gz $validationroot/DiSCo3_DWI_RicianNoise-snr30.nii.gz $testroot/DiSCo2_DWI_RicianNoise-snr30.nii.gz)
export rois=($trainroot/DiSCo1_ROIs.nii.gz $validationroot/DiSCo3_ROIs.nii.gz $testroot/DiSCo2_ROIs.nii.gz)
export scheme=$coderoot/gradient_files/DiSCo_gradients.scheme
export bvals=$coderoot/gradient_files/DiSCo_gradients.bvals
export bvecs=$coderoot/gradient_files/DiSCo_gradients_fsl.bvecs
export connectomes=($trainroot/DiSCo1_Connectivity_Matrix_Cross-Sectional_Area.txt $validationroot/DiSCo3_Connectivity_Matrix_Cross-Sectional_Area.txt)

# mrtrix3 transformations
parallel --plus '
# signal corrections
mrconvert {1} - -fslgrad $bvecs $bvals | dwidenoise - {1..}_denoised.mif -noise {1..}_noise.mif -force
mrdegibbs {1..}_denoised.mif - | mrcalc - 2 -pow $(mrcalc {1..}_noise.mif -finite {1..}_noise.mif 0 -if -) 2 -pow -sub -abs -sqrt - | mrcalc - -finite - 0 -if {1..}_rc.mif -force

# fod estimation
dwi2response {3} -lmax 0,10,10,10,10 -sfwm 1 -gm 0.5 -csf 5 {1..}_rc.mif {1..}_sfwm_{3}.txt {1..}_gm_{3}.txt {1..}_csf_{3}.txt
dwi2fod -lmax 10,0 msmt_csd {1..}_rc.mif {1..}_sfwm_{3}.txt {1..}_wmfod_{3}.nii.gz {1..}_csf_{3}.txt {1..}_csffod_{3}.mif

# masking
mrconvert {1..}_wmfod_{3}.mif -coord 3 0 -axes 0:2 - -nthreads $nthreads | mrcalc - {4} -ge 1 0 -if {2} -subtract 0 -ge 1 0 -if {1..}_wmfod_dc_mask_sub_roi_{3}.nii.gz
connectedcomp {1..}_wmfod_dc_mask_sub_roi_{3} {1..}_wmfod_dc_mask_sub_roi_conncomp_{3}
mrcalc {1..}_wmfod_dc_mask_sub_roi_conncomp_{3}.nii.gz $(LabelGeometryMeasures 3 {1..}_wmfod_dc_mask_sub_roi_conncomp_{3}.nii.gz | awk -F " " '\''{print $1,$2}'\'' | sed 1d | sort -k2 -g -r | head -n1 | awk '\''{print $1}'\'') -eq 1 0 -if {1..}_wmfod_dc_mask_sub_roi_conncomp_{3}_filtered.nii.gz' ::: ${dwi[@]} :::+ ${rois[@]} ::: dhollander :::+ 0.02

# trekker + tcksift2 + tck2connectome
parallel '
# trekker
trekker -fod {1}_{3}.nii.gz -seed_image {1}_dc_mask_sub_roi_conncomp_{3}_filtered.nii.gz -seed_count 1000000 -pathway_A=require_entry {2} -pathway_B=require_entry {2} -enableOutputOverwrite -output {1}_{3}_{4}_{5} -writeColors -numberOfThreads $nthreads -writeInterval 100 -minLength {5} -verboseLevel 0

# mrtrix3
tckconvert {1}_{3}_{4}_{5}.vtk {1}_{3}_{4}_{5}.tck
tcksift2 {1}_{3}_{4}_{5}.tck {1}_{3}.nii.gz {1}_{3}_{4}_{5}_weights.txt
tck2connectome {1}_{3}_{4}_{5}.tck {2} {1}_{3}_{4}_{5}_connectome.csv -tck_weights_in {1}_{3}_{4}_{5}_weights.txt -symmetric
sed -i "s:,: :g" {1}_{3}_{4}_{5}_connectome.csv' ::: ${dwi[@]/.nii.gz/_wmfod} :::+ ${rois[@]} ::: dhollander ::: trekker ::: 0

# scoring (UWMadison_Submission_01.txt)
parallel --plus -j1 -k 'echo {1/},{2/..},{3},{4},{5},$(python $coderoot/compute_correlation.py {1} {2/..}_wmfod_{3}_{4}_{5}_connectome.csv | sed "s:.* ::")' ::: ${connectomes[@]} :::+ ${dwi[@]:0:2} ::: dhollander ::: trekker ::: 0

#--------------------------------------------------------
# 10M version with tcksift2 x tck2connectome combinations
#--------------------------------------------------------
parallel '
trekker -fod {1}_{3}.nii.gz -seed_image {1}_dc_mask_sub_roi_conncomp_{3}_filtered.nii.gz -seed_count 1000000 -pathway_A=require_entry {2} -pathway_B=require_entry {2} -enableOutputOverwrite -output {1}_{3}_{4}_{5}_run_{6} -minLength {5}
tckconvert {1}_{3}_{4}_{5}_run_{6}.vtk {1}_{3}_{4}_{5}_run_{6}.tck' ::: ${dwi[@]/.nii.gz/_wmfod} :::+ ${rois[@]} ::: dhollander ::: trekker ::: 0 ::: {1..10}
parallel tckedit {1}_{3}_{4}_{5}_run_*.tck {1}_{3}_{4}_{5}_10M.tck ::: ${dwi[@]/.nii.gz/_wmfod} :::+ ${rois[@]} ::: dhollander ::: trekker ::: 0

parallel '
tcksift2 {1}_{3}_{4}_{5}_10M.tck {1}_{3}.mif {1}_{3}_{4}_{5}_10M_weights.txt
tck2connectome {1}_{3}_{4}_{5}_10M.tck {2} {1}_{3}_{4}_{5}_10M_connectome.csv -tck_weights_in {1}_{3}_{4}_{5}_10M_weights.txt -symmetric
sed -i "s:,: :g" {1}_{3}_{4}_{5}_10M_connectome.csv' ::: ${dwi[@]/.nii.gz/_wmfod} :::+ ${rois[@]} ::: dhollander ::: trekker ::: 0

#--------------------------------------------------------
# scoring 10M version (UWMadison_Submission_02.txt)
#--------------------------------------------------------
parallel --plus -j1 -k '
echo {1/},{2/..},{3},{4},{5},10M,$(python $coderoot/compute_correlation.py {1} {2..}_wmfod_{3}_{4}_{5}_10M_connectome.csv | sed "s:.* ::")' ::: ${connectomes[@]} :::+ ${dwi[@]:0:2} ::: dhollander ::: trekker ::: 0

# 81 combinations (tcksift2 (3 x 3 x 3 = 27) x tck2connectome (3)) with 10M tracts
# tcksift2
parallel --plus '
tcksift2 {1}_{3}_{4}_{5}_10M.tck {1}_{3}.nii.gz {1}_{3}_{4}_{5}_10M_weights_{6}_{7}{8/default/}.txt -min_td_frac {6} -min_iters {7} $(echo {8/default/})' ::: ${dwi[@]/.nii.gz/_wmfod} :::+ ${rois[@]} ::: dhollander ::: trekker ::: 0 ::: 0.05 0.1 0.2 ::: 10 20 30 ::: default -no_dilate_lut "-no_dilate_lut -make_null_lobes"

# tck2connectome
parallel --plus '
tck2connectome {1}_{3}_{4}_{5}_10M.tck {2} {1}_{3}_{4}_{5}_10M_connectome_{6}_{7}{8/default/}{9/default/}.csv -tck_weights_in {1}_{3}_{4}_{5}_10M_weights_{6}_{7}{8/default/}.txt -symmetric $(echo {9/default/})
sed -i "s:,: :g" {1}_{3}_{4}_{5}_10M_connectome_{6}_{7}{8/default/}{9/default/}.csv' ::: ${dwi[@]/.nii.gz/_wmfod} :::+ ${rois[@]} ::: dhollander ::: trekker ::: 0 ::: 0.05 0.1 0.2 ::: 10 20 30 ::: default -no_dilate_lut "-no_dilate_lut -make_null_lobes" ::: default -assignment_end_voxels "-assignment_radial_search 2"

#--------------------------------------------------------
# scoring the 81 combinations (UWMadison_Submission_03-83.txt)
#--------------------------------------------------------
parallel --plus -j1 -k '
echo {1/},{2/..},{3},{4},{5},10M,{6},{7},{8},{9},$(python $coderoot/compute_correlation.py {1} {2..}_wmfod_{3}_{4}_{5}_10M_connectome_{6}_{7}{8/default/}{9/default/}.csv | sed "s:.* ::")' ::: ${connectomes[@]} :::+ ${dwi[@]:0:2} ::: dhollander ::: trekker ::: 0 ::: 0.05 0.1 0.2 ::: 10 20 30 ::: default -no_dilate_lut "-no_dilate_lut -make_null_lobes" ::: default -assignment_end_voxels "-assignment_radial_search 2"