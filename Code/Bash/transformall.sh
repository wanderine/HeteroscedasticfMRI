#!/bin/bash

dataset=MixedGambles
analysis=Hetero
ten=10

for subject in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
do

	echo "Transforming subject $subject"

	for regressor in 1 2 3 4
	do

		if [ "$regressor" -lt "$ten" ]; then

			# Apply combined transformation
			flirt -in Results/${dataset}${analysis}/bold_${subject}_mc_sm_beta_fullposterior_regressor000${regressor}.nii.gz -ref MNI152_T1_2mm_brain.nii.gz -applyxfm -init Results/${dataset}${analysis}/subject${subject}_fMRI_to_MNI.mat -out Results/${dataset}${analysis}/bold_${subject}_mc_sm_beta_fullposterior_regressor000${regressor}_MNI.nii.gz

			flirt -in Results/${dataset}${analysis}/bold_${subject}_mc_sm_beta_regressor000${regressor}.nii.gz -ref MNI152_T1_2mm_brain.nii.gz -applyxfm -init Results/${dataset}${analysis}/subject${subject}_fMRI_to_MNI.mat -out Results/${dataset}${analysis}/bold_${subject}_mc_sm_beta_regressor000${regressor}_MNI.nii.gz
			
		else

			echo "5"

		fi
	
	done

done


