#!/bin/bash

dataset=MixedGambles
datadir=/home/andek/Research_projects/HeteroGauss/Data/OpenfMRI/$dataset
ten=10

for subject in {1..13}
do

	echo "Registering subject $subject"

	if [ "$subject" -lt "$ten" ]; then

		cp $datadir/sub00$subject/BOLD/task001_run001/bold.nii.gz .

		cp $datadir/sub00$subject/anatomy/highres001_brain.nii.gz T1_brain_${subject}.nii.gz
		cp $datadir/sub00$subject/anatomy/highres001.nii.gz T1_${subject}.nii.gz

	else

		cp $datadir/sub0$subject/BOLD/task001_run001/bold.nii.gz .

		cp $datadir/sub0$subject/anatomy/highres001_brain.nii.gz T1_brain_${subject}.nii.gz
		cp $datadir/sub0$subject/anatomy/highres001.nii.gz T1_${subject}.nii.gz

	fi

	# T1 to MNI
	flirt -in T1_brain_${subject}.nii.gz -ref MNI152_T1_2mm_brain.nii.gz -omat subject${subject}_T1_to_MNI.mat -out subject${subject}_T1_MNI.nii.gz
	
	echo -e "\nSegmenting T1 volume\n"
	fast T1_brain_${subject}.nii.gz
	mv T1_brain_${subject}_pve_1.nii.gz T1_brain_${subject}_gmseg.nii.gz
	mv T1_brain_${subject}_pve_2.nii.gz T1_brain_${subject}_wmseg.nii.gz

	echo -e "\nExtracting first fMRI volume\n"	
	3dTcat 'bold.nii.gz[0]' -prefix bold_${subject}_first.nii.gz  

	echo -e "\nRegistering T1 volume to fMRI\n"
	epi_reg --epi=bold_${subject}_first.nii.gz --t1=T1_${subject}.nii.gz --t1brain=T1_brain_${subject}.nii.gz --out=fMRI_to_T1_${subject}

	# Combine transformations
	convert_xfm -omat subject${subject}_fMRI_to_MNI.mat -concat  subject${subject}_T1_to_MNI.mat fMRI_to_T1_${subject}.mat

	cp subject${subject}_fMRI_to_MNI.mat Results/${dataset}Hetero

	cp subject${subject}_fMRI_to_MNI.mat Results/${dataset}Homo

done


