#!/bin/bash

dataset=RhymeJudgment
datadir=/home/andek/Research_projects/HeteroGauss/Data/OpenfMRI/$dataset
ten=10
eleven=11

for subject in {1..13}
do

	if [ "$subject" -lt "$ten" ]; then

		cp $datadir/sub00$subject/BOLD/task001_run001/bold.nii.gz bold_${subject}.nii.gz

		cp $datadir/sub00$subject/anatomy/highres001_brain.nii.gz T1_brain_${subject}.nii.gz
		cp $datadir/sub00$subject/anatomy/highres001.nii.gz T1_${subject}.nii.gz
		
		cp $datadir/sub00$subject/model/model001/onsets/task001_run001/cond001.txt .
		cp $datadir/sub00$subject/model/model001/onsets/task001_run001/cond002.txt .

	else

		cp $datadir/sub0$subject/BOLD/task001_run001/bold.nii.gz bold_${subject}.nii.gz

		cp $datadir/sub0$subject/anatomy/highres001_brain.nii.gz T1_brain_${subject}.nii.gz
		cp $datadir/sub0$subject/anatomy/highres001.nii.gz T1_${subject}.nii.gz
		
		cp $datadir/sub0$subject/model/model001/onsets/task001_run001/cond001.txt .
		cp $datadir/sub0$subject/model/model001/onsets/task001_run001/cond002.txt .

	fi

	#------------------

	events1=`cat cond001.txt | wc -l`
	events2=`cat cond002.txt | wc -l`

	# Add NumEvents as first line to task files
	sed -i "1s/^/NumEvents $events1 \n\n/" cond001.txt
	sed -i "1s/^/NumEvents $events2 \n\n/" cond002.txt

	#------------------

	echo -e "\nPerforming motion correction using BROCCOLI \n"
	MotionCorrection bold_${subject}.nii.gz

	echo -e "\nPerforming smoothing using BROCCOLI \n"
	Smoothing bold_${subject}_mc.nii.gz

	echo -e "\nSegmenting T1 volume using FSL\n"
	fast T1_brain_${subject}.nii.gz
	mv T1_brain_${subject}_pve_1.nii.gz T1_brain_${subject}_gmseg.nii.gz
	mv T1_brain_${subject}_pve_2.nii.gz T1_brain_${subject}_wmseg.nii.gz

	echo -e "\nExtracting first fMRI volume\n"	
	cp bold_${subject}.nii.gz bold.nii.gz
	3dTcat 'bold.nii.gz[0]' -prefix bold_${subject}_first.nii.gz  

	echo -e "\nRegistering T1 volume to fMRI using FSL\n"
	epi_reg --epi=bold_${subject}_first.nii.gz --t1=T1_${subject}.nii.gz --t1brain=T1_brain_${subject}.nii.gz --out=fMRI_to_T1_${subject}

	convert_xfm -omat T1_to_fMRI_${subject}.mat -inverse fMRI_to_T1_${subject}.mat

	echo -e "\nTransforming gray matter mask to fMRI space using FSL\n"
	flirt -in T1_brain_${subject}_gmseg.nii.gz -ref bold_${subject}_first.nii.gz -applyxfm -init T1_to_fMRI_${subject}.mat -out bold_${subject}_graymatter.nii.gz
	flirt -in T1_brain_${subject}_wmseg.nii.gz -ref bold_${subject}_first.nii.gz -applyxfm -init T1_to_fMRI_${subject}.mat -out bold_${subject}_whitematter.nii.gz

	echo "Thresholding gray matter segmentation to get a binary mask"
	3dcalc -a bold_${subject}_graymatter.nii.gz -expr 'ispositive(a-0.5)'  -prefix bold_${subject}_graymatter_binary.nii.gz

	1d_tool.py -infile bold_${subject}_motionparameters.1D -set_nruns 1 -derivative -write bold_${subject}_motionparameters_deriv.1D

	# Heteroscedastic analysis

	./HeteroGLM bold_${subject}_mc_sm.nii.gz -designfiles regressors.txt -contrasts contrasts.txt -ontrialbeta allbetamd.txt -ontrialgamma allgammamd.txt -ontrialrho allrho.txt -mask bold_${subject}_graymatter_binary.nii.gz -regressmotion bold_${subject}_motionparameters.1D -regressmotionderiv bold_${subject}_motionparameters_deriv.1D  -saveoriginaldesignmatrix -savedesignmatrix -draws 1000 -verbose -savefullposterior -burnin 1000

	mv bold_${subject}* Results/RhymeJudgmentHetero

	# Homoscedastic analysis

	#./HeteroGLM bold_${subject}_mc_sm.nii.gz -designfiles regressors.txt -contrasts contrasts.txt -ontrialbeta allbetamd.txt -gammaregressors gammaintercept.txt -ontrialrho allrho.txt -mask bold_${subject}_graymatter_binary.nii.gz -regressmotion bold_${subject}_motionparameters.1D -regressmotionderiv bold_${subject}_motionparameters_deriv.1D  -saveoriginaldesignmatrix -savedesignmatrix -draws 1000 -verbose -savefullposterior -burnin 1000

	#mv bold_${subject}* Results/RhymeJudgmentHomo

done





