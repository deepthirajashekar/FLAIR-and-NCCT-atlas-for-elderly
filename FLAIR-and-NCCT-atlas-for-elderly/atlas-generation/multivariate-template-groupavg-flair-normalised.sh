#!/bin/bash 
working_dir='/media/deepthi/DEEPTHI/miplab-atlas/ct/ct-1mm-aspects10'
cd ${working_dir}
data_dir='niftiis' 
target_dir='interimsHU'
results_dir='templatesHU'

#reference_img=${working_dir}/${results_dir}/ESCAPE_site07_011_Admission_Masked.nii.gz

cd ${working_dir}/${data_dir}
# ${ANTSPATH}/antsMultivariateTemplateConstruction2.sh -d 3 -c 2 -j 8 -i 4 -k 1 -y 0 -r 1 -a 0 -q 30x50x20 -m CC -o ctHU_ -b 1 *.nii.gz
${ANTSPATH}/ctMeanIntensityMultivariate.sh -d 3 -c 2 -j 8 -i 4 -k 1 -y 0 -r 1 -m 30x50x20 -s PR -o ctHU_ -b 1 *.nii.gz 

mkdir -p ${working_dir}/${target_dir}
mkdir -p ${working_dir}/${results_dir}

mv ${working_dir}/${data_dir}/ctHU_*.nii.gz ${working_dir}/${target_dir}/.
cp ${working_dir}/${target_dir}/ctHU_*template.nii.gz ${working_dir}/${results_dir}/.


