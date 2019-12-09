# MIPLAB brain atlases

This repository includes the high-resolution, age-specific FLAIR and non-contrast CT atlases of the elderly built from clinical images. The atlases provided are compatible with the standard [MNI atlases](http://www.bic.mni.mcgill.ca/ServicesAtlases/ICBM152NLin2009). 

Keywords: Fluid-attentuated inversion recovery, Non-contrast computed tomography, structure-preserved interpolation, group atlas. 

## Data provided 
* FLAIR atlases
  * Asymmetric 
    * Head (with skull) 
    * Brain (skull stripped)
    * Mask 
  * Symmetric
    * Head (with skull) 
    * Brain (skull stripped)
    * Mask   
* NCCT atlases
  * Asymmetric brain (skull stripped)
  * Symmetric brain (skull stripped) 
* Deformations 
  * MIPLAB to MNI
  * MNI to MIPLAB
  
## Usage Notes 
The atlases are modality specific for stroke research. They are best used in neuroimaging analysis of the elderly cohort (mean:69, std: 13.5). Since both the NCCT and FLAIR atlases are in the same physical space, they allow researchers to include both a multi-modal dataset if the study design permits the same. The atlases are compatible to use with the follwoing neuroimaging packages: 

* Registration:
  * [ANTS](https://github.com/ANTsX/ANTs)
  * [NiftyReg](http://cmictig.cs.ucl.ac.uk/wiki/index.php/NiftyReg)  
* Analysis:
  * [ANTS](https://github.com/ANTsX/ANTs)
  * [FSLutils](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Fslutils)
  * [Brain extraction toolkit](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BET)
  
## Citation 
  






