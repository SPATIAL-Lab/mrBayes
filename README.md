# mrBayes
Code for manuscript on blue-green water partitioning of the Mississippi River Basin (MRB). 

### code/
#
- **FINAL_mrb_ssn_object.R** Builds spatial stream network for d-excess 
- **FINAL_d_SSN_mrbcode.R** Builds spatial stream network model for d-excess 
- **FINAL_mrbayes_code.R** Code for Bayesian isotope mass balance model 
- **MRBayes_sumstats.R** Code for calculating ecohydrologic component summary statistics
- **FINAL_blockstats_mrbayes.R** Code tests if spatial patterns persist if posterior uncertainty is propagated
### data/ 
#
Data file must be downloaded from doi:10.5281/zenodo.17545916
- **dstreams_bay.RDS** spatial stream network for d-excess
- **GP_raw.RDS** Great Plains polygon [EPA's SDAM GP](https://sdam-for-great-plains-eprusa.hub.arcgis.com](https://sdam-for-great-plains-eprusa.hub.arcgis.com/))
- **lowmrb.RDS** Lower MRB watershed polygon [HydroBASINS](https://www.hydrosheds.org/products/hydrobasins#downloads)
- **midmrb.RDS** Middle MRB watershed polygon [HydroBASINS](https://www.hydrosheds.org/products/hydrobasins#downloads)
- **ohiobas.RDS** Ohio River watershed polygon [HydroBASINS](https://www.hydrosheds.org/products/hydrobasins#downloads)
- **uppermrbbas.RDS** Upper MRB watershed polygon [HydroBASINS](https://www.hydrosheds.org/products/hydrobasins#downloads)
- **midlowmrb.RDS** Middle-Lower MRB watershed polygon [HydroBASINS](https://www.hydrosheds.org/products/hydrobasins#downloads)
- **mergedMRBpoly.RDS** MRB watershed polygon [HydroBASINS](https://www.hydrosheds.org/products/hydrobasins#downloads)
- **mrb.rivFcopy.zip** File for data used to build spatial stream network 
- **posterior_matrix.rds** Posterior matix for MRB Bayesian isotope mass balance model
- **mrbevap.ssn.zip** Final MRB Landscape Stream Network file for spatial stream network model
