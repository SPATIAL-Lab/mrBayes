# mrBayes
Code for manuscript on blue-green water partitioning of the Mississippi River Basin (MRB). 

### code/
#
- **FINAL_mrb_ssn_object.R** Builds spatial stream network for d-excess 
- **FINAL_d_SSN_mrbcode.R** Builds spatial stream network model for d-excess 
- **FINAL_mrbayes_code.R** Code for Bayesian isotope mass balance model 
- **MRBayes_sumstats.R** Code for calculating ecohydrologic component summary statistics
- **FINAL_blockstats_mrbayes.R** Code tests if spatial patterns persist if posterior uncertainty is propagated
- **FINAL_blue-green_bi-plots.R** Code builds regional headwater vs downstream blue-green bivariate plots
- **FINAL_ScalingLaw_mrb.R** Code for fig 4 and scaling law climate coupling, aridity, effects on blue-green cycling
- 
### data/ 
#
All data files can be downloaded from doi:10.5281/zenodo.17545916

Data & folders only on doi:10.5281/zenodo.17545916
- **mrb.rivFcopy.zip** File for data used to build spatial stream network
- **mrbevap.ssn.zip** Final MRB Landscape Stream Network file for spatial stream network model
- **posterior_matrix.rds** Posterior matrix from MCMC sample parameters via MRB Bayesian isotope mass balance model
- **dstreams_bay_extended.RDS** MRB stream network with Bayesian isotope mass balance model component summary statistics


Data on mrBayes GitHub repository and doi:10.5281/zenodo.17545916
- **dstreams_bay.RDS** spatial stream network for d-excess
- **GP_raw.RDS** Great Plains polygon [EPA's SDAM GP](https://sdam-for-great-plains-eprusa.hub.arcgis.com](https://sdam-for-great-plains-eprusa.hub.arcgis.com/))
- **lowmrb.RDS** Lower MRB watershed polygon [HydroBASINS](https://www.hydrosheds.org/products/hydrobasins#downloads)
- **midmrb.RDS** Middle MRB watershed polygon [HydroBASINS](https://www.hydrosheds.org/products/hydrobasins#downloads)
- **ohiobas.RDS** Ohio River watershed polygon [HydroBASINS](https://www.hydrosheds.org/products/hydrobasins#downloads)
- **uppermrbbas.RDS** Upper MRB watershed polygon [HydroBASINS](https://www.hydrosheds.org/products/hydrobasins#downloads)
- **midlowmrb.RDS** Middle-Lower MRB watershed polygon [HydroBASINS](https://www.hydrosheds.org/products/hydrobasins#downloads)
- **mergedMRBpoly.RDS** MRB watershed polygon [HydroBASINS](https://www.hydrosheds.org/products/hydrobasins#downloads)
- **MRB_AI_mean_2013_2014.tif** MRB Aridity Index PET/P, annual average (years 2013 & 2014) from [TerraClimate](https://www.nature.com/articles/sdata2017191)
  
### out/
#
Figures generated from code

