# UVic2.9 with updated Marine Iron Biogeochemistry
This version contains the new marine iron biogeochemistry modifications from Somes et al., (2021; https://doi.org/10.1029/2021GB006948), which has now been implemented into MOBI2.1 updates level 08. The modified files are located in updates/latest directory. The new iron biogeochemistry modifications can be enabled with the following options in run/mk.in (O_mobi_iron_var_ligands, O_mobi_iron_sed_dale, O_mobi_iron_gesamp_atmfedep, O_mobi_iron_inscav_nonlinear, O_mobi_sedcox_flogel, O_mobi_omz_threshold_smooth). Please note some marine biogeochemical parameter changes in run/control.in. 

# UVic2.9
This is the base code of the University of Victoria (UVic) climate model version 2.9 used at OSU. The source directory should not be changed. It is the original code without updates. It resides in /usr/local/models/UVic_ESCM/2.9/. Changes should be made to "updates/latest", which contains the latest updates. It should reside in /usr/local/models/UVic_ESCM/2.9/updates/. 

## Further info
* [OSU-UVic2.9 webpage](https://osu-ceoas-schmittner.github.io/UVic2.9/)
* Model of Ocean Biogeochemistry and Isotopes [MOBI](https://github.com/andreasschmittner/UVic2.9/wiki/Model-of-Ocean-Biogeochemistry-and-Isotopes-(MOBI))
* How to use git and github with this code: [OSU-UVic2.9 wiki](https://github.com/OSU-CEOAS-Schmittner/UVic2.9/wiki)
