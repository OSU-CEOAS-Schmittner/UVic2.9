#if !defined O_carbon
#undef O_carbon_13
#undef O_carbon_14
#undef O_co2ccn_data
#undef O_save_carbon_carbonate_chem
#undef O_carbon_co2_2d
#undef O_carbon_uncoupled
#undef O_carbon_13_coupled
#undef O_carbon_14_coupled
#undef O_save_carbon_totals
#undef O_TMM_interactive_atmosphere
#endif

#if !defined O_carbon_13
#undef O_c13ccn_data
#undef O_carbon_13_coupled
#endif

#if !defined O_carbon_14
#undef O_c14ccn_data
#endif

#if !defined O_mobi
#undef O_mobi_nitrogen
#undef O_mobi_iron
#undef O_mobi_silicon
#undef O_mobi_caco3
#undef O_save_mobi_fluxes
#undef O_save_mobi_diagnostics
#endif

#if !defined O_mobi_nitrogen
#undef O_mobi_nitrogen_15
#endif

#if defined O_mobi_caco3
#undef O_kk_ballast
#endif

#if defined O_mobi_nitrogen && !defined O_mobi_o2
THIS is HERE TO GENERATE AN ERROR WHEN COMPILING: O_mobi_o2 must be ON when using O_mobi_nitrogen
#endif
 
#if defined O_mobi_iron && !defined O_mobi_o2
THIS is HERE TO GENERATE AN ERROR WHEN COMPILING: O_mobi_o2 must be ON when using O_mobi_iron
#endif
 

#if !defined O_PaTh
#undef O_save_PaTh_diagnostics
#endif