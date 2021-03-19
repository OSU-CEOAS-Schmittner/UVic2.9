Mar 18, 2021: 08 updates corresponds to MOBI2 with silicon cycling.
1) Si cycling simplified from Karin Kvale's model in that it does not include external sources and sinks. 
Opal has been included as an additional prognostic tracer.
Si is conserved as can be verified by adding o_sil and o_opl.
Option: O_mobi_silicon

2) Coccolithophores have been removed as a separate plankton functional type. CaCO2 production and dissolution now depend on the calciate saturation state omega.

3) Option names have been changed:
O_npzd -> O_mobi
O_save_npzd -> O_save_mobi_fluxes
O_npzd_extra_diagnostics and O_npzd_iron_diagnostics -> O_save_npzd_diagnostics

4) All half-saturation constants (N, P, Fe) are variable as in Aumont et al. (2015)

Nov 2, 2020: 08 updates from Samar that include Karin Kvales diatom and silicon code
Here are notes from Samar's email:
1) This introduces two new tracers silica and diatoms which are switched on via O_kk_si and 
O_kk_diat, respectively. They are nominally independent of each other and in principle one can 
be switched on and off without regard to the other. In practice, this is not applied uniformly in the 
code and I don't know what would happen if you were to do that.

2) Silica can be either prognostic or read in as a mask (by undefining O_kk_si).

3) Other related CPP options:
O_kk_variable_sipr
O_kk_si_compensating_sources

4) To look for modifications do a diff with 07 or grep for "kk_" in 08 or in Karin's version.

5) Where I've done something different from Karin's version, adapted MOBI code where there 
wasn’t an equivalent in her’s, or had a question I've marked those places with SPKKK. You can 
grep for this.

6) Unfinished business:
- mobi_src is not complete as I mentioned above.
- the river silica business won't currently work with the TMM but I can fix that later. Everything else 
should work offline.
- silica is not yet passed to co2calc as is done in Karin's version (and it seems used to be done in 
02 but we've removed that in MOBI). This is very easy to do but I first wanted to get the silica tracer 
working before we mess with carbon.

# updates
Jun 19, 2019: Andreas added updates 04-06 and copied 06 to latest
The 06 updates were created by Samar Khatiwala who added modifications to the light limitation from Juan Muglia and sends these notes:
I’m attaching the latest version of MOBI that I hope we can all use from now on. This is created as a new updates 06 that includes 
everything from 05 + Juan’s latest npzd_src.F + my modifications. The latter includes the TMM interface, hooks to extract TMs, a 
bug fix in setmom.F, and a new subroutine mom/set_sbc.F to set boundary conditions. As I said in a previous email this avoids 
having to modify the code with lines like "if (issalk .ne. 0) …” every time a new tracer is added. The way I’ve done this is to add a 
new variable trsbcindex which maps the tracer index to the sbc index for that tracer. Whenever you add a new tracer, you should 
add a corresponding line in subroutine sbc_init, e.g.:

#if defined O_npzd_alk
     call set (issalk, m, mapsbc(m), 'ssalk', m)
     trsbcindex(ialk) = issalk     <—————————————————  LIKE THIS
     call set (ialkflx, m, mapsbc(m), 'alkflx', m)
#endif

Long sequences of statements like "if (issalk .ne. 0) …” can then be replaced with a single do-enddo loop over the tracers 
(see setmom.F and tracer.F for examples).

I have verified that running this with the online model reproduces (to within the initialization bug in setmom.F) Juan’s configuration.
------------
- `01/` is copied from `UVic_ESCM.2.9.updated.tar.gz`, which is not avaiable anymore
- `02/` is copied from `UVic_ESCM.2.9.updated.tar.gz`, which is not avaiable anymore
- `03/` is copied from [`http://kelvin.earth.ox.ac.uk/spk/Research/TMM/Kiel_Jan_2017_updates_to_UVIC2.9.tar.gz`](http://kelvin.earth.ox.ac.uk/spk/Research/TMM/Kiel_Jan_2017_updates_to_UVIC2.9.tar.gz), which is part of [`https://github.com/samarkhatiwala/tmm`](https://github.com/samarkhatiwala/tmm)

