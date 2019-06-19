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

