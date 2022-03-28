THIS BRANCH

In this branch I am introducing the proper change of stellar spin due to the single and binary evolution processes.
I have introduced AngMomSpin evolution in:
- SSE property AngMomSpin, the spin evolve considering angular momentum mass loss through winds and magnetic braking
- BSE WindMass accretion, angular momentum accreted through the wind mass accretion
- BSE Tides, Spin evolution due to the tides, the evolution is limited by the reaching of the equilibrium spin
- BSE RLO, Donor loses AngMomSpin while accretor accretes it
- During change of track, the AngMomSpin is conserved
At the end of the evolution, OmegaSpin (the angular velocity) and Spin (Omega/Omegacrit) are updated. The AngMomSpin is limited
to the value AngMom_crit = Inertia * Omega_crit.

We have added the effective radius Radx=min(Radius, RL), therefore during the RLO, RL is used as an effective radius.
The effective radius is used in
- SSE/BSE Wind: The angular momentum lost by the star is equal to 2/3*Radx*Radx*DM. Since this value is lost in the SSE
evolution that does know nothing about the binary, in the BSE Windaccretion process the value lost during SSE (2/3*Radius*Radius*DM) is
corrected.
- Tides: for the tides the effective radius is of paramount importance, indeed da/dt, de/de, dSpin/dt in the Hut equations depends
on (R/a)**8, where a is the semimajor axis. However, since these equationa are bases on perturbation theory they are consistent
only when R/a<1. During the RLO, the stellar radius can become even ten times larger the semimajor axis, therefore the effective radius
is fundamental. In order to estimate the total angular momentum lost we also rescale the Inertia by a factor Radx^2/Radius^2, indeed
Inertia propto R^2.
- RLO: In the RLO the angular momentum lost by the primary is estimated as DM*Radx*Radx*OmegaSpin

NOTICE: We do not take into account Radx for the magnetic braking, but this is exactly what is done in BSE/MOBSE



**3.0alpha ongoing**
- Fixed a typo in Orbital_change_RL::kelvin_helmotz_tscale
- Collision during the RLO can now be disabled or enabled at runtime (disabled by default)
- Added parameter rlo_max_nuclearmt to set the maximum on the normalisation of the nuclear mass transfer
- Added SNmodel: directcollapse
- Added an optional parameter to allow for unstable MS and pureHEMS mass transfer (rlo_mtstable_ms, ture by default)
- Set a bug on the output of the event "Swallowed"
- New qcrit from startrack
- Added a parameter to set a minimum kick velocity (sn_min_vkick set to 0 by default)
- Added Sn kicks models CC15 and EC15CC265
- Fixed a typo during the merger of a pureHE star with an Hstar  + refactoring of Mix process
- Added the Process Kollision to handle collision at periastron (hurley and disabled are the possibile options), during the RLO
the collision at periastron is checked only if the new parameter rlo_enable_collision is set to true (it is false by default)
- Added the parameter sn_kick_velocity_stdev to set the standard deviation of the Maxwellian used to draw random kicks in some
SN Kicks models (e.g. Hobbs and Unified)
- Added the parameter star_tshold_WR_envelope to check WR condition, a star is considered WR if MHE>=(1-star_tshold_WR_envelope)*Mass
- A WR star (see above) is now automatically forced to jump to pureHe tracks
- In SEVN to BSE phase translation now nakedHe in Hurley are equivalent to pureHe in SEVN (i.e. stars that are following pureHe tracks)
- Added parameter ts_check_NSspin to enable or disable the check on the adaptime timestep for NS spin during mass accretion
- Added new qcrit implementation radiative_stable, the mass transfer is always stable for case BB, MS stars and HG stars.
- Added the runtime parameter optimistic_scenario_hg, it True the HG can survive the CE
- Change of runtime parameter names:
	- sn_compact_Mremnant_average_NS -> sn_Mremnant_average_NS
	- sn_comapact_Mremnant_std_NS -> sn_Mremnant_average_std_NS
	Previously the parameter were only related to the SN compact model, now they are shraed between all the SN models that need it.
- Added new sn models:
	- delayed_gauNS: same as delayed but the mass of the NS are drawn from a Gaussian with parameter mean: sn_Mremnant_average_NS and std: sn_Mremnant_std_NS (minimum NS mass 1.1 Msun), ECSN also have random Mass
	- rapid_gauNS: same as rapid but the mass of the NS are drawn from a Gaussian with parameter mean: sn_Mremnant_average_NS and std: sn_Mremnant_std_NS (minimum NS mass 1.1 Msun) ECSN also have a random mass
	- deathmatrix: implementation of the Death Matrix by Woosley+20 (https://arxiv.org/pdf/2001.10492.pdf, Tab. 2), the mass of the Remannt depends on the preSN MHE and it is takedn from a look-up table



-minor bug fix

**2.15alpha**
- New method in Star to check for Wolf-Rayet, the condition is MHE>0.98*Mass
- Change in SEVN-BSE phase conversion, now the BSE phase 7 is set if the star is a WR and
in SEVN phase 4, phase 8 is set if the star is a WR in SEVN phase 6
- Kicks Refactored, now the method apply is a wrapper defined in the abstract class and the _apply
is a pure virtual method to override.
- Added some checks in Supernova
- Changes on Lambda, the switch to star_lambda_pureHe is done when BSE phase is 7,8,9 not just when star is on  pureHe tracks.
  This is what is done in BSE/MOBSE where a stars can become 7 even through SSE.
- Mix Process refactored
- BUGFIX: Before when star_lambda>0, the lambda was always estimated using star_lambda_pureHe
- BUGFIX: Solved a bug triggered in the merging of a pureHE with a Hstar without He core

**2.14alpha**
- BUGFIX: Fixed a bug in compactness::explosion, before there was not a check for remannt with 0 mass due to the ppisn
- BUGFIX: Fixed a bug in the script run_sse.sh
- BUGFIX: Fixed a bug that could draw the same natal kick for both stars of a system. This was due to a problem with the random number
          generator that was reset each time a change of track was triggered

**2.13alpha*
- BUGFIX: we fix a bug causing the stars to jump to tracks with different metallicities and/or sn options.
    This was caused by  a mismatch between Z and snmode set in the input list and the same options set in runtime parameters (Z and snmode)
    When a star jump to a new track it gets the properties as Z and snmode from old track, but due to the bug the parameters were not set
    considering the effective Z and snmode as set from the runtime option but the ones in the input list were used.
- Added the runtime option spin to overwrite the initial spin of all the data
- Made a change in SEVN-BSE conversion, now SEVN phase>4 are not automatically transformed to BSE phase 5. Rather
  the fraction of convective envelope is checked, if >0.4, the BSE phase is 5 otherwise it remains 4.
- Minor change on the Nanjing lambda to consider the BSE phases instead of the SEVN phases


**2.12alpha*
- Added the option to use lambda estimate from Xu&Li10

**2.11alpha*

- Added the option to use lambda estimate from Klencki+21
- Added the new parameter star_lambda_pureHe to set a constant lambda for pureHe star
  (notice that the available lambda options do not have fit for pureHe stars)
- Fixed a bug on retrieving the path of the SEVN folder, it failed if the main folder name was different
  with respect SEVN


**2.10alpha*

- Total refactoring of qcrit, now MTStability it is used instead, the input parameter rlo_qcrit
is changed with rlo_stability
- Fixed some typos in Lambda ad added new options

**2.9alpha*

- Added the eddington accretion limit on degenate remnants  (WD, NS, BH) during the Wind mass accretion
- Changed make file and added the function make_unique in utility to restore compatibility with gnu compiler version 4.8.
- Solved a bug in correct_interpolation_error MHE, MCO

**2.8alpha**

- Complete refactoring of the remnants. Now we introduce a new class Staremnan, for each remnant there is an inherited subclass. These classes contains all the
info to evolve the remnants. Each obejct of class star as a pointer to a class remnant, this is equal to a nullptr at the beginning
at it is set after a call to SN. Then evolve remnant and set remnant get the quantities directly from the remnant object inside star.

- Added the property OmegaRem (angular velocity for the remnants), Bmag (magnetic field)

- Added a new qcrit option: Hurley_Webbink_Shao: classical qc from Hurley_Webbink with a stability criterium based on Shao+21 (https://arxiv.org/pdf/2107.03565.pdf)
when the accretor is a BH.

- Added a new proposed timestep condition to reduce the timestep when a star is becoming nakedHE or nakedCO

- Added two new boolean parameters: ev_set_maxCO and ev_set_minHE (false by default). If ev_set_maxCO=true
the first time a star develops a CO core, we set the maximum CO core Mass for SSE as the last value of the interpolating tracks.
If ev_set_minHE=true the first time a star develops a CO core, we set the minimum HE core Mass for SSE as the last value of the interpolating tracks.

**2.7alpha**

- Changed the function gen_rseed to solve the problems of duplicates seeds

- Added the properties Event and BEvent to store sse and bse events (e.g. phase change, merger, rlo etc). Enabled
the option for dtout "events". It this option is chosen, the output will contains only the timestep with a triggered event

- Spin evolve following the conservation of angular momentum and taking into account the wind mass loss.

- Tides have been refactored to take into account the info of the convective layer in the tables (Qconv, Depthconv, Tconv)
If the tables are not present the evolution of Qconv and Depthconv is now consistent with the BSE assumptions

- Now the error during initilisation are properly handled.  Inside chunck dispatcher during the initialisations of systems inside the vector
there is try-catch statement of sevnio error, if an error is catched the system is not inserted in the vector and the input properties are output in failed_initialisation.dat.
The parameter initerror_stop has been handled to control the initilisation errors: if true halt the programme when a initialisation error has been thrown, otherwise neglect the sytem and continue. The file failed_initialisation will contains all the sytems for which the initialisation failed.

**2.6alpha**

- Added a lot of new log prints

- Added Merger log print (it is used both in Mix and CommonEnvelope after coalescence)

- Added new qcrit prescriptions (Neijssel10, COMPAS, COSMIC_Claeys)

- Transformation from SEVN to BSE phase now takes into account BSE phase 0

- New conditions to match the core when chaing tracks:
  - First we try to match the innermost core.
  - If the innermost core was CO and the match fails we try with MHE
  - If also the second try to match the core fails and only if we are looking for a new track after a merge we try to match the total mass

- Number of bug solved

- Improved handlinf of NakedCO and Naked heilum

- Fixed a bug adding the new concept of auxiliary_star

  The bug was caused when we added the possibility to overwrite the input
  parameters tini and Z using global options at runtime. This  caused a bug
  when we created auxiliary stars (e.g. during the change of tracks on in the
  Lambda estimate) because  indepently by the choice of Z and tini in the
  constructor, these values were always overwritten by the global value.
  In order to solve this issue, we added the auxiliary_star paradigm:
  - We add a boolean Star member auxiliary_star (initiliased to false), if it is true the class
  instance should be considered auxiliary and used only for special purpose.
  - The only way to create an auxiliary star is using the two special constructor
  that have among parameters a pointer to another star. In these constructor
  auxiliary_star is set to true.
  - The methods inspect_param_tini and inspect_param_Z have been modified with
  an additional parameters ignore_global (default False), it is true the global
  paramet is ignored.
  - If the star is auxiliary in init1 and init2 the local variable ignore_global is set to true
  and passed to inspect_param_tini and inspect_param_Z
  - TODO: This new paradigm is OK and works, but it is a bit involved a hard to mantain, it will be
  much better to create a inherited class Star_auxiliary that can directly handle the problem (we add a placeholder in Star.h).

- Added placeorder options for Z (Z), SN (snmode), Tstart (tini), Tend (tf) and Dtout (dtout)
In the list these values can be set to xxx, nut in this case is mandatory to not use list as option
for the global parameters

**2.5alpha** 29/03/2021
- Fixed a bug on random kick and reproducibility
- Added the possibility to have optional tables. In particular now Inertia, RHE and RCO are optional
- Added the possibility to initialise a compact object, it is only needed to had a label after the mass in the input tables:
e.g. 30BH will initiliase a BH with 30 Msun, the other label are NS for neutron stars, NSEC, for Electorn Capture neutron star
HEWD for Helium White Dwarf, COWD for carbon-oxygen white dwarf, ONEWD for Oxygen-Neon white dwarf.
- Added the possibility to initialise an object with 0 Mass, it will be initiliased as an Empty remnant.
- Added tracks for Mist and parsec from Spera+19
- Minor bugs fixed


**2.4alpha**  27/01/2021
- Added new pureHE tables
- Cleaned by warning messages
- New compile and run scripts

**2.3alpha** 20/01/2021
- Added  the option to use a list of values in dtout.
- Now the evolution is forced to evalue the state of the star/binary at the time chosen using the parameter dtout

**2.2alpha** 28/12/2020
- Added functors EvolveBinaryCompact and EvolveBLC
- New header evolve.h. It contains functions to evolve a list of stars.
- The evolve functions are implemented using functor so that they can be passed to a general chunk evolve function.


**2.1alpha** 26/12/2020
- Added the option to set a global value for the SN model (-snmode ..) overwriting the option in the list
- Memory leaks fixed
- Added the option to evolve stars and binaries in chunk
- Added new SN models:
	- rapid
	- compactness (Mapelli+20)
	- compactness with probabilistic choices for explosion/implosion  (Mapelli+20 and Patton&Sukbhold20)