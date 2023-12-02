# ICDM
Ideal Conductor / Dielectric Model for finite size corrections in ion transport as described in [Haji-Akbari A, Shoemaker B. Ideal conductor/dielectric model (ICDM): A generalized technique to correct for finite-size effects in molecular simulations of hindered ion transport. ChemRxiv. Cambridge: Cambridge Open Engage; 2023 doi:10.26434/chemrxiv-2023-1mgbh]

-----------------------------------------------

NOTE: This project requires access to the file "dielectricInterfacesGeneralized.py" which can be downloaded here:
https://github.com/shoemb/dielectricInterfacesGeneralized

-----------------------------------------------

The correction is run by calling the function "runCorrection" in the file main.m. For demonstration purposes, the correction is run for two cases: the initial passage of chloride as the first ion to leave the feed, and the subsequent passage of a sodium ion. "runCorrection" takes two arguments: an input filename, and an output filename. 

The input file consists of a series of variable definions of the form "varName = value" with one definition on each line. They can appear in any order. Input files for the example cases are in the directory "sampleInputs". The meaning of the inputs is as follows:

**totalSystemCharge** is the net charge of the entire simultion box (typically zero)

**chargeTransitingIon** is the charge of the transiting ion (in units of fundamental charge) on which the correction is being applied

**startingZ** and **feedStartingZ** are the positions (in Angstrom) where the correction is first applied and where the transiting ion has fully left the feed such that it produces its own induced charge in the feed. For passage of the initial ion these will typically be the same value.

**free-energy-profile** is the path to a txt file which contains the uncorrected free energy profile. This file contains two columns separated by a comma where the first column is the position of the transiting ion (in Angstrom) and the second column is the free energy (in units of kT). If two or more independent free energy profiles are available through different trajectories/configurations, the keyword can be used multiple times. In this case, error bars in the free energy profile can be determined.

**OPValues** gives values of the order parameter (in Angstrom) at which empirical measurements of the charge density of non-transiting ions outside the feed are available. If the transiting ion is the only ion outside the feed compartment, this line can be omitted from the input file.

**chlorideDensityFilename** path to a csv file which contains the empirical charge density of non-transiting ions outside the feed. The first two columns give the lower and upper bound of each observation bin (in Angstrom), the third column gives the midpoint of each bin (in Angstrom), and each subsequent column gives the charge density at particular values of the order parameter (in e/Angstrom^3). These milestones must align with those provided in **OPValues**. The first row must contain non-numeric entries as headers.  If the transiting ion is the only ion outside the feed compartment, this line can be omitted from the input file.

**endMembrane** is the location (in Angstrom) of the end of the membrane on the filtrate side. This is used to determine whether bins in the empirical charge density of non-transiting ions correspond to points within the membrane (and thus point charges) or points within the filtrate (and thus charged slabs). If the transiting ion is the only ion outside the feed compartment, this line can be omitted from the input file.

**z_c** is the location of the conducting surface on the feed side (in Angstrom)

**z_c_err** is the standard error in the measurement of z_c.

**shiftZc** can be set to "up", "down", or "none". The first two options will increase or decrease z_c by the value specificed in z_c_err in order to assess the sensitivity of the correction to measurement uncertainty in z_c.

**shiftUncorrected** can be set to "up", "down", or "none". The first two options will increase or decrease the uncorrected free energh profile by the standard error. Note that this requires providing multiple independent free energy profiles with the **free-energy-profile** tag.

**L_x** is the length of the simulation box (in Angstrom) in the x-direction.

**L_y** is the length of the simulation box (in Angstrom) in the y-direction.

**id_min** is the line number (numbering begins at 1) in the file provided by **free-energy-profile** which corresponds to the order parameter value where the correction is first applied.

**id_max** is the line number (numbering begins at 1) in the file provided by **free-energy-profile** which corresponds to the order parameter value up to which the correction is applied. This can be set to the maximum line number but must at least extend beyond the transition state.

**z_piston** is the position (in Angstrom) of the piston on the feed side. 

**boundaries** gives a list of z-positions (in Angstrom) of the interfaces between different conducting/dielectric regions.

**feedRegionIndices** gives a list of indices for regions defined by **boundaries**. In the provided example, boundaries is set to [-2,8.7]. This means there are three regions ([-inf,-2],[-2,8.7],[8.7,inf]) which are identified by indices 0,1,2 respectively. This input is necessary for computing the total induced charge in the feed.

**epsilon** gives a list for the dielectric constants in each region. This list will have one more entry than **boundaries**.

**numImageChargeIterations** is the number of iterations of the algorithm which generates charges using the method of images. Larger values give better accuracy but also increase the runtime.

**XImagesOfPointCharges** and **YImagesOfPointCharges** describe the number of periodic replicates of each image charge to include. The periodic replicates will form an array of extent [-XImagesOfPointCharges * L_x,+XImagesOfPointCharges * L_x] in the x-direction and  [-YImagesOfPointCharges * L_y,+YImagesOfPointCharges * L_Y] in the y-direction.

**targetXBinSize** and **targetYBinSize** set the bin size for numerical integration of a charged finite slab to determine the potential it generates. If necessary, the actual bin size will be adjusted to be slightly smaller than the target values so that an integer number of bins spans the simulation box.


The output file contains four rows. The first two rows give the order parameter values (in Angstrom) and free energy values (in kT) prior to the correction being applied. The third and fourth rows give the order parameter values (in Angstrom) and corrected free energy values (in kT) after applying the method.
