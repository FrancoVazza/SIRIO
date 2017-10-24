# SIRIO

Sirio (SImulator of Radio Interferometic Observations) is an  IDL package to practice with the logic of radio imaging
First developed in ~2012 by F. Vazza & A. Bonafede @Uni Hamburg as a tool to practice with the principles of radio imaging with simple IDL routines. 

It is not meant for real production of mock radio observations, but to illustrate all most important steps in radio imaging, giving qualitatively correct answer. It is still in beta-version and being tested, but feel free to have a try it on it and report problems to franco.vazza2@unibo.it

SIRIO includes (not in order of appearance):
- reading of 2D sky model (single frequency);
- creation of UV coverage of a given antenna dataset;
- beam formation and visibilities;
- formation of the "dirty" image;
- thermal noise contribution;
- baseline interpolation;
- cleaning by iterations.


Compilation & execution

Download everything
IDL

>.r needs_good

>.r needs_good

>.r sirio

>.r sirio 

(both twice)

>sirio    

or 

>sirio,file_sky_model='input.fits',number_point_sources=20,live_plot='y', sigma_rms=1e-6     etc...to change default parameters


Default values of parameters and meaning
   file_sky_model='relic1.fits'  ;...input sky model 
   folder='/Users/francovazza/Downloads/SIRIO/'   ;...main folder containing SIRIO
   file_antenna='VLA.reg'   ;...file with atenna positions
   number_point_sources=2   ;...additional pointlike sources to be generated
   live_plot='y'  ;...plots on 'x' device while computing 
   max_iter=1e4   ;..max iterations in cleaning
   sigma_rms=8e-5 ;...Jy/beam
   hour=10  ;...hours of integration
   

Multiple windows should appear during the execution, and intermediate and final outputs will appear in /output.
The live plotting procedure gets slow for >256^2 datasets

Suggested tests:
- lowering the hour of exposition;
- lowering/increasing the noise per beam;
- removing pointsources;
- imposing less iterations;
- test different input files available in /input: 3 are 200x200 radio relic maps, one is a larger 1200x1200 dataset.

Known issue: some incompatibilities across different IDL versions have been reported, and can be fixed by a simple switch of functions. For example, if "gauss_smooth" (within cleaning_iterations()) does not exist on your version, replace with "smooth". Feel free to report any other problem. 


...Enjoy! 


