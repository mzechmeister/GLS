# Download all files

svn export https://github.com/mzechmeister/GLS/trunk/fortran GLS
cd GLS

# Compile with

make

# or if it doesn't work, try e.g.: f77 GLS.f -o GLS.e
# or whatever your fortran compiler is (g77, gfortran, ...).
# This creates the executable file GLS.e
# (The make command tries also to create a command GLS in your ~/bin that
# allows to execute GLS.e from every working directory.)

# For a short overview of the command line options type:

GLS -?

# Example to start the program:

GLS

# or (in the program directory) ./GLS.e, respectively
# and you will be asked for the inputs.
# For the start and end frequency there will be suggestions, which are
# a frequency close to zero (a period longer than the time base) and
# the averaged Nyquist frequency. Just press enter to use these values.
# Other frequency limits can be also specified here.

# You can also pipe the inputs as arguments via the command line like:
# GLS file freqbeg freqend starmass
# Example for a program start with piped arguments:

GLS GJ1046.dat 0.001 0.008

# which is probably be more convenient, when the program is often used.
# Best fit values are outputted on the screen and in the GLS.log file.
# The graphical output should work well with gnuplot 4.2.3. It can be suppressed with -s (silent) option.
# Unfortunately with gnuplot v4.2, only mouse zoom-in and keypress un-zoom ('u') is
# interactively possible in the X11-plots which run via a pipe in background.
# However, the windows can be closed very fast by pressing any key.
# kpdf or okular is recommended to view the pdf-files (GLS.pdf,GLSRV.pdf), because they can watch for updated files.

GLS -r

# recalls the last plots.

# The GLS and the Keplerian periodogram are stored in GLS.plt and GLSKep.plt
# (If you want for the GLS the analogon to the normalisation that is commonly used in the
# Lomb-Scargle just multiply p with the number of data points N: P = p*(N-1)/2.)
# If you are also interested in the Keplerian periodogram, which can take an extensive amount of computing time,
# you can use the -k option.

# The example of GJ1046 is taken from Zechmeister & Kuerster (2009), if you wish to compare.

# The data file should be an ascii file of the form:
# JD RV RVerr
# The program assumes that the time is given in JD, i.e. days and the data are radial velocities in m/s.
# Otherwise labels and companion masses may be wrong.
# If the data have no errors or errors, e.g.
# JD RV
# it will be recognised by the program. Use the option -v to perform unweighted fitting
# and to suppress the warning about missing error values.

# The program looks in the working directory for the parameter file GLS.par
# to adjust user specified step sizes for e, T0, and freq.
# If not found, default values are used.
# Note: For (very) eccentric orbits the default values may be not sufficient! Then adjust a denser grid size in GLS.par.

# Tips:
# 1) The residuals of the best fits you find in RVSinRes.plt and RVKepRes.plt
#     (iRVSinRes.plt and iRVKepRes.plt resp. in refine mode).
#    You can use (a duplicate of) them again as input data file
#    for a periodogram analysis of the residuals.
# 2) Use GLS to replot and refine ("improve") the parameters with a Marquardt-Levenberg algorithm.


# The description of the generalised Lomb-Scargle periodogram you find in Zechmeister & Kuerster (2009).
# Please include this reference when you publish results with it.


# Have a lots of fun!


# by Mathias Zechmeister (2012-08-03)

# file descriptions:
# GLS.f        source code
# GLS.e        exe-file
# ~/bin/GLS    executing main script to handle the input and graphical outputs
# GJ1046.dat   example data input file
# GLS.par      input parameter file
# GLSplot.gnu  gnuplot script to plot periodogram
# GLSRV.gnu    gnuplot script to plot RV time series
# GLSRVandPhased.gnu gnuplot script to plot RV time series
# fitsin.gnu   gnuplot script to refine the sinus parameters
# fitkep.gnu   gnuplot script to refine the Keplerian parameters

# generated outputs
# GLS.log      fit results of GLS
# GLS.plt      table GLS periodogram (including windowfunction and LS periodogram)
# GLSKep.plt   table Keplerian periodogram
# GLSFit.plt   table phase RV curve (GLS and Keplerian)
# RVSinRes.plt time series RV residuals of the sinusoid
#              (including the phase of the data points)
# RVKepRes.plt time series RV residuals of Keplerian orbit
#              (including the phase of the data points)
# iRVSinRes.plt time series RV residuals of the sinusoid (refine mode)
#              (including the phase of the data points)
# iRVKepRes.plt time series RV residuals of Keplerian orbit (refine mode)
#              (including the phase of the data points)
# GLS.ps       plot GLS/Keplerian periodogram
# GLS.pdf      plot GLS/Keplerian periodogram
# GLSRV.ps     plot RV time series
# GLSRV.pdf    plot RV time series
# fit.log      a gnuplot fitting log-file
# GLSfit.log   log file for refined parameters
# GLSEinit.tab interpolation table to solve fastly Kepler's equation
# zoom.gnu     small gnuplot script to scroll and zoom with keyboard keys
# pausepatch.gnu small gnuplot patch for v4.2.3

