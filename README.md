# advectus
Molecular movement quantified - it not just diffusion anymore....

The transport of molecules within the cell occurs not just by diffusion but also by flow-mediated processes (i.e., advection)


This repository contains Matlab scripts* used in "Compartmentalized Cytoplasmic Flows Direct Protein Transport to the Cell's Leading Edge."
Catherine G. Galbraith, Brian P. English, Ulrike Boehm, and James A. Galbraith
bioRxiv [Preprint]. 2024 May 14:2024.05.12.593794. 
doi: 10.1101/2024.05.12.593794.PMID: 38798549


The code has been tested on Matlab versions 2021 and above on both Windows and Mac OS.  Install and run times are minimal on either operating system.  

Supplied as is with no guarantees as to coding efficiency.
MIT license - Copyright © 2024  [GalbraithLab - JA Galbraith, CG Galbraith]


1) TransportRatio.m
   
 To emulate the hydrodynamic experiment where a dye is continuously injected into a liquid medium, a diffraction-limited beam was used to continuously photoactivate photoactivatable (PA) fluorescent proteins, and the dispersal of fluorescence was monitored over time. The fluorescence intensity along a line through the activation (injection) point was recorded to determine the amount of material transported in each direction away from this point parallel and perpendicular to the cell edge.  

This script identifies the center peak of the fluorescent intensity curve and calculates the area under the curve on each side of the peak.   This is done for both the raw data curve and a spline-smoothed trace.  If dispersion is only due to diffusion, the curve is symmetric, and the areas on either side of the peak will be equal, with a ratio of one. However, if there is an advective component to the transport, more material will accumulate on the "downstream" side, leading to a ratio greater than one.  Input is either a .xlsx or .txt file with data in columns.  The first row of each column is the experiment identifier.  Output is the calculated ratio.  

 
2) MSScalc.m

This script calculates the Moment Scaling Spectrum (MSS) from single molecule track data.  The analysis is based on the works of Ferrari, R. et al. (2001), Physica D 154, 111–137, and Ewers, H et al.  (2005), PNAS  102 (42), 15110-15115.  The MSS calculates scaling exponents for the higher orders of the mean displacement curve, with the commonly used mean square displacement (MSD) being the second order.  The slope of the exponent curve versus the moment order is used to characterize the type of motion, with a slope >0.5 being directed motion, 0.5 being diffusion, and <0.5 being restricted motion.  

Tracking data from single molecule imaging experiments was obtained using the Matlab program TrackIt (https://gitlab.com/GebhardtLab/TrackIt), but any program that generates XY coordinates can be used.  MSScalc.m uses polyfitZero.m, a function by Mark Mikofski, and can be obtained at: https://www.mathworks.com/matlabcentral/fileexchange/35401.  The input data set is a .mat file containing two variables:  "TrackItData" - Trajectory ID Number, X, Y, Step #, ROI (not presently used) and "file" – the name of the dataset which is used for labeling output files.  Negative diffusion coefficients are ignored, and a cutoff for the coefficient of regression can be specified to exclude poor fits.  The output consists of three graphs: a histogram of the MSS slope, a histogram of the diffusion coefficients, and a scatter plot of the diffusion coefficient versus the MSS slope.  The scaling exponents for each moment and the slope of the MSS curve are saved as .mat and .csv files.  
  

3) FCSfit.m and TwoSpeciesFlow.m

This script fits the autocorrelation data supplied by the Leica TCS SP8 Falcon FLIM/FCS to a model that incorporates two species with advective flow from the work of Köhler et al. (2000)  Journal of Cell Science 113 (22): 3921–3930.  The triplet state is included in addition to the diffusive and advective components.  FCSfit.m requires the Matlab Curve fit toolbox and the TwoSpeciesFlow.m function, which provides the equation for the model.  FCSfit.m reads a .xlsx file containing the Leica autocorrelation data where Column A is time (ms), and the others are the data traces.  The first two rows (formatted as text variables) are the experiment identifiers - row 1, trial/run number, and row 2, experimental date.  The sheet name contains the experimental condition/treatment and is added to the filenames and graph titles.  The user can choose the desired tab to analyze if multiple sheets are present (for different conditions).  The output consists of graphs (saved as tiffs) for each trace, with the model fit overlayed on the raw data and a plot of the residuals.   A .mat file and separate CSV files contain the data, model fit, seeds, limits, residuals, and std of residual.  The FCSdata array contains the fitted parameters: Trace number, regression coefficient, number of molecules, pct of species A, Deq species1, Deq species2, flow, triplet pct, and triplet time constant.



* To determine the dependencies and necessary toolboxes for any Matlab script, run:
[fList, pList] = matlab.codetools.requiredFilesAndProducts('ScriptFileName.m') from the command line.  
fList is a cell array with the function dependencies, and pList is a structured array with the required Matlab products.
