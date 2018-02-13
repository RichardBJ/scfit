Single channel mechanism fitting
================================

12 February 2018


Description
--------------------------------

This application takes an idealized list of openings and closings from single channel recordings and fits rates in a single channel gating mechanism using an exact correction for missed events. The application is based largely on the work of Hawkes, Jalali, and Colquhoun (1990), Colquhoun and Hawkes (1987) and Hawkes and Jalali (1992).


Getting Started
--------------------------------

A standalone GUI application can be [downloaded](https://github.com/ogdenkev/scfit/releases/download/v0.1.0-alpha/SCFit_Installer.zip) and installed by running SCFitInstaller_web.exe.

To use the source code, add all the folders to the Matlab path. This can be done in the Current Folder window, by selecting all the folders, right clicking and then selecting Add to Path > Selected Folders and Subfolders. You can also enter `addpath DIRECTORY` from the Matlab command line. For example, if the path to the `SC_ML_fit directory` is `C:\MATLAB\SC_ML_fit` then 

```
addpath C:\MATLAB\SC_ML_fit\src C:\MATLAB\SC_ML_fit\SCAN C:\MATLAB\SC_ML_fit\Qmatrix C:\MATLAB\SC_ML_fit\exp_mix_dist
```

will add the appropriate directories. 


Building MEX files -- compiled Matlab code
------------------------------------------

The application consists of a number of Matlab functions and will run after adding the `src`, `SCAN`, `Qmatrix`, and `exp_mix_dist` folders to the Matlab path. However, considerable increases in the speed of fitting can be achieved by compiling a few of the bottleneck functions.

To compile these functions, add the `build` directory to the Matlab path and then run the `generate_Qmatrix_C_code_v2.m` script. This will compile three functions: Qmatrixloglik, `Qloglik_bursts`, and hjcdist. It will also generate a `codegen` folder that contains all the code needed to compile these functions, but all that's needed are the three MEX files.


Quick start
------------------------------------------

This single channel fitting mechanism application comes with several demo datasets and mechanisms. Here we will fit the channel mechanism in Chapter 19 of Single Channel Recording.

1. Add the scripts directory to the Matlab path
    `addpath scripts`
2. Run the HJCFIT_demo.m script
    `HJCFIT_demo`

This demo will simulate openings and closings from a five state mechanism with two open states and three shut states. After simmulating the openings and closigns, rates from the mechanism will be fit.  A window will pop up showing the open and closed time histograms and a line showing the predicted pdfs from the rates of the current iteration. A number of metrics will also be printed on the Matlab command line, including the opposite of the log-likelihood. As Matlab optimization functions only minimize functions, so to get the maximum log-likelihood, we minimize the opposite.

The fitting should complete in less than a minute, and the results will be saved to the Matlab workspace. Of note are the variables `rates`, which contains the fitted rates, `ll`, which contains the log-likelihood, and `qnew`, which contains the fitted Q matrix. After fitting completes, the open and shut state histograms will be displayed with the fitted pdfs and the predicted pdf if no evenets were missed.

### Fitting your own data

Try running the script `runFitting_QuB_model.m` or `runFitting_Excel_model.m` for mechanisms in QuB qmf format and Excel format, respectively.


Loading in idealized openings and closings
------------------------------------------

Idealized openings and closings can be imported from QuB DWT files and from SCAN SCN files.

### QuB DWT Import

```
[dwells, states, ndwells, amps] = dwtread;
```

### SCAN SCN Import

```
[durations, amps, ~, cal] = scanread("scansamp.scn");
```


Specifying a mechanism
----------------------

Single channel gating mechanisms consist of the number of states in the mechanism, the connectivity of those states, which states are open (i.e.conducting), the rates of transition between states, and which rates are contrained, for example by microscopic reversibility.

Mechanisms can be specified manually in Matlab code, by importing a QuB qmf file, by importing a model in Excel, or by importing a ChanneLab MDL file.

### QuB qmf model

```
[q, A, F, idxall, idxvary, gamma, xi, numstates, ~, ~, ~, ~, filename] = qmfread(qmffile);
q = q*1e-3;
```

### Excel model

See the test file `qtest.xlsx` for the format of Excel model files

```
[q, A, F, idxall, idxvary, gamma, xi, fname] = qmatxlsread("qtest.xlsx")
```


Imposing open and shut time resolutions
---------------------------------------

```
%% Impose the resolution
open_deadtime = 0.04; % milliseconds
shut_deadtime = 0.04;
[resolved_dwells,resolved_states] = imposeres (durations, amp, open_deadtime, shut_deadtime);
```


