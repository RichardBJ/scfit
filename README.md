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

An ion channel gating mechanism can be specified in an Excel workbook (.xlsx). The workbook must contain at least three sheets, named `States`, `Rates`, and `Constraints`. All data must start in the first column, i.e. cell A1, in each tab.

#### States Worksheet

| State | Conductance | Name  | Notes |
|-------|-------------|-------|-------|
| 1     | 0           | R     | shut  |
| 2     | 0           | RA    | shut  |
| 3     | 0           | RA2   | shut  |
| 4     | 0           | RA2f  | shut  |
| 5     | 0           | RA2s  | shut  |
| 6     | 1           | RA2fs | open  |
| 7     | 0           | D1    | shut  |
| 8     | 0           | D2    | shut  |

The States tab specifies how many states the model will have. It must contain at least two columns named `State` and `Conductance`. An optional third column named `Name` may also be present.  Any other columns are ignored, but may be useful for keeping notes on the states in the model.

The State column contains the number of the state in the model. State numbers must start at 1 (for the first state) and increase by 1 for each additional state. Thus, this column should simply be number from 1 to `N` where `N` is the number of states in the model.

The Conductance column should contain the conductance of the state. Any conductance equal to zero (0) will be considered a shut state and any state with a conductance not equal to zero will be considered open.

The Name column, if present, can be a name for the state.  It is best not to leave any state unnamed. Rather, simply omit the Name column, or create a filler name, like `State3`.

#### Rates Worksheet

| State1 | State2 | Value  | Note   |
|--------|--------|--------|--------|
| 1      | 2      | 0.0632 | 2*kon  |
| 2      | 1      | 1.01   | koff   |
| 2      | 3      | 0.0316 | kon    |
| 3      | 2      | 2.02   | 2*koff |
| 3      | 4      | 3.14   | kfp    |

The Rates tab specifies the connectivity between states in the model. It must contain at least three columns named `State1`, `State2`, and `Value`. It can also contain other columns, but these will be ignored. Hopefully this is straightforward, but every rate in the mechanism corresponds to a transition from one state to another. Hence, `State1` is the number of the first state in the transition and `State2` is the second state. `Value` is the rate constant for that transition. 

#### Constraints

| State1 | State2 | Type      | SourceRate | Value |
|--------|--------|-----------|------------|-------|
| 2      | 3      | constrain | 1,2        | 0.5   |
| 3      | 2      | constrain | 2,1        | 2     |
| 4      | 6      | constrain | 3,4        | 1     |
| 6      | 4      | constrain | 3,5        | 1     |
| 5      | 6      | constrain | 4,3        | 1     |
| 6      | 5      | fix       |            | 0.05  |

The Constraints tab specifies any constraints in the model. Currently only 2 types of constraints are accepted in the Constraints table: `fix` and `constrain`. These should be sufficient to specify any physical constraint on the model (but what about [microscopic reversibility](#microscopic-reversibility)?).

The Constraints tab must have 5 columns named `State1`, `State2`, `Type`, `SourceRate`, and `Value`. Other columns may be present, but they will be ignored. As in the Rates worksheet, `State1` and `State2` specify the state numbers for the starting and ending states of the transition whose rate should be constrained. The `Type` for each constraint must be either `fix` or `constrain`. `SourceRate` is used only for the `constrain` type (see [below](#constrain-rates)).

##### Fix Rates

The `fix` type of constraint will simply set the rate from `State1` to `State2` to `Value`. The `SourceRate` column is ignored for this type of constraint.

##### Constrain Rates

The `constrain` type of constraint will set the rate from `State1` to `State2` to some multiple of `SourceRate`. The multiplier is given by the `Value` column. `SourceRate` must be a comma-separated string of the form `N1,N2` where `N1` is the number of the first state in the source transition rate and `N2` is the number of the second state in the source transition.

Imposing open and shut time resolutions
---------------------------------------

```
%% Impose the resolution
open_deadtime = 0.04; % milliseconds
shut_deadtime = 0.04;
[resolved_dwells,resolved_states] = imposeres (durations, amp, open_deadtime, shut_deadtime);
```

### Microscopic Reversibility

Microscopic reversibility is enforced on all single channel gating mechanisms
using the minimum spanning tree method of Colquhoun et al. (2004). Each mechnism
is treated as an undirected graph with the states equivalent to
nodes and the transitions between states equivalent to edges in the graph.
Then a set of rates is found that will be constrained to enforce microscopic
reversibility.

Essentially, a minimum spanning tree is found in the graph that represents the gating model.
Rates in the mechanism with physical or theoretical contraints are assigned
a weight in the undirected graph to ensure they would be included in the minimum
spanning tree, if possible.
Once the minimum spanning tree is found, graph edges not part of the MST are
selected to be constrained.  Edges correspond to transitions between states (i.e. rates) in
the gating mechanism, and because each transition is reversible there are two
rates associated with each edge. Therefore, given an edge not in the MST, there are
two choices for which rate to constrain, and one of the rates is arbitrarily
selected.

For each rate selected to be constrained so that the mechanism obeys microscopic
reversibility, one of the cycles in the mechanism that contains the rate is
found by determining the shortest path between the pair of states connected by
the rate's corresponding transition. This cycle was used to set the contraint.

After constraining rates to enfoce microscopic reversibilty, physical or theoretical constraints are added.

### References

1. Colquhoun D, Hatton CJ, Hawkes AG (2003) The quality of maximum likelihood estimates of ion channel rate constants. J Physiol 547:699–728 
1. Colquhoun D, Hawkes AG, Srodzinski K (1996) Joint Distributions of Apparent Open and Shut Times of Single-Ion Channels and Maximum Likelihood Fitting of Mechanisms. Philosophical Transactions: Mathematical, Physical and Engineering Sciences 354:2555–2590 
1. Hawkes AG, Jalali A, Colquhoun D (1990) The Distributions of the Apparent Open Times and Shut Times in a Single Channel Record when Brief Events Cannot Be Detected. Philosophical Transactions: Physical Sciences and Engineering 332:511–538 
1. Hawkes AG, Jalali A, Colquhoun D (1992) Asymptotic Distributions of Apparent Open Times and Shut Times in a Single Channel Record Allowing for the Omission of Brief Events. Philosophical Transactions: Biological Sciences 337:383–404 
1. Colquhoun D, Dowsland KA, Beato M, Plested AJR (2004) How to Impose Microscopic Reversibility in Complex Reaction Mechanisms. Biophysical Journal 86:3510–3518
1. Jalali A, Hawkes AG (1992) Generalised Eigenproblems Arising in Aggregated Markov Processes Allowing for Time Interval Omission. Advances in Applied Probability 24:302–321 
1. Qin F, Auerbach A, Sachs F (1996) Estimating single-channel kinetic parameters from idealized patch-clamp data containing missed events. Biophys J 70:264–280 
1. Golub, G. H., and C. F. Van Loan (1989) Matrix Computations. Johns Hopkins University Press, Baltimore.
