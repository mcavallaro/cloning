
## Direct evaluation of the SCGF of the current of the semi-Markovian ASEP on a single site.

As imput we need to provide only the simulation time, the size of the ensemble, *s*,  and the destination folder:
```
one-node-ASEP-CTRW-clon.exe [# simulation time] [DIM. ENSEMBLE] [s] [folder]
```
The files `update.cpp` and `time_jumps.cpp` contain the definitions of functions that define the dynamics of the system.
In particular, we made use of the following parameters:
```
/*For the ``Gamma'' case:
double k01L = 1.5;
double k10L = 1.5;
double k01R = 0.5;
double k10R = 0.5;*/

/* Erlang case:*/
double k01L = 2;
double k10L = 2;
double k01R = 2;
double k10R = 2;

//In order to compare it to the Exponential case, we first fix the effective rates
double m01L = 0.1;
double m10L = 0.2;
double m01R = 0.3;
double m10R = 0.2;

//then we can derive the Scale parameters.
double alpha = m01L/k01L;
double beta  = m10L/k10L;
double gamm  = m01R/k01R;
double delta = m10R/k10R;
```

Use several values of *s*. The results are aggregated and plotted in a `.pdf` file by the script 
```
bash parse.sh [input_prefix] [output file without extension]
```
which also remove all the files with prefix `[input prefix]`.




