

# The onset of self-excited radial pulsation in evolved red giant stars

The Fortran program `intp.f` can be used to compute the critical luminosity corresponding to the onset of self-excited radial pulsation for an arbitrary combination of physical and chemical stellar parameters by interpolation in the grid of models published in [Trabucchi & Pastorelli, 2025](https://ui.adsabs.harvard.edu/abs/2025ApJ...978...30T/abstract) (their table 1).

## Installation

Download the source code `intp.f` and the data file `table1.dat`. The former can be compiled using your preferred Fortran-compatible compiler, for instance using `gcc`:

```
gfortran -o intp.exe intp.f
```

The compiled file has to be made executable. During execution, the program will try to access input unit 1 to get the grid data, so it is convenient to create a symbolic link to the file `table1.dat`. Assuming the latter is in the same directory as the executable `intp.exe`:

```
ln -s ./table1.dat fort.1
```

## Use

Enter the directory where the executable `intp.exe` is located then run like this:

```
./intp.exe <mode> <inputfile> <outputfile>
```

Here `<inputfile>` is the path to a file containing a list of combinations of the following parameters:

1. $X$: hydrogen abundance in the envelope by mass fraction
2. $Z$: envelope metallicity by mass fraction
3. $M$: total stellar mass in solar units
4. $\alpha_{\nu}$: adopted value of the turbulent viscosity parameter
5. $\log(T_{\rm eff})$: decimal logarithm of the effective temperature in ${\rm K}$

Each line must give each of these quantities, and nothing else, in this specific order, separated by spaces or tabs. Each line cannot be longer than 256 characters, and the file should have no headers nor footers. The results will be written to the output file `<outfile>` (make sure there is no such file, the program will stop if there is to avoid overwriting it). The output file will be essentially a copy of the input file, except 7 values are appended to each line:

1. $\log(L_{\rm c})$: the critical value of luminosity, for the specified combination of $X$, $Z$, $M$, $\alpha_{\nu}$, and $\log(T_{\rm eff})$, at which the envelope becomes unstable to self-excited radial pulsation
2. $\log(P_{\rm c})$: the period of the dominant mode the star pulsates in at the very onset of self-excited pulsation.
3. `fX`: linear interpolation flag for parameter $X$
4. `fZ`: linear interpolation flag for parameter $Z$
5. `fM`: linear interpolation flag for parameter $M$
6. `fAN`: linear interpolation flag for parameter $\alpha_{\nu}$
7. `fT`: quadratic interpolation flag for parameter $\log(T_{\rm eff})$

**WARNING**! The quantity $L_{\rm c}$ provides a criterion to assess whether or not a red giant star is unstable to radial pulsation: it is if its luminosity is $L\geq L_{\rm c}$, otherwise it is not. The quantity $P_{\rm c}$ is the period the star pulsates in at the very onset of pulsation. $P_{\rm c}$ is **NOT** the pulsation period for the input combination of parameters, except in special case in which the input luminosity is such that $L=L_{\rm c}$.

The interpolation flags allow to keep track of whether, instead interpolating along a given dimension (the corresponding flag has value 0), extrapolation below or above the grid bounds occurred (the corresponding flag has values 1 or 2, respectively). Extrapolation is allowed as it is not expected to result in bad behavior as long as it does not happen too far from the grid boundaries. Yet such cases should be treated with care.

The option `<mode>` must be either `s` (as in "static") or `d` (as in "dynamic"). It is used to specify whether, for the purpose of interpolation, you want to adopt the values of luminosity and effective temperature taken from an hydrostatic envelope model or derived from a time average of the corresponding hydrodynamic time series. The former is recommended, especially when one aims to compare with results from "standard" stellar evolution codes that enforce hydrostatic equilibrium at any time.

## Examples

Example input/output files are provided with the files `test.in` and `test.out`. They can be used as references or templates for building your own input files.

## Notes

The program interpolates the quantity $\log(L_{\rm c}/{\rm L}_{\odot})$ and the quantity $\log(P_{\rm c}/{\rm d})$ quadratically in $\log(T_{\rm eff}/{\rm K})$ and linearly in $X$, $Z$, $M/{\rm M}_{\odot}$, and $\alpha_{\nu}$.

Note that $X$, $Z$, and $M$ are the **current** values, not the initial values! In other words, if you use this program for a stellar evolutionary track, you need to give it the values at a given time step, not the ones representing the zero-age main sequence.

The choice of the value to adopt for the turbulent viscosity parameter $\alpha_{\nu}$ is a bit critical, as the parameter is uncalibrated. Users without specific requirements should simply use $\alpha_{\nu}=0$, which corresponds to assuming no turbulent viscosity at all. In this case, they should keep in mind that the resulting value of $L_{\rm c}$ represents the "earliest" possible onset of self-excited radial pulsation for the chosen combination of parameters, i.e., it is formally a lower limit.

## Disclaimer

This software and the associate data are provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement. In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the software or the use or other dealings in the software, as well as the associated data.
