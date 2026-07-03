# Macro for Event Analysis

This directory is to analyze the events generated under `macro_gen/`.


## Creation of Histograms from DST Files

You first analyze the DST files of each event set to extract histogmras.
You can use `ana-one-set.sh` as follows:

```
./ana-one-set.sh jpsi
./ana-one-set.sh psip
./ana-one-set.sh dy
```

Outputs of each event set will appear in `jpsi/`, `psip/` or `dy/`.
You might comment out some lines in `ana-one-set.sh`,
because for example you need not execute `make-dst-list.sh` when the input DST files haven't changed.


## Computation of Anticipated Statistics

`DrawStatExp.C` is to draw the expected accuracy of J/psi TSSA.

```
root -b DawStatExp.C
```


## Drawing the mass distributions

This step is optional.
`DrawMassDist.C` is to draw the mass distributions of J/psi, psi', D-Y and their totals.

```
root -b DawMassDist.C
```
