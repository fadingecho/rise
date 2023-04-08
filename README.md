# Influence Maximization based on Reverse Influence Sampling and Evolutionary Algorithm (RISE)

## Building

```
make
```

## Running
Here is an example. Of course you can try other parameters.

`./moga_demo -i ./CA-GrQc.txt.encoded.bin -psize 150 -maxgen 500 -k 10 -m IC -crate 0.75 -rpratio 0.1 -numrf 30000 -amode 10000 -rrate 0.75`

## Input
Edges lists format supported, see [SSA Instructions](./SSA/Instructions.pdf). For convinence, we provide CA-GrQc binary file.

## Output
Ignore the first line, it is deprecated now. And the second line goes the parameter settings : ls_ratio, maxage, size_of_hypergraph, popsize, maxgen, ls_rate. Following lines show the results, the last percentage number is the the proportion of the most individuals.

## Verification
You can use the tool in SSA folder, see [SSA Instructions](./SSA/Instructions.pdf). If you do that, do not use multi-thread under IC model, there seems to be a bug. In the paper, we use monto-carlo simulation results.

## Other Information
I built this project based on https://github.com/KRM7/genetic-algorithms and
  https://github.com/hungnt55/Stop-and-Stare. And I am still working on this project.
