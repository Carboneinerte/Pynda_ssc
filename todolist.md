# Spatial Analysis
- [X] Cell neighbors (AD project, squidpy)
- [ ] Cell to cell communication (L-R) : ongoing (AG)
    - [ ] Spacia was tested but long run and suboptimal result
    - [ ] CellChat tested, promising, allow condition comparisons


# Metacycle and CircaComapre
- [X] Implement Circacompare, compare with metacycle
- [ ] Implement possibility to continue an interrupted analysis (by scanning which analysis file already exist and continue the analysis once a file is missing)
- [ ] Add possiblity to run the summary files analysis separately
- [ ] Optimize memory (RAM) usage: check if large object/variable can be deleted if no more used in the analysis

- [ ] Include circacompare analysis in notebook
- [ ] Save logs of output for each run as separate file. ~~R function "sink()" should do the trick. Save in a subfolder "log"~~ "Logr" package is way better.


# Misc.
- [X] Same scale polar plot
