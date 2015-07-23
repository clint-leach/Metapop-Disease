Metadata for generated data
==========================

## pathgrid.RData

R environment containing output of `pathgrid.shell.r` for fully connected metapopulation.  Model output is stored in the `out` object, a 60000 by 12 dataframe  with columns:

* `treatment`: either `+high` or `+low` indicating either high or low quality habitat distribution, respectively
* `longevity`: the pathogen's environmental longevity (defined as the pathogen's half-life and measured in units of the expected occupancy time of hosts on healthy unit quality patch)
* `delta`: direct transmission probability
* `nu`: infectious survival
* `S`: final susceptible occupancy
* `I`: final infectious occupancy
* `R`:  final reservoir/empty occupancy
* `S.pop`: mean total susceptible population size over last 500 time steps
* `I.pop`: mean total infectious population size over last 500 time steps
* `maxI`: maximum infectious occupancy observed
* `tmaxI`: time to maximum infectious occupancy
* `quality0`: quality of randomly selected initial infectious patch

## pathgrid(lattice).RData

Same as above, but for metapopulation with patches arranged in a lattice (on a torus), so that each patch is only accessibly from its four nearest neighbors.

## trap.RData

R environment containing output of `trap.shell.r` for fully connected metapopulation with fixed direct transmission probability (`delta = 0.3`) and infectious survival (`nu = 0.2`).  Model output is stored in the `out` object, a 400000 by 19 data frame containing patch-level simulation results with columns:

* `quality`: quality of individual patch
* `tinf`: time until patch first infected
* `I`: the proportion of time patch spent occupied by infectious hosts
* `S`: the proportion of time patch spent occupied by susceptibly hosts
* `E`: proporion of time patch spent unoccupied/home/clint/Downloads/prsb.bst/home/clint/Downloads/prsb.bst
* `inf.events`: total number of infection events on the patch (S to I transitions)
* `susc.col`: total number of susceptible colonizations (E to S)
* `inf.col`: total number of infectious colonizations (E to I)
* `susc.ex`: total number of susceptible extinctions (S to E)
* `inf.ex`: total number of infectious extinctions (I to E)
* `Sfin`: final susceptible occupancy in metapopulation
* `Ifin`: final infectious occupancy in metapopulation
* `Efin`: final empty occupancy in metapopulation
* `maxI`: maximum infectious occupancy reached in metapopulation
* `quality0`: quality of randomly selected initial infectious patch
* `repID`: unique replicate ID
* `longevity`: the pathogen's environmental longevity (defined as the pathogen's half-life and measured in units of the expected occupancy time of hosts on healthy unit quality patch)
* `xi_im`: value of parameter controlling how immigration rate scales with patch quality/population size
* `xi_em`: value of parameter controlling how emigration rate scales with patch quality/population size
