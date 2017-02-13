# Presentation

Code of the "BCM" model, used in 

- Gallée H., J. P. Ypersele, T. Fichefet, C. Tricot and A. Berger (1991), Simulation of the last glacial cycle by a coupled, sectorially averaged climate-ice sheet model. Part I : The Climate Model., Journal of Geophysical Research, (96) 13139-13161 doi:10.1029/91JD00874

- Gallée H., J. P. Ypersele, T. Fichefet, I. Marsiat, C. Tricot and A. Berger (1992), Simulation of the last glacial cycle by a coupled, sectorially averaged climate-ice sheet model. Part II : Response to insolation and CO$_2$ variation, Journal of Geophysical Research, (97) 15713-15740 doi:10.1029/92JD01256

Code legacy by authors (c) 1991 , with some minor reorganisation + introduction of a stochastic parameterisation by M. Crucifix (2008 - 2016) 

Acknowledge the authors of the original publications publications

# Running

```
./configure
make
qsub
```

inspect qsub and see what it does ! 

# Directory structure and files

## Inititialisation files  [init]

- cal5,6,7.start : initial ice sheet profiles with 0.5d resolution
- cal.start : net accumulation - ablation balance over ice sheets with 5 degree resolution


## Configuration files  [etc]


- bcmin contains (line 4) the start and the end of the experiment, in ka 

## External data files [data] 

- contains among others the co2 data used by the model. Do not change without inspecting the code and make sur you understand how it is read. As it is the fact that co2 data contains data over the last 800 ka is hardwired. See function  `palco2 (apal)` in bcm.f 

- the bcm.. files  contain spectral data for infrared and solar radiation. 

# source in [src] 

(src/bcm.f)[src/bcm.f] is the atmosphere model
(src/bcm.f)[src/cal.f] is the ice sheet model, runnig 1000 years (or see `calctr` to change this)

Note: (src/bcm.f)[bcm.f] includes the function `co2inter` for a possible inclusion of interactive CO2 as use in the master thesis of Maxime Antoine. Not used at the moment.