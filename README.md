# Random-Field-Generation
This matlab codes are used to generate 2-dimensional random fields using Local Average Subdivision (LAS) method.

There is an example: `example.dat` contains all the input parameters. `getinp.m` reads the the specific example input file, so it can be adapted for your own input file. `wrtinp.m` write results into `example.out` which contains all the simulation result statistics.

`rfgen.m` is the MAIN file.

In this example, 100 realizations are applied: `example_random field 1.fig` shows the generated random field 1 (the 1st realisation). `example_random field 2.fig` shows the generated random field 2 (the 1st realisation). The two random fields have a cross-correlation coefficient of 0.5 (see `example.dat`). `example_correlation structure_random field 1.fig` is the averaged auto correlation structure of all the 100 generated random field 1, compared to the prescribed theoretical structure. `example_correlation structure_random field 1.fig` is for random field 2.


Author:        Dr Pengpeng He    
Organisation:  University of Dundee    
Email:         phe001@dundee.ac.uk    
Website:       http://discovery.dundee.ac.uk/en/persons/pengpeng-he
