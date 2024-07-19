# Random-Field-Generation
This matlab codes are used to generate 2-dimensional random fields using Local Average Subdivision (LAS) method.

If you are using MATLAB File Exchange site, you cannnot see the examples described below. You need to click onto my GitHub site where you can find examples.

`rfgen.m` is the MAIN file. It requires a single input of the input file name, having this format:
----------------------------------`example.dat`--------------------------------
1.  mean and standard deviation of random field #1  . . . 10.0 1.0 2 0 0 0 0
2.  mean and standard deviation of random field #2  . . . 20.0 2.0 2 0 0 0 0
3.  correlation lengths of random fields  . . . . . . . . 0.5 0.5
4.  cross-correlation between the two fields  . . . . . . 0.5
5.  element sizes in x- and y-direction . . . . . . . . . 0.01953125 0.01953125
6.  number of elements in x- and y-direction  . . . . . . 256 256
7.  initial seed number . . . . . . . . . . . . . . . . . 37
8.  total number of realizations  . . . . . . . . . . . . 100
9.  plot a random field?  . . . . . . . . . . . . . . . . 1 1 1
10. debug on? . . . . . . . . . . . . . . . . . . . . . . 1
----------------------------------`example.dat`--------------------------------

This is exactly the same as `example.dat`, which contains all the input parameters. `getinp.m` reads the the specific example input file, so it can be adapted for your own input file. `wrtinp.m` write results into `example.out` which contains all the simulation result statistics. You can adapt `rfgen.m`, `example.dat`, `getinp.m`, and `wrtinp.m` according to what you need.

In this example, 100 realizations are applied: `example_random field 1.fig` shows the generated random field 1 (the 1st realisation). `example_random field 2.fig` shows the generated random field 2 (the 1st realisation). The two random fields have a cross-correlation coefficient of 0.5 (see `example.dat`). `example_correlation structure_random field 1.fig` is the averaged auto correlation structure of all the 100 generated random field 1, compared to the prescribed theoretical structure. `example_correlation structure_random field 1.fig` is for random field 2.


Author:        Dr Pengpeng He    
Organisation:  University of Dundee    
Email:         phe001@dundee.ac.uk    
Website:       http://discovery.dundee.ac.uk/en/persons/pengpeng-he
