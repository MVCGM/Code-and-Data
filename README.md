# Code-and-Data
Code and data to reproduce NeurIPS submission results.  All code is written in Matlab.

# Simulation results
CGM_test() will generate a new simulated data set with two views, and run both 2 CGM and Multiview CGM algorithms on the data.

Usage: CGM_test(n,t,iters,method) 
n - size of location grid, computed as n x n
t - number of time steps
iters - maximum number of iterations for each algorithm
method - Method for simulating data.  Must be values 0, 1, or 2.
  0: Attract-Same method
  1: Attract-Diff method
  2: Attract/Repel method
  
Output: vector mets
  mets[1]: Total difference of incoming and outgoing flow for predicted values of 2 CGM algorithm.
  mets[2]: Total difference of incoming and outgoing flow for predicted values of Multiview CGM algorithm.
  mets[3]: NAE view 1 for 2 CGM algorithm.
  mets[4]: NAE view 2 for 2 CGM algorithm.
  mets[5]: NAE of both views combined for 2 CGM algorithm.
  mets[6]: NAE view 1 for Multiview CGM algorithm.
  mets[7]: NAE view 2 for Multiview CGM algorithm.
  mets[8]: NAE of both views combined for Multiview CGM algorithm.
  mets[9]: Runtime for 2 CGM algorithm.
  mets[10]: Runtime for Multiview CGM algorithm.
  
script(reps,n,t,iters,method) will make 'reps' number of calls to CGM_test with the given parameters.  
script() returns a vector m that is the average value of 'mets' for each call, and a matrix M that contains all of the individual 'mets' vectors from each call.

To reproduce paper simulation experiments (with different randomized datasets):
[m5_0,M5_0] = script(10,5,100,500,0)
[m5_1,M5_1] = script(10,5,100,500,1)
[m5_2,M5_2] = script(10,5,100,500,2)
[m10_0,M10_0] = script(10,10,100,500,0)
[m10_1,M10_1] = script(10,10,100,500,1)
[m10_2,M10_2] = script(10,10,100,500,2)

# Japan population flow data
run_Japan_data() will load the pre-processed Nightly Inc data, and run the Single CGM, 2 CGMs, and Multiview CGM algorithms on it.  Uses parallel processing.

Usage: run_Japan_data(pref,Fout,max_iter)
  pref - Determines which city data to use.  Must be the strings 'Tokyo', 'Osaka', or 'Nagoya'.
  Fout - Base string for output files.  Should end in '.m'.
  max_iter - Maximum number of iterations for each algorithm.
  
Output: The script will save 3 output files to the current directory.  Each file will end in the string specified by 'Fout', and start with the strings 'base_', 'base_all_', and 'multi_'.  These correspond to the 2 CGMs, Single CGM, and Multiview CGM algorithms, respectively.
Each of these output files will have the predicted values of M and Z, indexed by view, the NAE values of each view, and a runtime t.  The Single CGM algorithm will only have outputs for a single view, which is the combination of both views.

To reproduce paper experiments:
run_Japan_data('Tokyo','Tokyo.m',500)
run_Japan_data('Osaka','Osaka.m',500)
run_Japan_data('Nagoya','Nagoya.m',500)

# US job and migration data
run_US_data() will load the pre-processed census data and run the 2 CGMs and Multiview CGM algorithms on it.

Usage: run_US_data(Fout,max_iter)
  Fout - Base string for output files.  Should end in '.m'.
  max_iter - Maximum number of iterations for each algorithm.
  
Output: The script will save 2 output files to the current directory.  Each file will end in the string specified by 'Fout', and start with the strings 'base_' and 'multi_'.  These correspond to the 2 CGMs and Multiview CGM algorithms, respectively.
Each of these output files will have the predicted values of M and Z, indexed by view, the NAE values of each view, and a runtime t.  View 1 is the flow of jobs between states, view 2 is the migration of people between states.

To reproduce paper experiments:
run_US_data('US.m',500)
