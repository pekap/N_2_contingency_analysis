## Fast N-2 contingency analysis algorithm

Basic implementation of fast N-2 contingecy analysis algorithm developed by Prof. Konstya Turitsyn MIT research group
The algorithm aims to find all dangerous N-2 contingencies in the grid. 
We call dangerous N-k contingency a set of k elements which tripping leads to violation of at least one constraint in 
power grid. 

In this implementation we consider only line as contingency elements and their power limits as constraints.

### How to run
Download this repo and add path to it to matlab path collection. 
You can see example of usage in `code/tests/test_grid_class.m` file. In that file we

1. Load a MATPOWER grid case from `/code/cases/` folder
  
  ```
  runcase = loadcase('case300');
  ```
2. Initialize grid object using constructor of the Grid_class
  
  ```
  grid = Grid_class(runcase,'case300_experiment');
  ```

3. Run N-0, N-1 and N-2 analysis. 

  ```
  grid.N_1_analysis();
  grid.N_2_analysis('fast');
  grid.N_2_analysis('bruteforce');
  ```

After running this, inside MATLAB command window we should see

```
OPF was successfuly solved
***********************************************
************   Start N-1 analysis  ************
***********************************************
	Processed 322/322. Number of dangerous N-1 contingecies is 0 : 

	N-1 analysis was performed,0 dangerous contigencies were found, 0 lines are violated 

***********************************************
************Start fast N-2 analysis************
***********************************************
	 0 iteration: number of potential confingecies::51558, B::103362
	 1 iteration: number of potential confingecies::4, B::291
	 2 iteration: number of potential confingecies::1, B::6

	Bruteforce enumeration over 1 pairs 
	Processed 100 percent. Number of contingencies 1; fake 0
	Running time for fast algorithm is 7.487977e-02 sec 

***********************************************
************Start bruteforce N-2 analysis************
***********************************************

	Bruteforce enumeration over 51842 pairs 
	Processed 100 percent. Number of contingencies 1; fake 123
	Running time for brute force algorithm is 1.420059e+01 sec 
```

In this case running time of our algorithm was almost 200 _~O(Number of branches)_ times less direct bruteforce enumeration.

### What happens behind scenes

1. We create grid object. During the process of the creation we
  1. Remap odd numeration. (some matpower cases have buses numerated in odd order)
  2. Delete branches that have 'off' status.
  3. Set shift angle to 0 and tap coefficient to 1 (generalization of our algorithm without this is trivial, but for demonstration we are left with this normalization)
  4. Remove parallel lines
  5. Run DC OPF
  6. Aggregate leafes (they lead to trivial N-1 contingencies which we are not interested in)

2. We run N-1 contingency analysis
  1. Switch every single line 'off' one by one
  2. Look at response on other lines
  3. Create matrix of LODF
  4. Create array of max margins, to understand if grid is N-1 safe
  5. Save results to file in '/results/' folder

3. We run N-2 contingency analysis
  1. It checks if N-1 analysis have already been run.
  2. If not runs it / else loads N-1 analysis data
  3. Starts appropriate analysis method
    1. 'bruteforce' just enumerates all possible pairs
    2. 'fast' uses developed algorithm to shrink the search space and then enumerates over it
  4. Results are being recorded to file in '/results/' folder. The results file contains all dangerous contingencies 

For further explanation see [description of the code](https://github.com/pekap/N_2_contingency_analysis/tree/master/code) and comments inside class files

### Links

We use dcopf solver and grid cases from [MATPOWER](http://www.pserc.cornell.edu/matpower/) MATLAB package, which is included in our repo.
  
The algorithm presented here was developed in Prof. Kostya Turitsyn research group at MIT. Links to related works:  
[Original paper](http://arxiv.org/pdf/1211.0728.pdf)  
[Elaborated approach for generators as contingencies](http://arxiv.org/pdf/1303.3938.pdf)

### For any suggestions please contact us at

Petya Kaplunovich pekap@mit.edu  
Kostya Turitsyn   turitsyn@mit.edu
