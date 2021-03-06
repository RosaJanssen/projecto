# Protein Pow(d)er 
## Code 

### Algorithms 
Four different algorithms will be explained.

#### Random
* The random algoritm first asks the user for a protein structure and a number of times to run the algoritm. 
* It starts by folding the given protein in a random manner: the first aminoacid is placed in the middle of a grid. Then the algorithm checks whether the neighbours of the aminoacid (up, down, left and right) in the grid are free. The algorithm randomly chooses one of the free neighbours to place the succeeding aminoacid. The latter process is repeated until all aminoacids from the protein structure have been placed. If protein folding is valid, its score is being calculated by finding (HH and HC) bonds. The algoritm updates the best route and its score. The whole latter process is repeated as many times as indicated by the user. 
* The script visualizes the best protein route outcome in a grid and saves the characteristics of the best found protein route in a seperate file (resultsfile). The resultsfile can be used to create graphs.


#### Greedy
* The random algoritm first asks the user for a protein structure and a number of times to run the algoritm. 
* The algoritm starts by folding the given protein. The first aminoacid is placed in the middle of a grid. Then the algorithm checks whether the neighbours of the aminoacid (up, down, left and right) in the grid are free. In addition, it checks whether these free neighbours have H- or C-neighbours. If the latter is true, the algorithm randomly chooses from one or more of the neighbours that have H- Or C-neighbours. If this is not the case, the algorithm randomly chooses from the available neighbours. If all aminoacids can be placed and the route is valid, the score is being calculated by finding (HH and HC) bonds. The algoritm updates the best route and its score. The whole latter process is repeated as many times as indicated by the user. The whole latter process is repeated as many times as indicated by the user. 
* The random algoritm is split into two scripts. The first script (greedy_plotprotein.py) visualizes the best protein route outcome in a grid, while the second script (greedy_collectdata.py) saves the characteristics of the best found protein route in a seperate file (resultsfile). The resultsfile can be used to create graphs.

#### Breadthfirst 
* The breadthfirst algoritm first asks the user for a protein structure. 
* The algoritm creates a queue that lists all possible unique routes. Subsequently, the algoritm translates the routes into coordinates as to check whether the routes are actually valid (i.e. a route cannot pass the same spot twice). The algoritm calculates the score for only those routes that are valid, and from these returns the best protein with its attributes. Finally, the algoritm plots the best protein folding. 

#### Breadthfirst with look-ahead 
* The breadthfirst look-ahead algoritm first asks the user for a protein structure a number of times to run the algoritm.  
* The algorithm let's the user give a number for which the programm will optimize the next amino acids, let's call this number x for now. It does that by splitting the protein structure into different parts. Firstly it will take the coming x amino acids and decide which route will generate the best stability. Based on that the algorithm will fix the first move in the route, and adds it to the constant string. The algorithm repates this principle over the whole protein structure length that is given by the user. When there are more solutions with the same score in stability the algorithm makes a random choice between those.
* The breadthfirst look-ahead algoritm is split into two scripts. The first script (lookahead_plotprotein.py) visualizes the best protein route outcome in a grid, while the second script (lookahead_collectdata.py) saves the characteristics of the best found protein route in a seperate file (resultsfile). The resultsfile can be used to create graphs.

#### Histograms
The create_histogram.py code creates histograms based on the resultfiles made in some of the previous codes (_collectdata) and is meant to give insights in the results of algorithms. An output example can be found in the results folder. 

### Classes
In the classes folder, one can find the classes used in the algoritms.