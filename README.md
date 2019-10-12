# LCM-game

Python2 code to perform a simple structural simulation of TM plane in L-iCo-Mn oxide and to analyze the percolation of the resulting structure. This is the code I used when working on https://doi.org/10.1021/cm200250a

In a nutshell, the problem was apparently contradictory experimental data. EXAFS showed clear Co-Mn segregation, while HAADF microscopy showed local (nanometer scales) structure variations consistent with segregation. However, EELS showed that the Co/Mn ratio was the same in both kinds of regions. The idea was to perform a structure simulation imposing the known Co-Mn simulation on the EXAFS scale (a cluster of a few atoms), on models consistent with microscopy scales (tens to thousands of atoms wide.

First, a hexagonal supercell, size atoms wide is created. Then the supercell is filled by 25 atom big Co or Mn tiles, consistent with the bulk coordination of LiCoO2 and Li2MnO3. At each simulation step, a tile is placed in a random position within the supercell. The tile composition, Co or Mn, is selected with probability given by the intended overall model Co/Mn ratio; except when the tile touches an existing domain, in which case the domain is extended. This approach ignores energy of interaction between domains, therefore is purelly driven by enthropy.

![model step](/img/figure__II_step_II.bmp)

When the supercell is full, the game counts the number of disconnetced Co and Mn domains in the model, applying periodic boundary conditions. The idea is that two not connected Mn domains have only a 1/3 probability of having their Li columns alligned; if they are not, the HAADF projection averages out to look like Co.

![model full](/img/figure__II_CoMn.bmp)

By running the game repeatedly for different sizes and Co/Mn ratios, the percolation threshold can be estimated. \src\game.py contains the routines to run the game and count disconnect domains. \src\analysis.py contains code to analyze the model outputs and reformat the data for further plotting. The ```if __name__ == '__main__':``` section at the end of each file contains examples of how I run it for testing and experiemnting. Different runs are commented out. The outputs of those runs are stored under \old_runs. Feel free to explore and have fun!

![model stats I](/img/Model_clusters_1col.tif)

![model coordination](/img/Model_coordination_1col.tif)
