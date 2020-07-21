# Zeolite Distance Matrix
## Information
This python script takes a zeolitic unit cell and outputs a distance matrix corresponding to the number of tetrahedra links between every pair of aluminums.

The input file must be in .xyz (cartesian coordinate) format. Additionally, the current format only supports rhombohedral lattice systems, meaning it only takes in one angle and one side length.

The file within this repository titled "isomer-0.xyz" can be used as an example. It is one of the many isomers of Faujasite-5.

## How it works
The script works by first making a 3x3x3 lattice of the initial unit cell that removes all atoms except Al and Si (these are what make up distinct linkages in zeolites like faujasite). The goal is to find the shortest path between the aluminums in the original unit cell to every other aluminum in the original unit cell. However, the original unit cell doesn't accomadate for periodic conditions; hence, the need for the 3x3x3 lattice to check periodic conditions. The script looks for an aluminum atom found in the original unit cell, then it utilizes the Dijkstra algorithm to find the shortest path to another aluminum atom. From periodic conditions, there exists 26 other possibilities for the location of the second aluminum, so the Dirjkstra algorithm is ran for each of the next 26 total possibilites. At this point all 27 possibilites (the original aluminum and following 26 possibilites) are checked and the smallest path is returned. The script repeats this process for all possible pairs of aluminums and then outputs a .csv file to view the distance matrix.
