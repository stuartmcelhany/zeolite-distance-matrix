# Zeolite Distance Matrix
## Information
This python script takes a zeolitic unit cell and outputs a distance matrix corresponding to the number of tetrahedra links between every pair of aluminums.

The input file must be in .xyz (cartesian coordinate) format. Additionally, the current format only supports rhombohedral lattice systems, meaning it only takes in one angle and one side length.

The file within this repository titled "isomer-0.xyz" can be used as an example. It is one of the many isomers of Faujasite-5.

# How it works
The script works by first making a 3x3x3 lattice of the initial unit cell that removes all atoms except Al and Si (these are what make up distinct linkages in zeolites like faujasite).
