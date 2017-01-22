# Floyd–Warshall: All-Pairs Shortest Paths

<br>
Implementation of the [Floyd–Warshall](https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm) algorithm using [Open MPI](https://www.open-mpi.org/) to develop a parallel programming for distributed memory environments.

The problem to solve consists in finding the shortest path between all pairs of nodes in a given graph ***G = (V, E)*** where ***V*** represent the vertexes and ***E*** represents the edges. The idea is for every pair of vertexes ***(Vi, Vj)*** find the size of the shortest path possible that connects ***Vi*** to ***Vj***.

<br>
#### The input file must have the following format:

* The first line must have the parameter **N** which is the number of rows and columns of the matrix that represents the graph
* The **N** following lines must have also **N** columns with numbers where each one represents the distance between the vertex of the line ***i*** and the vertex of the column ***j***


<br><br>
## Dependencies
- [Open MPI](https://www.open-mpi.org/software/ompi/v2.0/)


<br><br>
## How to Build
```bash
$ make build
```


<br><br/>
## How to Run
```bash
$ mpirun -np ${N_PROCESSES} -hostfile ${HOSTFILE} floyd < ${INPUT_FILE}
$ # Example: mpirun -np 4 -hostfile hosts.txt floyd < inputs/input300x300.txt
```

<br><br>
## Run Tests
```bash
$ make tests
```


<br><br>
## Clean Build Files
```bash
$ make clean
```
