There are two versions available. One is for OpenMP shared memory parallelism useful for running on a single machine. The other version is
for MPI to use on a cluster utilizing multiple nodes. The input file in the MPI folder takes a very long time to run, so I recommend changing
it before attempting to run it. 

I use the clang compiler for the OpenMP version when on MacOS, and it seems to run considerably faster than the GCC counterpart. On Windows,
I use GCC. I compile with O3 optimization and for older versions of GCC, you will have to use the 2011 c++ standard flag, along with fopenmp.
On the cluster, I use the intel compiler. I have found it is much better with this code and I recommend using it if you can. For OpenMP, use
qopenmp and the 2011 C++ standard. For MPI, be sure to compile with impi compatibility if using Intel. A slurm batch file is included.
