# Thomas_Fermi_Scattering
Thomas-Fermi Scattering Cross Section Numerical Solver

This program is being written to practice better C++ programming and for my graduate
level computational physics course and is not meant to be a stable release of some
simulation software that is intended for end-users.

Required files are:
`ThomasFermi.cpp`
`Globals.h`
`Besselj.h`
`Bessely.h`
`Legendre.h`
`Potentials.h`
`Numerov.h`
`Derivatives.h`
`Periodic.h`


On MacOS Mojave with clang++ and libomp, I compile with:

`clang++ -Xpreprocessor -fopenmp -lomp -o Scattering ThomasFermi.cpp`
