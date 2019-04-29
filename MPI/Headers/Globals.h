extern const int threads; //number of processor threads to compute with
extern const double m;     //mass eV/c^2
extern const double PI;    //pi
extern const double h2m;   //hbar^2/m eV-A^2
extern const double ee;    //electron charge in eV-A
extern const double aBohr; //Bohr Radius in A
extern const double hbar;  //hbar eV-s
extern const double ch;    //electron charge C
extern const double SOL;   //speed of light in m/s
extern const double hbar2m;//hbar^2/2m

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <complex>
#include <string>
#include <cstring>
#include <string.h>
#include <cstdlib>
#include <chrono>
#include <ctime>

#include "Besselj.h"
#include "Bessely.h"
#include "LegendreP.h"
#include "Numerov.h"
#include "Potentials.h"
#include "Derivatives.h"
#include "Periodic.h"
#include "mpi.h"