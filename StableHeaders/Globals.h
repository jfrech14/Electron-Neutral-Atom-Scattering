extern int MeshSize; //mesh size
extern int nMax;//How far to calculate Bessels out to
extern int Z;//atomic number
extern int Energy;  //energy level to calculate (number of nodes for numerov) in eV
extern double Rmax; //Mesh radius limit. angstroms
extern double h;  //step size. angstroms
extern double k;  //angular wavenumber rad/m
extern const double m;  //mass. eV/c^2
extern const double PI;  //pi
extern const double h2m; //hbar^2/m eV-A^2
extern const double ee; //electron charge in eV-A
extern const double aBohr; //Bohr Radius in A
extern const double hbar; //hbar eV-s
extern const double ch; //electron charge C
extern const double SOL; //speed of light in m/s
extern std::vector<std::vector<double > > Potential;  //Potential Vector
