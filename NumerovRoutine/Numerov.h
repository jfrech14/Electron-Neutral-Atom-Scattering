#ifndef Numerov_H
#define Numerov_H

  double g(double V, double E, double r, int lam)
  {
      return((2*m/hbar/hbar/SOL/SOL)*( V + hbar*hbar*SOL*SOL/2/m*lam*(lam+1)/r/r - E));
  }

  std::vector<double> ForwardNumerov(std::vector<std::vector<double> >& Vpot,std::vector<std::vector<double> >& LegenP, double dh, int lmax, double scale, double eEnergy, double k)  
  {
    
    int i, l, i1, i2, n;
    double ddh12, norm, kappa, r1, r2, tandelta, delta, TotalCross ;
    
    int size=Vpot.size()-1;
    int psize=LegenP.size()-1;
   
    /* allocate arrays */

    std::vector<double> cross (lmax+1,0.0);
    std::vector<double> u (size,0.0);
    std::vector<double> p (size,0.0);
    std::vector<double> f (size,0.0);
    std::vector<double> DiffCross (181,0.0);
    std::vector<std::vector<double > > DiffCross1 (181,std::vector<double>(lmax+1,0));


      //omp_set_num_threads(threads);
      //#pragma omp parallel for shared(lmax, k,eEnergy,Vpot,LegenP,DiffCross1)
      for ( l = 0; l <= lmax; l++ ) 
      {
        cross[l] = 0.0 ;
        ddh12 = dh*dh/12.0;	
        i2 = size-1;
        r2 = Vpot[i2][0];
        i1 = int (0.75*i2);
        r1 = Vpot[i1][0];

        /*  set up the f-function used by the Numerov algorithm */

        for (int i = 0; i <= size ; i++ ) 
        {
          f[i] = 1.0 - ddh12*g(Vpot[i][1],eEnergy,Vpot[i][0],l);
          u[i] = 0.0 ;
        }
   
        u[0] = 0;
        u[1] = dh;
        u[2]=dh*dh*g(Vpot[2][1],eEnergy,Vpot[2][0],l)*u[1] +2*dh;
        /* outward integration */

        for (int i = 2; i < size; i++ ) 
        {
          u[i+1] = ((12.0-10.0*f[i])*u[i]-f[i-1]*u[i-1])/f[i+1];
        }

        /*  normalization through summation integral */

        norm = 0.0 ;
        for (int i = 0; i <= size; i++ ) { norm += u[i]*u[i]*dh ;}
        for (int i = 0; i <= size; i++ ) { u[i] = u[i]/sqrt(norm);}

        /* Calculating delta */

        kappa = r1*u[i2]/(r2*u[i1]) ;
        tandelta = (kappa*Bessely(l,k*r1)-Bessely(l,k*r2)) / (kappa*Bessely(l,k*r1)-Bessely(l,k*r2));
        delta = atan(tandelta);
        cross[l] = cross[l] + 4*PI/(k*k) * (double)(2*l+1)*sin(delta)*sin(delta);

        /*  calculation of asymptotic wavefunction p  */

        for (int i = 0; i <= size; i++ ) 
        {
          p[i]= sin (k*Vpot[i][0]- (double)l*PI/2.0 + delta); 
        }
        /* normalize p so that it is equal to u for Vpot=r2 */
        for (int i = 0; i <= size; i++ ) 
        {
          p[i] = p[i] / (p[i2]/u[i2] );
        }
        for (int theta=0;theta<=psize;theta++)
        {
          DiffCross1[theta][l]=1.0/k/k*(2*l+1)*(2*l+1)*LegenP[theta][l]*LegenP[theta][l]*sin(delta)*sin(delta);
          //std::cout<<LegenP[theta][l]<<std::endl;
        }
      }

      omp_set_num_threads(threads);
      #pragma omp parallel for shared(lmax, DiffCross1)
      for(int theta=0;theta<=180;theta++)
      {
          double sum=0;
          for (int n=0; n<=lmax; ++n)
          {
              sum+=DiffCross1[theta][n];
          }
          DiffCross[theta]=sum;
      }

      TotalCross =0.0;
      for (int a= 0; a <= lmax; a++ ) 
        { 
          TotalCross += cross[a];
        }
        std::cout<<TotalCross<<std::endl;

      return(DiffCross);
  }

#endif