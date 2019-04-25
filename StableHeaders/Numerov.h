#ifndef Numerov_H
#define Numerov_H

  double g(double V, double E, double r, int lam)
  {
      return((1.0/hbar2m)*( V + hbar2m*lam*(lam+1)/r/r - E));
  }

  std::vector<double> ForwardNumerov(double Ewant, std::vector<std::vector<double> > &PartialWaves,std::vector<std::vector<double> > &CrossSec, std::vector<std::vector<double> >& Vpot, double dh, int lmax, double scale, double eEnergy, double k, int iter)
  {

      int i1, i2,check=0;
      double ddh12, norm, kappa, r1, r2, TotalCross ;
      int size=Vpot.size()-1;

      /* allocate arrays */
      std::vector<long double> tandelta (lmax+1,0.0);
      std::vector<long double> delta (lmax+1,0.0);
      std::vector<double> cross (lmax+1,0.0);
      std::vector<long double> u (size+1,0.0);
      std::vector<long double> p (size+1,0.0);
      std::vector<long double> f (size+1,0.0);
      std::vector<double> DiffCross (181,0.0);
      std::vector<std::vector<double > > DiffCross1 (181,std::vector<double>(lmax+1,0));

      ddh12 = dh*dh/12.0;
      i1 = (int) ((scale-3)*((size-1)/scale/2.0));
      r1 = Vpot[i1][0];
      i2 = (int) ((scale-1)*((size-1)/scale/2.0));;
      r2 = Vpot[i2][0];

      for (int l = 0; l <= lmax; l++ )
      {
          cross[l] = 0.0 ;

          //  set up the f-function used by the Numerov algorithm
          for (int i = 0; i <= size ; i++ )
          {
            f[i] = 1.0 - ddh12*g(Vpot[i][1],eEnergy,Vpot[i][0],l);
          }

          u[0] = 0;
          u[1] = dh/(1e50);
          u[2]=dh*dh*g(Vpot[1][1],eEnergy,Vpot[1][0],l)*u[1] +2.0*u[1];

          /* outward integration */
          for (int i = 2; i < size; i++ )
          {
            u[i+1] = ((12.0-10.0*f[i])*u[i]-f[i-1]*u[i-1])/f[i+1];
            if (u[i+1]==0.0) return DiffCross;
          }

          //  normalization through summation integral
          norm = 0.0 ;
          for (int i = 0; i <= size; i++ ) { norm += u[i]*u[i]*dh ;}
          for (int i = 0; i <= size; i++ ) { u[i] = u[i]/sqrt(norm);}

          // Calculating phase shift
          kappa = r1*u[i2]/(r2*u[i1]);
          tandelta[l] = (kappa*Besselj(l,k*r1)-Besselj(l,k*r2)) / (kappa*Bessely(l,k*r1)-Bessely(l,k*r2));
          delta[l] = atan(tandelta[l]);
          if(std::isnan(delta[l])==1){delta[l]=0;}
          cross[l] = cross[l] + 4*PI/(k*k) * (double)(2*l+1)*sin(delta[l])*sin(delta[l]);

          // calculation of asymptotic wavefunction p
          for (int i = 0; i <= size; i++ )
          {
            p[i]= sin (k*Vpot[i][0]- (double)l*PI/2.0 + delta[l]);
          }
          // normalize p so that it is equal to u for Vpot=r2
          for (int i = 0; i <= size; i++ )
          {
            p[i] = p[i] / (p[i2]/u[i2]);
          }

          //Calculate Differential Cross Section - first of two summations
          for (int theta=0;theta<=180;theta++)
          {
            for (int ll=0;ll<=lmax;ll++)
              {
                DiffCross1[theta][l]+=1.0/k/k*(2.0*l+1)*(2.0*ll+1)*LegendreP(ll,cos(theta*PI/180.0))*LegendreP(l,cos(theta*PI/180.0))*sin(delta[l])*sin(delta[ll]);
              }
          }
      }

      //Calculate Differential Cross Section - last of two summations
      for(int theta=0;theta<=180;theta++)
      {
          for (int n=0; n<=lmax; ++n)
          {
              DiffCross[theta]+=DiffCross1[theta][n];
              if ((int)CrossSec[iter][0]==(int)Ewant&&check!=1){PartialWaves[theta][n]=DiffCross1[theta][n];}
          }
      }check=1;

      //"integrate" total cross section of all partial waves
      TotalCross =0.0;
      for (int a= 0; a <= lmax; a++ )
      {
        TotalCross += cross[a];
      }
      CrossSec[iter][1]=TotalCross;

      return(DiffCross);
  }

#endif
