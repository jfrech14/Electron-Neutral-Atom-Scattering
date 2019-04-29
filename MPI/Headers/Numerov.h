#ifndef Numerov_H
#define Numerov_H

  double g(double V, double E, double r, int lam)
  {
      return((1.0/hbar2m)*( V + hbar2m*lam*(lam+1)/r/r - E));
  }

  std::vector<double> ForwardNumerov(double Ewant,std::vector<std::vector<double> > &CrossSec, std::vector<std::vector<double> >& Vpot, double dh, int lmax, double scale, double eEnergy, double k, int iter)
  {

      int i1, i2,check=0;
      double ddh12, norm, kappa, r1, r2, TotalCross ;
      int arrsize=Vpot.size()-1;

      /* allocate arrays */
      std::vector<long double> tandelta (lmax+1,0.0);
      std::vector<long double> delta (lmax+1,0.0);
      std::vector<long double> cross (lmax+1,0.0);
      std::vector<long double> u (arrsize+1,0.0);
      std::vector<long double> p (arrsize+1,0.0);
      std::vector<long double> f (arrsize+1,0.0);
      std::vector<double> DiffCross (181,0.0);
      std::vector<std::complex<double> > Amplitude(181,std::complex<double> (0.0,0.0));

      ddh12 = dh*dh/12.0;
      i1 = (int) ((scale-3)*((arrsize-1)/scale/2.0));
      r1 = Vpot[i1][0];
      i2 = (int) ((scale-1)*((arrsize-1)/scale/2.0));;
      r2 = Vpot[i2][0];

      for (int l = 0; l <= lmax; l++ )
      {
          cross[l] = 0.0 ;

          //  set up the f-function used by the Numerov algorithm
          for (int i = 0; i <= arrsize ; i++ )
          {
            f[i] = 1.0 - ddh12*g(Vpot[i][1],eEnergy,Vpot[i][0],l);
            u[i] = 0.0 ;
          }

          u[0] = 0;
          u[1] = dh/(1e50);
          u[2]=dh*dh*g(Vpot[1][1],eEnergy,Vpot[1][0],l)*u[1] +2.0*u[1];

          /* outward integration */
          for (int i = 2; i < arrsize; i++ )
          {
            u[i+1] = ((12.0-10.0*f[i])*u[i]-f[i-1]*u[i-1])/f[i+1];
            if (u[i+1]==0.0) return DiffCross;
          }

          //  normalization through summation integral
          norm = 0.0 ;
          for (int i = 0; i <= arrsize; i++ ) { norm += u[i]*u[i]*dh ;}
          for (int i = 0; i <= arrsize; i++ ) { u[i] = u[i]/sqrt(norm);}

          // Calculating phase shift
          kappa = r1*u[i2]/(r2*u[i1]);
          tandelta[l] = (kappa*Besselj(l,k*r1)-Besselj(l,k*r2)) / (kappa*Bessely(l,k*r1)-Bessely(l,k*r2));
          delta[l] = atan(tandelta[l]);
          if(std::isnan(delta[l])==1){delta[l]=0;}
          cross[l] = cross[l] + 4*PI/(k*k) * (double)(2*l+1)*sin(delta[l])*sin(delta[l]);

          // calculation of asymptotic wavefunction p
          for (int i = 0; i <= arrsize; i++ )
          {
            p[i]= sin (k*Vpot[i][0]- (double)l*PI/2.0 + delta[l]);
          }
          // normalize p so that it is equal to u for Vpot=r2
          for (int i = 0; i <= arrsize; i++ )
          {
            p[i] = p[i] / (p[i2]/u[i2]);
          }

          std::complex<double> cossin (cos(delta[l]),sin(delta[l]));
          //Calculate Differential Cross Section - first of two summations
          for (int theta=0;theta<=180;theta++)
          {
                std::complex<double> real((2*l+1)/k*sin(delta[l])*LegendreP(l,cos(theta*PI/180.0)),0.0);
                Amplitude[theta]+=real*cossin;
          }
      }

      //Calculate Differential Cross Section - last of two summations
      if (CrossSec[0][iter]>=Ewant-Ewant/20.0&&CrossSec[0][iter]<=Ewant+Ewant/20.0){check=1;}
      else{check=0;}
      for(int theta=0;theta<=180;theta++)
      {
          DiffCross[theta]=std::norm(Amplitude[theta]);
      }check=1;

      //"integrate" total cross section of all partial waves
      TotalCross =0.0;
      for (int a= 0; a <= lmax; a++ )
      {
        TotalCross += cross[a];
      }
      CrossSec[1][iter]=TotalCross;

      return(DiffCross);
  }

#endif
