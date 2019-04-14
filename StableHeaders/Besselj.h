#ifndef Besselj_H
#define Besselj_H 

    double Besseljup (int l, double x)
    {
        double j0=fabs(x)>1e-15 ? sin(x)/x : 1;
        if (l<=0) return j0;
        double j1=fabs(x)>1e-15 ? j0/x-cos(x)/x : x/3.;
        if (fabs(x)<1e-20&&l<=0) return 1;
        if (fabs(x)<1e-20&&l>0) return 0;
        double j2=j1;
        for (int i=2; i<=l; i++)
        {
            j2=j1*(2*i-1)/x - j0;
            j0=j1;
            j1=j2;
        }
        return j2;
    }

    double Besselj (int l, double x)
    {
        if (fabs(x)>l) return Besseljup(l,x);
        if (fabs(x)<1e-20&&l<=0) return 1;
        if (fabs(x)<1e-20&&l>0) return 0;
        int lstart = l + static_cast<int>(sqrt(40*l)/2.);
        double j2=0, j1=1;
        double j0,jl, x1=1/x;
        for (int i=lstart; i>=0; i--)
        {
            j0=(2*i+3.)*x1*j1-j2;
            if (i==l) jl=j0;
            {
                j2=j1;
                j1=j0;
            }

        }
        double true_j0=sin(x)/x;
        return jl*true_j0/j0;
    }

#endif