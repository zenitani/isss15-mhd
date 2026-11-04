#include <math.h>
inline double max(double a, double b)
{
  return( (a > b)?a:b );
}
inline double min(double a, double b)
{
  return( (a > b)?b:a );
}

#define EPS (1e-8)

void calc_flux_hlld(double rol, double vnl, double vtl, double vul, double btl,
                    double bul, double prl,
                    double ror, double vnr, double vtr, double vur, double btr,
                    double bur, double prr,
                    double bnc, double gamma,
                    double *fro, double *fmn, double *fmt, double *fmu, double *fbt,
                    double *fbu, double *fen)
/* Calculate HLLD fluxes */
/* rol,vnl,vtl,vul,btl,bul,prl: input primitive variables at the left side */
/* ror,vnr,vtr,vur,btr,bur,prr: input primitive variables at the right side */
/* bnc: input normal magnetic field at the interface */
/* gamma: specific heat ratio */
/* fro,fmn,fmt,fmu,fbt,fbu,fen: output HLLD fluxes at the interface*/
{
  double gammam1i=1.0/(gamma-1.0);
  /* Bn at the interface */
  double bnc2=bnc*bnc;
  int sgn=(bnc > 0)?(1):(-1);
  /* Variables at the left-face */
  double roli=1.0/rol;
  double pml=0.5*(btl*btl+bul*bul);
  double ptl=prl+pml;
  double enl=gammam1i*prl+pml+0.5*rol*(vnl*vnl+vtl*vtl+vul*vul);
  double vbl=vtl*btl+vul*bul;
  /* Variables at the right-face */
  double rori=1.0/ror;
  double pmr=0.5*(btr*btr+bur*bur);
  double ptr=prr+pmr;
  double enr=gammam1i*prr+pmr+0.5*ror*(vnr*vnr+vtr*vtr+vur*vur);
  double vbr=vtr*btr+vur*bur;
  /* Maximum/minimum wave speeds */
  double cl2=gamma*prl*roli;
  double cr2=gamma*prr*rori;
  double cal2=bnc2*roli;
  double car2=bnc2*rori;
  double cbl2=cl2+cal2+2.0*pml*roli;
  double cbr2=cr2+car2+2.0*pmr*rori;
  double cfl2=0.5*(cbl2+sqrt(fabs(cbl2*cbl2-4.0*cl2*cal2)));
  double cfr2=0.5*(cbr2+sqrt(fabs(cbr2*cbr2-4.0*cr2*car2)));
  double cfl=sqrt(cfl2);
  double cfr=sqrt(cfr2);
  double sl=min(0.0,min(vnl,vnr)-max(cfl,cfr));
  double sr=max(0.0,max(vnl,vnr)+max(cfl,cfr));
  /* HLL average of the normal velocity and the total pressure */
  double slvl=sl-vnl;
  double srvr=sr-vnr;
  double rslvl=rol*slvl;
  double rsrvr=ror*srvr;
  double drsvi=1.0/(rsrvr-rslvl);
  double vnc=(rsrvr*vnr-rslvl*vnl-ptr+ptl)*drsvi;
  double ptc=(rsrvr*ptl-rslvl*ptr+rsrvr*rslvl*(vnr-vnl))*drsvi;
  /* Variables of the outer sides in the Riemann fan */
  double slvc=sl-vnc;
  double srvc=sr-vnc;
  double ro2l=rslvl/slvc;
  double ro2r=rsrvr/srvc;
  double rhdl=rslvl*slvc-bnc2;
  double rhdr=rsrvr*srvc-bnc2;
  double vt2l,vu2l,bt2l,bu2l;
  double vt2r,vu2r,bt2r,bu2r;
  if (fabs(rhdl) > EPS){
    double rhdli=1.0/rhdl;
    double rhnvl=(vnl-vnc)*bnc;
    double rhnbl=rslvl*slvl-bnc2;
    vt2l=vtl+rhnvl*rhdli*btl;
    vu2l=vul+rhnvl*rhdli*bul;
    bt2l=rhnbl*rhdli*btl;
    bu2l=rhnbl*rhdli*bul;
  } else{
    vt2l=vtl;
    vu2l=vul;
    bt2l=btl;
    bu2l=bul;
  }
  if (fabs(rhdr) > EPS){
    double rhdri=1.0/rhdr;
    double rhnvr=(vnr-vnc)*bnc;
    double rhnbr=rsrvr*srvr-bnc2;
    vt2r=vtr+rhnvr*rhdri*btr;
    vu2r=vur+rhnvr*rhdri*bur;
    bt2r=rhnbr*rhdri*btr;
    bu2r=rhnbr*rhdri*bur;
  } else{
    vt2r=vtr;
    vu2r=vur;
    bt2r=btr;
    bu2r=bur;
  }
  double vb2l=vt2l*bt2l+vu2l*bu2l;
  double vb2r=vt2r*bt2r+vu2r*bu2r;
  double en2l=(slvl*enl-ptl*vnl+ptc*vnc+bnc*(vbl-vb2l))/slvc;
  double en2r=(srvr*enr-ptr*vnr+ptc*vnc+bnc*(vbr-vb2r))/srvc;
  /* Variables of the inner sides in the Riemann fan */
  double rro2l=sqrt(ro2l);
  double rro2r=sqrt(ro2r);
  double rro2i=1.0/(rro2r+rro2l);
  double vt3m=(rro2r*vt2r+rro2l*vt2l+(bt2r-bt2l)*sgn)*rro2i;
  double vu3m=(rro2r*vu2r+rro2l*vu2l+(bu2r-bu2l)*sgn)*rro2i;
  double bt3m=(rro2l*bt2r+rro2r*bt2l+rro2r*rro2l*(vt2r-vt2l)*sgn)*rro2i;
  double bu3m=(rro2l*bu2r+rro2r*bu2l+rro2r*rro2l*(vu2r-vu2l)*sgn)*rro2i;
  double vb3m=vt3m*bt3m+vu3m*bu3m;
  double en3l=en2l-rro2l*(vb2l-vb3m)*sgn;
  double en3r=en2r+rro2r*(vb2r-vb3m)*sgn;
 /* Variables at the interface */
  double rou,vtu,vuu,btu,buu,enu;
  if (vnc-fabs(bnc)/rro2l > 0){
    rou=ro2l;
    vtu=vt2l;
    vuu=vu2l;
    btu=bt2l;
    buu=bu2l;
    enu=en2l;
  } else{
    if (vnc >= 0){
      rou=ro2l;
      vtu=vt3m;
      vuu=vu3m;
      btu=bt3m;
      buu=bu3m;
      enu=en3l;
    } else{
      if (vnc+fabs(bnc)/rro2r >= 0){
       rou=ro2r;
       vtu=vt3m;
       vuu=vu3m;
       btu=bt3m;
       buu=bu3m;
       enu=en3r;
      } else{
       rou=ro2r;
       vtu=vt2r;
       vuu=vu2r;
       btu=bt2r;
       buu=bu2r;
       enu=en2r;
      }
    }
  }
  /* HLLD fluxes */
  *fro=rou*vnc;
  *fmn=rou*vnc*vnc+ptc-bnc2*0.5;
  *fmt=rou*vtu*vnc-bnc*btu;
  *fmu=rou*vuu*vnc-bnc*buu;
  *fbt=btu*vnc-bnc*vtu;
  *fbu=buu*vnc-bnc*vuu;
  *fen=(enu+ptc)*vnc-bnc*(vtu*btu+vuu*buu);
}
