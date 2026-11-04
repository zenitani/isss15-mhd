!-------------------------------------------------------------------------------
subroutine calc_flux_hlld(rol,vnl,vtl,vul,btl,bul,prl, &
                          ror,vnr,vtr,vur,btr,bur,prr, &
                          bnc,                         &
                          fro,fmn,fmt,fmu,fbt,fbu,fen)
!-------------------------------------------------------------------------------
  real(DP), intent(IN)  :: rol,vnl,vtl,vul,btl,bul,prl
  real(DP), intent(IN)  :: ror,vnr,vtr,vur,btr,bur,prr
  real(DP), intent(IN)  :: bnc
  real(DP), intent(OUT) :: fro,fmn,fmt,fmu,fbt,fbu,fen
  real(DP), parameter   :: eps=1.D-40
  real(DP)              :: bnc2,sgn
  real(DP)              :: roli,pml,ptl,enl,vbl
  real(DP)              :: rori,pmr,ptr,enr,vbr
  real(DP)              :: cl2,cal2,cbl2,cfl2,cfl
  real(DP)              :: cr2,car2,cbr2,cfr2,cfr
  real(DP               :: sr,sl
  real(DP)              :: slvl,srvr,rslvl,rsrvr,drsvi
  real(DP)              :: vnc,ptc
  real(DP)              :: slvc,rhdl,rhdli,rhnvl,rhnbl
  real(DP)              :: srvc,rhdr,rhdri,rhnvr,rhnbr
  real(DP)              :: ro2l,vt2l,vu2l,bt2l,bu2l,vb2l,en2l
  real(DP)              :: ro2r,vt2r,vu2r,bt2r,bu2r,vb2r,en2r
  real(DP)              :: rro2l,rro2r,rro2i
  real(DP)              :: vt3m,vu3m,bt3m,bu3m,vb3m,en3l,en3r
  real(DP)              :: rou,vtu,vuu,btu,buu,enu
!...............................................................................
! Bn at the interface
!...............................................................................
  bnc2 = bnc**2
  sgn  = sign(1.0_DP,bnc)
!...............................................................................
! variables at the left-face
!...............................................................................
  roli = 1.0_DP/rol
  pml  = 0.5_DP*(btl**2+bul**2)
  ptl  = prl+pml
  enl  = gammam1i*prl+pml+0.5_DP*rol*(vnl**2+vtl**2+vul**2)
  vbl  = vtl*btl+vul*bul
!...............................................................................
! variables at the right-face
!...............................................................................
  rori = 1.0_DP/ror
  pmr  = 0.5_DP*(btr**2+bur**2)
  ptr  = prr+pmr
  enr  = gammam1i*prr+pmr+0.5_DP*ror*(vnr**2+vtr**2+vur**2)
  vbr  = vtr*btr+vur*bur
!...............................................................................
! maximum / minimum wave speeds
!...............................................................................
  cl2  = gamma*prl*roli
  cr2  = gamma*prr*rori
  cal2 = bnc2*roli
  car2 = bnc2*rori
  cbl2 = cl2+cal2+2.0_DP*pml*roli
  cbr2 = cr2+car2+2.0_DP*pmr*rori
  cfl2 = 0.5_DP*(cbl2+sqrt(abs(cbl2**2-4.0_DP*cl2*cal2)))
  cfr2 = 0.5_DP*(cbr2+sqrt(abs(cbr2**2-4.0_DP*cr2*car2)))
  cfl  = sqrt(cfl2)
  cfr  = sqrt(cfr2)
  sl   = min(0.0_DP,min(vnl,vnr)-max(cfl,cfr))
  sr   = max(0.0_DP,max(vnl,vnr)+max(cfl,cfr))
!...............................................................................
! HLL average of the normal velocity and the total pressure
!...............................................................................
  slvl  = sl-vnl
  srvr  = sr-vnr
  rslvl = rol*slvl
  rsrvr = ror*srvr
  drsvi = 1.0_DP/(rsrvr-rslvl)
  vnc   = (rsrvr*vnr-rslvl*vnl-ptr+ptl)*drsvi
  ptc   = (rsrvr*ptl-rslvl*ptr+rsrvr*rslvl*(vnr-vnl))*drsvi
!...............................................................................
! variables of the outer sides in the Riemann fan
!...............................................................................
  slvc  = sl-vnc
  srvc  = sr-vnc
  ro2l  = rslvl/slvc
  ro2r  = rsrvr/srvc
  rhdl  = rslvl*slvc-bnc2
  rhdr  = rsrvr*srvc-bnc2
  if(abs(rhdl) > eps) then
    rhdli = 1.0_DP/rhdl
    rhnvl = (vnl-vnc)*bnc
    rhnbl = rslvl*slvl-bnc2
    vt2l = vtl+rhnvl*rhdli*btl
    vu2l = vul+rhnvl*rhdli*bul
    bt2l =     rhnbl*rhdli*btl
    bu2l =     rhnbl*rhdli*bul
  else
    vt2l = vtl
    vu2l = vul
    bt2l = btl
    bu2l = bul
  end if
  if(abs(rhdr) > eps) then
    rhdri = 1.0_DP/rhdr
    rhnvr = (vnr-vnc)*bnc
    rhnbr = rsrvr*srvr-bnc2
    vt2r = vtr+rhnvr*rhdri*btr
    vu2r = vur+rhnvr*rhdri*bur
    bt2r =     rhnbr*rhdri*btr
    bu2r =     rhnbr*rhdri*bur
  else
    vt2r = vtr
    vu2r = vur
    bt2r = btr
    bu2r = bur
  end if
  vb2l = vt2l*bt2l+vu2l*bu2l
  vb2r = vt2r*bt2r+vu2r*bu2r
  en2l = (slvl*enl-ptl*vnl+ptc*vnc+bnc*(vbl-vb2l))/slvc
  en2r = (srvr*enr-ptr*vnr+ptc*vnc+bnc*(vbr-vb2r))/srvc
!...............................................................................
! variables of the inner sides in the Riemann fan
!...............................................................................
  rro2l = sqrt(ro2l)
  rro2r = sqrt(ro2r)
  rro2i = 1.0_DP/(rro2r+rro2l)
  vt3m  = (rro2r*vt2r+rro2l*vt2l+            (bt2r-bt2l)*sgn)*rro2i
  vu3m  = (rro2r*vu2r+rro2l*vu2l+            (bu2r-bu2l)*sgn)*rro2i
  bt3m  = (rro2l*bt2r+rro2r*bt2l+rro2r*rro2l*(vt2r-vt2l)*sgn)*rro2i
  bu3m  = (rro2l*bu2r+rro2r*bu2l+rro2r*rro2l*(vu2r-vu2l)*sgn)*rro2i
  vb3m  = vt3m*bt3m+vu3m*bu3m
  en3l  = en2l-rro2l*(vb2l-vb3m)*sgn
  en3r  = en2r+rro2r*(vb2r-vb3m)*sgn
!...............................................................................
! variables at the interface
!...............................................................................
  if(vnc-abs(bnc)/rro2l > 0.0_DP) then
    rou = ro2l
    vtu = vt2l
    vuu = vu2l
    btu = bt2l
    buu = bu2l
    enu = en2l
  else if((vnc-abs(bnc)/rro2l <= 0.0_DP).and.(vnc >= 0.0_DP))then
    rou = ro2l
    vtu = vt3m
    vuu = vu3m
    btu = bt3m
    buu = bu3m
    enu = en3l
  else if((vnc < 0.0_DP).and.(vnc+abs(bnc)/rro2r >= 0.0_DP))then
    rou = ro2r
    vtu = vt3m
    vuu = vu3m
    btu = bt3m
    buu = bu3m
    enu = en3r
  else
    rou = ro2r
    vtu = vt2r
    vuu = vu2r
    btu = bt2r
    buu = bu2r
    enu = en2r
  end if
!...............................................................................
! HLLD fluxes
!...............................................................................
  fro = rou*vnc
  fmn = rou*vnc*vnc+ptc-bnc2*0.5_DP
  fmt = rou*vtu*vnc-bnc*btu
  fmu = rou*vuu*vnc-bnc*buu
  fbt = btu*vnc-bnc*vtu
  fbu = buu*vnc-bnc*vuu
  fen = (enu+ptc)*vnc-bnc*(vtu*btu+vuu*buu)
!-------------------------------------------------------------------------------
end subroutine calc_flux_hlld
!-------------------------------------------------------------------------------