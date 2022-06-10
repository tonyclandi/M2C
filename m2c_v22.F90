!
! Version 22.2
!
! Input data: moloch mhf data file, file with geografical coordinates (lat and
! lon) in 2 ascii format file (chimere_lat.dat and chimere_lon.dat),
! and auxilary input namelist file param_geo.inp with point number of
! geografical coordinates
!
! In output: soilt, that are soil temperature at the levels
! 3,5, 17,5, 64,0, 194,5 cm under surface;
! and winv, that is vertical velocity at integer moloch levels.
!
! ### Sezione dei moduli ###
 module mod_moloch

 integer :: nfdr(50)
 real    :: pdr(100)
 integer, parameter :: nlev_snow=11

 integer, dimension(5,100) :: idatec
 integer :: idate0(5), iperiod(3)
 integer :: nlon, nlat, nlev, nlevp1, nlevg, njump=1, ierr_read1=0, ierr_read2=0
 real, dimension(:), allocatable :: hxv, hxt
 real, dimension(:,:), allocatable :: alatt, alont, alatu, alonu, alatv, alonv, work
 real, dimension(3) :: lat_min, lat_max, lon_min, lon_max
 real :: dtstep

 !physical constants

 real :: rd, ra, cp, cpv, rcp, omd, eps, ep, g, we, rv, alp0, pi, ttr, ak, h, b0

 real, dimension(:,:,:), allocatable :: u, v, wh, t, tvirt, p, q, qcw, qci, tg, qg, qgmax, qgmin, fice_soil, &
 snow_lev, snow_t, snow_fice, snow_age, snow_melt_age, snow_dens, &
 zeta, pre
 real, dimension(:,:), allocatable :: prectot, snow, precconv, albedo, snow_albedo, emis1, emis2, &
 rgm, rgq, fmask, precsolid, runoff, runoff_tot, tskin, tgsurf, qgsurf, fice_soil_surf, &
 cswfl, clwfl, chflux, cqflux, cwvflux, csoilh_bott_flux, csoilw_bott_flux, &
 t2min, t2max, ws10max, fice, iceth, phig, ps, &
 hflux, qflux, qskin, cloudt, t2, q2, q2rel, u10, v10, ustar, pbl, richs, lspre
 real, dimension(:), allocatable :: soillev

 contains

    function gzita (zita)    ! Decay function
    gzita = 1. -3.*(zita/h)**2 +2.*(zita/h)**3
    end function gzita

    function bzita (zita)    ! Stretching function
    bzita = b0 + (1.-b0)*(zita/h)
    end function bzita

 end module mod_moloch

 module mod_chimere

 integer,parameter   :: ndstr = 19, nsoillev_out=4
 real, dimension(nsoillev_out) :: soillev_out=(/0.035, 0.175, 0.64, 1.945/)
 character(len=19),dimension(100)   :: Times ! la var tempo in formato "CHIMERE"
 character(len=4)    :: varnames

 integer :: nx, ny

 ! id files
 integer :: ncout, ncstat
 ! dimensioni
 integer :: nlon_did, nlat_did, lev_did, time_did, dstr_did
 ! variabili 2D
 integer :: time_vid, psfc_vid, swrd_vid, lwrd_vid, sshf_vid, slhf_vid, usta_vid, tem2_vid, soim_vid, hght_vid, rh2m_vid, &
                                          lspc_vid, copc_vid, u10m_vid, v10m_vid, weas_vid, stl1_vid 
 ! variabili 3D
 integer :: temp_vid, cliq_vid, rain_vid, sphu_vid, pres_vid, alti_vid, winz_vid, winm_vid, winw_vid

 !character(len=*), parameter :: meteo_out = "/data/landi/exdomout_moloch.nc"
 character(len=*), parameter :: meteo_out = "/opt/M2C_2022/exdomout_moloch_22.nc"
 real, dimension(:,:),     allocatable :: lon, lat, lonr, latr, uf, vf
 real, dimension(:,:,:),   allocatable :: psfc, swrd, lwrd, sshf, slhf, usta, tem2, soim, hght, rh2m, &
                                          lspc, copc, u10m, v10m, weas, stl1
 real, dimension(:,:,:,:), allocatable :: temp, cliq, rain, sphu, pres, alti, winz, winv, winm, soilt
 !real :: dx, dy, xmin, ymin

 namelist/outgrid/ nx, ny

 end module mod_chimere

! ### Fine sezione dei moduli ###

! ### Inizio programma principale ###
program m2c

 use mod_moloch
 use mod_chimere
 use netcdf

 integer :: nunic1=21, nunic2=22, ierr_open1=0, ierr_open2=0, ivalt=0
 integer :: nunic_grid=11, ierr_open_grid=0

  real, dimension(:,:), allocatable :: prectotr, precconvr, precsolidr, runoffr, runoff_tot_r, &
 cswflr, clwflr, chfluxr, cqfluxr, t2minr, t2maxr, ws10maxr, &
 cwvfluxr, csoilh_bott_fluxr, csoilw_bott_fluxr
  real, dimension(:,:,:), allocatable :: w, tg_levout
  real, dimension(200) :: zaxis_inp, zvar_inp, zaxis_out, zvar_out

 ! tcl, 28/06/18: creato dallo script bash prepara_chimere_var.sh
 ! per adesso ci sono solo le variabili relative alla griglia di
 ! chimere (dato il nome del dominio di simulazione presente nel
 ! file domaninlist.nml di chimere, eg EUR20)
 ! include "chimere_vars.h"

! Define physical parameters

 rd     = 287.05
 rv     = 461.51
 ra     = 6371.e+3
 cp     = 1004.6
 cpv    = 1869.46
 zdelta = cpv/cp
 omd    = zdelta-1.
 eps    = rd/rv
 ep     = rv/rd -1.
 g      = 9.807
 pi     = abs(acos(-1.))
 alp0   = log(1.e5)
 zstef  = 5.67e-8
 ttr    = 273.15
 ak     = 0.4
 rcp    = rd/cp

 open (nunic1,file='moloch_atm.mhf',status='old',form='unformatted',iostat=ierr_open1)
 print *
 if (ierr_open1 /= 0) then
   print *,'Not found input file moloch_atm.mhf'
   stop
 else
   print *,'Input file moloch_atm.mhf opened on the unit ',nunic1
 endif
 open (nunic2,file='moloch_soil.mhf',status='old',form='unformatted',iostat=ierr_open2)
 print *
 if (ierr_open2 /= 0) then
   print *,'Not found input file moloch_soil.mhf'
   stop
 else
   print *,'Input file moloch_soil.mhf opened on the unit ',nunic2
 endif
 print *
 read (nunic1) nfdr
 read (nunic1) pdr

 rewind nunic1

 nlon =nfdr(2)
 nlat =nfdr(3)
 nlev =nfdr(4)
 nlevg=nfdr(15)
 nhist=nfdr(17)

 zdlat  = pdr(1)
 zdlon  = pdr(2)
 dtstep = pdr(3)
 zalat0 = pdr(4) ! zalat0 is the southernmost latitude (v-point)
 zalon0 = pdr(5) ! zalon0 is the westernmost longitude (v-point)
 if (zalon0 > 180.) zalon0=zalon0-360.
 x0 = pdr(39)    ! centro griglia ruotata in gradi
 y0 = pdr(38)
 h  = pdr(40)
 b0 = pdr(42)

 allocate(soillev(nlevg+1), stat=ierr)
 soillev(1)         = 0.
 soillev(2:nlevg+1) = pdr(6:5+nlevg)

 nlevp1 = nlev+1
 ztcu = dtstep*nhist*njump  ! tempo di accumulo
   print*, "ztcu = ", ztcu

 ! stampo parametri di griglia moloch
 print*, " >>> Parametri griglia MOLOCH"
 !write(*, '(a,3i5)'), " nlon, nlat, nlev:", nlon, nlat, nlev
 print*, " nlon, nlat, nlev:", nlon, nlat, nlev

! Dinamic array allocation:

 allocate(hxv(nlat), stat=ierr)
 allocate(hxt(nlat), stat=ierr)
 allocate(alatt(nlon,nlat), stat=ierr)
 allocate(alont(nlon,nlat), stat=ierr)
 allocate(alatu(nlon,nlat), stat=ierr)
 allocate(alonu(nlon,nlat), stat=ierr)
 allocate(alatv(nlon,nlat), stat=ierr)
 allocate(alonv(nlon,nlat), stat=ierr)
 allocate(work(nlon,nlat), stat=ierr)
 allocate(phig(nlon,nlat), stat=ierr)
 allocate(zeta(nlon,nlat,nlev), stat=ierr)
 allocate(t(nlon,nlat,nlev), stat=ierr)
 allocate(p(nlon,nlat,nlev), stat=ierr)
 allocate(tvirt(nlon,nlat,nlev), stat=ierr)
 allocate(pre(nlon,nlat,nlev), stat=ierr)
 allocate(u(nlon,nlat,nlev), stat=ierr)
 allocate(v(nlon,nlat,nlev), stat=ierr)
 allocate(wh(nlon,nlat,nlev+1), stat=ierr)
 allocate(q(nlon,nlat,nlev), stat=ierr)
 allocate(qcw(nlon,nlat,nlev), stat=ierr)
 allocate(qci(nlon,nlat,nlev), stat=ierr)
 allocate(ps(nlon,nlat), stat=ierr)
 allocate(t2(nlon,nlat), stat=ierr)
 allocate(q2(nlon,nlat), stat=ierr)
 allocate(q2rel(nlon,nlat), stat=ierr)
 allocate(u10(nlon,nlat), stat=ierr)
 allocate(v10(nlon,nlat), stat=ierr)
 allocate(ustar(nlon,nlat), stat=ierr)

 allocate(tg(nlon,nlat,nlevg), stat=ierr)
 allocate(qg(nlon,nlat,nlevg), stat=ierr)
 allocate(qgmax(nlon,nlat,nlevg), stat=ierr)
 allocate(qgmin(nlon,nlat,nlevg), stat=ierr)
 allocate(fice_soil(nlon,nlat,nlevg), stat=ierr)

 allocate(     snow_lev (nlon, nlat, nlev_snow), stat=ierr)
 allocate(       snow_t (nlon, nlat, nlev_snow), stat=ierr)
 allocate(    snow_fice (nlon, nlat, nlev_snow), stat=ierr)
 allocate(     snow_age (nlon, nlat, nlev_snow), stat=ierr)
 allocate(snow_melt_age (nlon, nlat, nlev_snow), stat=ierr)
 allocate(    snow_dens (nlon, nlat, nlev_snow), stat=ierr)

 allocate(prectot(nlon,nlat), stat=ierr)
 allocate(lspre(nlon,nlat), stat=ierr)
 allocate(snow(nlon,nlat), stat=ierr)
 allocate(runoff(nlon,nlat), stat=ierr)
 allocate(runoff_tot(nlon,nlat), stat=ierr)
 allocate(precconv(nlon,nlat), stat=ierr)
 allocate(albedo(nlon,nlat), stat=ierr)
 allocate(snow_albedo(nlon,nlat), stat=ierr)
 allocate(emis1(nlon,nlat), stat=ierr)
 allocate(emis2(nlon,nlat), stat=ierr)
 allocate(rgm(nlon,nlat), stat=ierr)
 allocate(rgq(nlon,nlat), stat=ierr)
 allocate(fmask(nlon,nlat), stat=ierr)
 allocate(precsolid(nlon,nlat), stat=ierr)
 allocate(tskin(nlon,nlat), stat=ierr)
 allocate(tgsurf(nlon,nlat), stat=ierr)
 allocate(qgsurf(nlon,nlat), stat=ierr)
 allocate(fice_soil_surf(nlon,nlat), stat=ierr)
 allocate(cswfl(nlon,nlat), stat=ierr)
 allocate(clwfl(nlon,nlat), stat=ierr)
 allocate(chflux(nlon,nlat), stat=ierr)
 allocate(cqflux(nlon,nlat), stat=ierr)
 allocate(cwvflux(nlon, nlat), stat=ierr)
 allocate(csoilh_bott_flux(nlon, nlat), stat=ierr)
 allocate(csoilw_bott_flux(nlon, nlat), stat=ierr)
 allocate(t2min(nlon,nlat), stat=ierr)
 allocate(t2max(nlon,nlat), stat=ierr)
 allocate(ws10max(nlon,nlat), stat=ierr)
 allocate(fice(nlon,nlat), stat=ierr)
 allocate(iceth(nlon,nlat), stat=ierr)
 allocate(hflux (nlon,nlat), stat=ierr)
 allocate(qflux (nlon,nlat), stat=ierr)
 allocate(qskin(nlon,nlat), stat=ierr)
 allocate(cloudt(nlon,nlat), stat=ierr)
 allocate(pbl(nlon,nlat), stat=ierr)
 allocate(richs(nlon,nlat), stat=ierr)

 allocate(prectotr   (nlon, nlat), stat=ierr)
 allocate(precconvr  (nlon, nlat), stat=ierr)
 allocate(precsolidr (nlon, nlat), stat=ierr)
 allocate(runoffr (nlon, nlat), stat=ierr)
 allocate(runoff_tot_r (nlon, nlat), stat=ierr)
 allocate( cswflr (nlon, nlat), stat=ierr)
 allocate( clwflr (nlon, nlat), stat=ierr)
 allocate(chfluxr (nlon, nlat), stat=ierr)
 allocate(cqfluxr (nlon, nlat), stat=ierr)
 allocate( t2minr (nlon, nlat), stat=ierr)
 allocate( t2maxr (nlon, nlat), stat=ierr)
 allocate(ws10maxr(nlon, nlat), stat=ierr)
 allocate(cwvfluxr(nlon, nlat), stat=ierr)
 allocate(csoilh_bott_fluxr(nlon, nlat), stat=ierr)
 allocate(csoilw_bott_fluxr(nlon, nlat), stat=ierr)

 allocate(w(nlon,nlat,nlev), stat=ierr)
 allocate(tg_levout(nlon,nlat,nsoillev_out), stat=ierr)

! Definition of real grid point latitudes and longitudes (in degrees) using double precision
! Definition of geographical factors hxv, hxt

 call rot_grid (pdr(39),pdr(38),zalon0,          zalat0+zdlat*.5, zdlon,zdlat,alont,alatt,nlon,nlat)
 call rot_grid (pdr(39),pdr(38),zalon0+zdlon*.5, zalat0+zdlat*.5, zdlon,zdlat,alonu,alatu,nlon,nlat)
 call rot_grid (pdr(39),pdr(38),zalon0,          zalat0          ,zdlon,zdlat,alonv,alatv,nlon,nlat)

!  Definition of point number of output grid

 open (nunic_grid, file='param_geo.inp', status='old', iostat=ierr_open_grid)
 if (ierr_open_grid /= 0) then
   print *,"Not found input file param_geo.inp with point number of output grid, stop"
   stop
 endif
 read (nunic_grid, outgrid)
 close (nunic_grid)

 print *
 Print *," >>> Point number of CHIMERE grid (nx,ny): ",nx, ny

! Geographical coordinates of CHIMERE grid (regolar lat and lon)

 allocate(lon(nx,ny), stat=ierr)
 allocate(lat(nx,ny), stat=ierr)

! Reading of input data about CHIMERE grib lon and lat in regolar non rotated coordinate system

 ierr_open_grid = 0
 ierr_read_grid = 0

 open (nunic_grid, file='chimere_lat.dat', status='old', iostat=ierr_open_grid)
 if (ierr_open_grid /= 0) then
   print *,"Not found input file chimere_lat.dat with latitude of output grid, stop"
   stop
 endif

 do j = 1,ny
 do i = 1,nx
   read (nunic_grid, *, iostat=ierr_read_grid) lat(i,j)
   if (ierr_read_grid /= 0) then
     print *,"End of file chimere_lat.dat during reading, stop"
     stop
   endif
 enddo
 enddo

 close (nunic_grid)

 ierr_open_grid = 0
 ierr_read_grid = 0

 open (nunic_grid, file='chimere_lon.dat', status='old', iostat=ierr_open_grid)
 if (ierr_open_grid /= 0) then
   print *,"Not found input file chimere_lon.dat with latitude of output grid, stop"
   stop
 endif

 do j = 1,ny
 do i = 1,nx
   read (nunic_grid, *, iostat=ierr_read_grid) lon(i,j)
   if (ierr_read_grid /= 0) then
     print *,"End of file chimere_lon.dat during reading, stop"
     stop
   endif
 enddo
 enddo

 close (nunic_grid)

! Geographical coordinates of CHIMERE grid (lat and lon rotated in Moloch coord. system)

 allocate(lonr(nx,ny), stat=ierr)
 allocate(latr(nx,ny), stat=ierr)

! Definizione di lon e lat di chimere in coordinate ruotate (in lonr, latr)

 call anti_rot_grid (x0, y0, lon, lat, lonr, latr, nx, ny)

 ! 2D
 allocate(uf(nx,ny), stat=ierr)
 allocate(vf(nx,ny), stat=ierr)
 allocate(psfc(nx,ny,100), stat=ierr)
 allocate(swrd(nx,ny,100), stat=ierr)
 allocate(lwrd(nx,ny,100), stat=ierr)
 allocate(sshf(nx,ny,100), stat=ierr)
 allocate(slhf(nx,ny,100), stat=ierr)
 allocate(usta(nx,ny,100), stat=ierr)
 allocate(tem2(nx,ny,100), stat=ierr)
 allocate(soim(nx,ny,100), stat=ierr)

 allocate(hght(nx,ny,100), stat=ierr)
 allocate(rh2m(nx,ny,100), stat=ierr)
 allocate(lspc(nx,ny,100), stat=ierr)
 allocate(copc(nx,ny,100), stat=ierr)
 allocate(u10m(nx,ny,100), stat=ierr)
 allocate(v10m(nx,ny,100), stat=ierr)
 allocate(weas(nx,ny,100), stat=ierr)
 allocate(stl1(nx,ny,100), stat=ierr)

 ! 3D
 allocate(temp(nx,ny,nlev,100), stat=ierr)
 allocate(cliq(nx,ny,nlev,100), stat=ierr)
 allocate(rain(nx,ny,nlev,100), stat=ierr)
 allocate(sphu(nx,ny,nlev,100), stat=ierr)
 allocate(pres(nx,ny,nlev,100), stat=ierr)
 allocate(alti(nx,ny,nlev,100), stat=ierr)
 allocate(winz(nx,ny,nlev,100), stat=ierr)
 allocate(winv(nx,ny,nlev,100), stat=ierr)
 allocate(winm(nx,ny,nlev,100), stat=ierr)

 allocate(soilt(nx,ny,nsoillev_out,100), stat=ierr)

call rd_param_const(nlon, nlat, nlevg, x0, y0, zalon0, zalat0+zdlat*0.5, zdlon, zdlat, phig, fmask, qgmax, qgmin)

! inizio loop temporale di lettura mhf di moloch

iist  = 0
iist2 =0

prectot(:,:) = 0.
precconv(:,:) = 0.
precsolid(:,:) = 0.
runoff(:,:) = 0.
runoff_tot(:,:) = 0.
cswfl(:,:)  = 0.
clwfl(:,:)  = 0.
chflux(:,:) = 0.
cqflux(:,:) = 0.
t2min(:,:)  = 999.
t2max(:,:)  = 0.
ws10max(:,:)= 0.
cwvflux(:,:)= 0.
csoilh_bott_flux(:,:)= 0.
csoilw_bott_flux(:,:)= 0.

do while (.true.)

! Read mhf file

  print *
  print *, 'Read ',iist+1,' instant from moloch_atm.mhf and moloch_soil.mhf'

  call rdmhf_atm(nunic1, nlon, nlat, nlev, nfdr, pdr, p, u, v, wh, t, q, qcw, qci)
  if (ierr_read1 /= 0) exit

  call rdmhf_soil(nunic2, nlon, nlat, nlevg, nlev_snow, nfdr, pdr, &
                  rgm, rgq, fice, iceth, albedo, emis1, emis2, &
                  cloudt, prectotr, precconvr, precsolidr, &
                  tskin, tgsurf, tg, qskin, qgsurf, qg, fice_soil_surf, fice_soil, &
                  snow, snow_lev, snow_t, snow_fice, snow_age, snow_melt_age, snow_dens, snow_albedo, &
                  cswflr, clwflr, chfluxr, cqfluxr, t2minr, t2maxr, ws10maxr, runoffr, runoff_tot_r, &
                  cwvfluxr, csoilh_bott_fluxr, csoilw_bott_fluxr)
  if (ierr_read2 /= 0) exit

  idate0(1:5)=nfdr(5:9)
  iperiod(1:3)=nfdr(10:12)
  ivalt = nfdr(10)*10000+nfdr(11)*100+nfdr(12)

  if (iist == 0) then
    if (ivalt == 0) print*, "The 1st processed instant is an initial condition"
    if (njump /= 1) print '(A8,I2,A53)'," njump =",njump,": note that the 1st forecast instant is not processed"
  endif

  iist = iist + 1

  do jlat=1,nlat
  do jlon=1,nlon
    prectot(jlon,jlat)=prectot(jlon,jlat)+prectotr(jlon,jlat)
    precconv(jlon,jlat)=precconv(jlon,jlat)+precconvr(jlon,jlat)
    precsolid(jlon,jlat)=precsolid(jlon,jlat)+precsolidr(jlon,jlat)
    runoff(jlon,jlat)=runoff(jlon,jlat)+runoffr(jlon,jlat)
    runoff_tot(jlon,jlat)=runoff_tot(jlon,jlat)+runoff_tot_r(jlon,jlat)
    cswfl(jlon,jlat)=cswfl(jlon,jlat)+cswflr(jlon,jlat)
    clwfl(jlon,jlat)=clwfl(jlon,jlat)+clwflr(jlon,jlat)
    chflux(jlon,jlat)=chflux(jlon,jlat)+chfluxr(jlon,jlat)
    cqflux(jlon,jlat)=cqflux(jlon,jlat)+cqfluxr(jlon,jlat)
    if (t2minr(jlon,jlat)<t2min(jlon,jlat)) t2min(jlon,jlat)=t2minr(jlon,jlat)
    if (t2maxr(jlon,jlat)>t2max(jlon,jlat)) t2max(jlon,jlat)=t2maxr(jlon,jlat)
    if (ws10maxr(jlon,jlat)>ws10max(jlon,jlat)) ws10max(jlon,jlat)=ws10maxr(jlon,jlat)
    cwvflux(jlon,jlat)=cwvflux(jlon,jlat)+cwvfluxr(jlon,jlat)
    csoilh_bott_flux(jlon,jlat)=csoilh_bott_flux(jlon,jlat)+csoilh_bott_fluxr(jlon,jlat)
    csoilw_bott_flux(jlon,jlat)=csoilw_bott_flux(jlon,jlat)+csoilw_bott_fluxr(jlon,jlat)
  enddo
  enddo

 if (iist==1) then
   if (ivalt==0) then
     iist0=-1 ! validity time of first field=0: analysis
   else
     iist0=0
   endif
 endif

 if (mod(iist+iist0,njump) == 0) then ! ---> (very long if)
   iist2=iist2+1 ! calculator of output instants

! Calculation and writing of forecast data

 call calendar (nfdr(5), nfdr(6), nfdr(7), nfdr(8), nfdr(9), nfdr(10), nfdr(11), nfdr(12), idatec(1,iist2), &
                idatec(2,iist2), idatec(3,iist2), idatec(4,iist2), idatec(5,iist2), idummy)

 print*, idatec(:,iist2)

    dz = h/nlev

    do jlat=1,nlat
    do jlon=1,nlon
    do jklev=1,nlev
    zita =(jklev-1)*dz+dz/2.
    zfz  = 1.-zita/h
    zeta(jlon,jlat,jklev) = max(phig(jlon,jlat)/g*gzita(zita) -h*bzita(zita)*log(zfz) &
                          - phig(jlon,jlat)/g, 0.) ! height above orography
    enddo
    enddo
    enddo

    tvirt = t * (1.+ep*q-qcw-qci)

    zz1 = -g*h*bzita(.5*dz)*log(1.-.5*dz/h)
    do jlat = 1, nlat
    do jlon = 1, nlon
    zdgz = phig(jlon,jlat)*(gzita(.5*dz)-1.) + zz1
    ps(jlon,jlat) = p(jlon,jlat,1)*exp(zdgz/(rd*tvirt(jlon,jlat,1)))
    enddo
    enddo

    call vtsurf (u10, v10, t2, q2, q2rel, hflux, qflux, ustar, nlon, nlat, nlev, fmask, rgm, rgq, &
                 richs, phig, ps, tskin, qskin, u, v, t, p, q, h, dz, b0)

    call cpbl (nlon, nlat, nlev, richs, zeta, u, v, t, p, q, qcw, qci, tvirt, pbl)

 nfilt = 3
 print*, nfilt
 do kf = 1, nfilt
 call filt (ps, nlon,nlat,1)
 call filt (t2, nlon,nlat,1)
 call filt (q2rel, nlon,nlat,1)
 call filt (ustar, nlon,nlat,1)
 call filt (pbl, nlon,nlat,1)
 call filt (t, nlon,nlat,nlev)
 call filt (p, nlon,nlat,nlev)
 call filt (q, nlon,nlat,nlev)
 call filt (qcw, nlon,nlat,nlev)
 call filt (qci, nlon,nlat,nlev)
 call filt (zeta, nlon,nlat,nlev)
 call filt (u, nlon,nlat,nlev)
 call filt (v, nlon,nlat,nlev)
 call filt (wh, nlon,nlat,nlev+1)
 enddo

  qcw = qcw + qci   ! to be improved....

! Vertical interpolation of w to integer levels from half-levels
! (interpolation along zita coordinate is simple linear combination)

  do k = 1, nlev
  do jlat = 1, nlat
  do jlon = 1, nlon
    w(jlon,jlat,k) = (wh(jlon,jlat,k+1) + wh(jlon,jlat,k))*0.5
  enddo
  enddo
  enddo

! Vertical interpolation of soil temperature to soil output level
! from input soil levels

  do jlat = 1, nlat
  do jlon = 1, nlon
    zvar_inp(1)         = tgsurf(jlon,jlat)
    zvar_inp(2:nlevg+1) = tg(jlon,jlat,1:nlevg)
    call interp (0.0, 1.0, 1.0, &
 nlevg+1, sqrt(soillev(1:nlevg+1)), zvar_inp(1:nlevg+1), &
 sqrt(soillev_out(1:nsoillev_out)), tg_levout(jlon,jlat,1:nsoillev_out), nsoillev_out)
  enddo
  enddo

! definizione di psfc
! tutti i flussi di moloch sono positivi se downward..

  call h_inter (ps, nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, psfc(:,:,iist2))
  work(:,:) = cswfl(:,:)/ztcu/(1.-albedo(:,:))
  call h_inter (work, nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, swrd(:,:,iist2))
  work(:,:) = clwfl(:,:)/ztcu + zstef*tskin(:,:)**4*emis1(:,:)  ! perche' chimere vuole quello downward...
  call h_inter (work, nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, lwrd(:,:,iist2))
  call h_inter (-chflux/ztcu, nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, sshf(:,:,iist2))
  call h_inter (-cqflux/ztcu, nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, slhf(:,:,iist2))
  call h_inter (qgsurf, nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, soim(:,:,iist2))

  call h_inter (snow, nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, weas(:,:,iist2))
  work(:,:) = prectot(:,:)/ztcu*3600.  ! total precip. in mm/h
  lspre(:,:) = min(work(:,:), 5./3600.) ! large scale precip.
  precconv(:,:)= work(:,:)-lspre(:,:) ! convective precip in mm/h
  call h_inter (lspre,  nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, lspc(:,:,iist2))
  call h_inter (precconv, nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, copc(:,:,iist2))

  call h_inter (ustar, nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, usta(:,:,iist2))
  call h_inter (q2rel, nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, rh2m(:,:,iist2))
  call h_inter (t2, nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, tem2(:,:,iist2))
  call h_inter (pbl, nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, hght(:,:,iist2))

  call h_inter (u10, nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, uf)
  call h_inter (v10, nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, vf)
  call anti_rot_wind (x0, y0, lon, lat, lonr, latr, nx, ny, uf, vf, u10m(:,:,iist2), v10m(:,:,iist2))

  do k = 1, nlev
    call h_inter (t(:,:,k), nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, temp(:,:,k,iist2))
    call h_inter (q(:,:,k), nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, sphu(:,:,k,iist2))
    call h_inter (.7*qcw(:,:,k), nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, cliq(:,:,k,iist2))
    call h_inter (.3*qcw(:,:,k), nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, rain(:,:,k,iist2)) ! parte della nube in pioggia
    call h_inter (p(:,:,k), nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, pres(:,:,k,iist2))
    call h_inter (zeta(:,:,k), nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, alti(:,:,k,iist2))
    call h_inter (w(:,:,k), nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, winv(:,:,k,iist2))

    call h_inter (u(:,:,k), nlon, nlat, zdlon, zdlat, zalon0+zdlon*.5, zalat0+zdlat*.5, lonr, latr, nx, ny, uf)
    call h_inter (v(:,:,k), nlon, nlat, zdlon, zdlat, zalon0, zalat0, lonr, latr, nx, ny, vf)
    call anti_rot_wind (x0, y0, lon, lat, lonr, latr, nx, ny, uf, vf, winz(:,:,k,iist2), winm(:,:,k,iist2))

  enddo

  do k = 1, nsoillev_out
    call h_inter (tg_levout(:,:,k), nlon, nlat, zdlon, zdlat, zalon0, zalat0+zdlat*.5, lonr, latr, nx, ny, soilt(:,:,k,iist2))
  enddo

  stl1(:,:,iist2) = soilt(:,:,1,iist2)

write (Times(iist2),'(I4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":",I2.2)') idatec(:,iist2), 0
print *, Times(iist2)

!if (iist2==2) then
!
!  z_lon_first = lon(1,1)
!  z_lat_first = lat(1,1)
!  dlon_c = lon(2,1) - lon(1,1)
!  dlat_c = lat(1,2) - lat(1,1)

! Surface --->
!!  call outgraph(80102,2,nlon,nlat,x0,y0,zalon0,zalat0+zdlat*0.5,zdlon,zdlat,&
!! 0, 0,  0,103,  2, 0,idate0(1:5),iperiod(1:3),t2(:,:),1.,-273.15)

!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 0,  0,103,  2, 0,idate0(1:5),iperiod(1:3),tem2(:,:,iist2),1.,-273.15)
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 1,  1,103,  2, 0,idate0(1:5),iperiod(1:3),rh2m(:,:,iist2),100.,0.)
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 2,  2,103, 10, 0,idate0(1:5),iperiod(1:3),u10m(:,:,iist2),1.,0.)
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 2,  3,103, 10, 0,idate0(1:5),iperiod(1:3),v10m(:,:,iist2),1.,0.)
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 4,  2,  1,  0, 0,idate0(1:5),iperiod(1:3),swrd(:,:,iist2),1.,0.)
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 5,  3,  1,  0, 0,idate0(1:5),iperiod(1:3),lwrd(:,:,iist2),1.,0.)
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 0, 11,  1,  0, 0,idate0(1:5),iperiod(1:3),sshf(:,:,iist2),1.,0.)
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 0, 10,  1,  0, 0,idate0(1:5),iperiod(1:3),slhf(:,:,iist2),1.,0.)
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 2, 0,  9,  1,  0, 0,idate0(1:5),iperiod(1:3),soim(:,:,iist2),1.,0.)
!
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 1, 13,  1,  0, 0,idate0(1:5),iperiod(1:3),weas(:,:,iist2),1.e3,0.)
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 1, 54,  1,  0, 0,idate0(1:5),iperiod(1:3),lspc(:,:,iist2),1.,0.)
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 1, 37,  1,  0, 0,idate0(1:5),iperiod(1:3),copc(:,:,iist2),1.,0.)
!
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 2, 30,  1,  0, 0,idate0(1:5),iperiod(1:3),usta(:,:,iist2),1.,0.)
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 3, 18,  1,  0, 0,idate0(1:5),iperiod(1:3),hght(:,:,iist2),1.,0.)
! <---
!
! Atm. levels --->
!  k=10
!  k=50
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 3,  6,105,  k, 0,idate0(1:5),iperiod(1:3),alti(:,:,k,iist2),1.,0.)
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 3,  0,105,  k, 0,idate0(1:5),iperiod(1:3),pres(:,:,k,iist2),1.e-2,0.)
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 2,  2,105,  k, 0,idate0(1:5),iperiod(1:3),winz(:,:,k,iist2),1.,0.)
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 2,  3,105,  k, 0,idate0(1:5),iperiod(1:3),winm(:,:,k,iist2),1.,0.)
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 0,  0,105,  k, 0,idate0(1:5),iperiod(1:3),temp(:,:,k,iist2),1.,-273.15)
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 1,  0,105,  k, 0,idate0(1:5),iperiod(1:3),sphu(:,:,k,iist2),1.e3,0.)
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 6,193,105,  k, 0,idate0(1:5),iperiod(1:3),cliq(:,:,k,iist2),1.e6,0.)
!  call outgraph(80102,2,nx,ny,0.,0.,z_lon_first,z_lat_first,dlon_c,dlat_c,&
! 0, 6,193,105,  k, 0,idate0(1:5),iperiod(1:3),rain(:,:,k,iist2),1.e6,0.)
! <--- 

!print *,minval(alti(:,:,50,iist2)),maxval(alti(:,:,50,iist2))
!print *,minval(pres(:,:,50,iist2)),maxval(pres(:,:,50,iist2))
!print *,minval(sphu(:,:,50,iist2)),maxval(sphu(:,:,50,iist2))
!print *,minval(cliq(:,:,50,iist2)),maxval(cliq(:,:,50,iist2))

!endif

  endif ! <--- mod(iist+iist0,njump) == 0

enddo ! while (.true.)


  !stop!!!!!!!!!!!!!!!!


 lwrd(:,:,1) = .5*(lwrd(:,:,1)+lwrd(:,:,2))

! inserire qui la conversione a netcdf

! Create output netCDF --> enter define mode
! ncstat=nf90_create(path=meteo_out, cmode=NF90_CLOBBER,ncid = ncout)  
 ncstat=nf90_create(path=meteo_out, cmode=IOR(NF90_NETCDF4, NF90_MPIIO),ncid = ncout)
 
 ! Definisco le dimensioni da usare per le variabili di moloch intepolato su chimere
 ! Time
 ncstat=nf90_def_dim(ncout,"Time",NF90_UNLIMITED,time_did)
 ! Date string lenght
 ncstat=nf90_def_dim(ncout,"DateStrLen",ndstr,dstr_did)
 ! Longitude
 ncstat=nf90_def_dim(ncout,"west_east",nx,nlon_did)
 ! Latitude
 ncstat=nf90_def_dim(ncout,"south_north",ny,nlat_did)
 ! Levels
 ncstat=nf90_def_dim(ncout,"bottom_top",nlev,lev_did)

 ! Scrive le varibili di griglia nel netcdf

 ! Longitude
 ncstat=nf90_def_var(                  &
      ncout,                           &
      'lon',                           &
      NF90_FLOAT,                      &
      (/nlon_did,nlat_did/),           &
      lon_vid)

 ncstat=nf90_put_att(                   &
       ncout,                           &
       lon_vid,                         &
       'units',                         &
       'degrees_east'                   &
       )

  ncstat=nf90_put_att(                  &
       ncout,                           &
       lon_vid,                         &
       'long_name',                     &
       'Longitude'                      &
       ) 

  ! Latitude 
  ncstat=nf90_def_var(                  &
       ncout,                           &
       'lat',                           &
       NF90_FLOAT,                      &
       (/nlon_did,nlat_did/),           &
       lat_vid)

  ncstat=nf90_put_att(                  &
       ncout,                           &
       lat_vid,                         &
       'units',                         &
       'degrees_north'                  &
       )
  ncstat=nf90_put_att(                  &
       ncout,                           &
       lat_vid,                         &
       'long_name',                     &
       'Latitude'                       &
       )


! definisco le variabili che sono funzione del tempo

 ! Times
  ncstat=nf90_def_var(                  &
       ncout,                           &
       'Times',                         &
       NF90_CHAR,                       &
       (/dstr_did,time_did/),           &
       time_vid) 

 ! 2D vars ------------------------------------------------
 ! psfc
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'psfc',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,time_did/),           &
         psfc_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       psfc_vid,                        &
       'units',                         &
       'Pa'                             &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      psfc_vid,                        &
      'long_name',                     &
      'Surface pressure'               &
       )

 ! swrd 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'swrd',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,time_did/),           &
         swrd_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       swrd_vid,                        &
       'units',                         &
       'W/m^2'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      swrd_vid,                        &
      'long_name',                     &
      'SW radiation down'               &
       )

 ! lwrd 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'lwrd',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,time_did/),           &
         lwrd_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       lwrd_vid,                        &
       'units',                         &
       'W/m^2'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      lwrd_vid,                        &
      'long_name',                     &
      'LW radiation down'               &
       )

 ! sshf 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'sshf',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,time_did/),           &
         sshf_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       sshf_vid,                        &
       'units',                         &
       'W/m^2'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      sshf_vid,                        &
      'long_name',                     &
      'Surface sensible heat flux'               &
       )

 ! slhf 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'slhf',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,time_did/),           &
         slhf_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       slhf_vid,                        &
       'units',                         &
       'W/m^2'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      slhf_vid,                        &
      'long_name',                     &
      'Surface latent heat flux'               &
       )

 ! usta 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'usta',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,time_did/),           &
         usta_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       usta_vid,                        &
       'units',                         &
       'm/s'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      usta_vid,                        &
      'long_name',                     &
      'Frictional velocity'            &
       )

 ! tem2 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'tem2',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,time_did/),           &
         tem2_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       tem2_vid,                        &
       'units',                         &
       'K'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      tem2_vid,                        &
      'long_name',                     &
      '2m air temperature'               &
       )

 ! soim 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'soim',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,time_did/),           &
         soim_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       soim_vid,                        &
       'units',                         &
       'm'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      soim_vid,                        &
      'long_name',                     &
      'Soil Moisture level 1'          &
       )

 ! stl1 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'stl1',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,time_did/),           &
         stl1_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       stl1_vid,                        &
       'units',                         &
       'K'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      stl1_vid,                        &
      'long_name',                     &
      'Soil temperature 0-7 cm'        &
       )


 ! hght 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'hght',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,time_did/),           &
         hght_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       hght_vid,                        &
       'units',                         &
       'm'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      hght_vid,                        &
      'long_name',                     &
      'PBL height from MOLOCH'         &
       )

 ! rh2m 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'rh2m',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,time_did/),           &
         rh2m_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       rh2m_vid,                        &
       'units',                         &
       'fraction'                       &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      rh2m_vid,                        &
      'long_name',                     &
      'Relative Humidity at 2m'               &
       )

 ! lspc 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'lspc',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,time_did/),           &
         lspc_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       lspc_vid,                        &
       'units',                         &
       'kg/m^2/h'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      lspc_vid,                        &
      'long_name',                     &
      'Large scale precipitation'               &
       )

 ! copc 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'copc',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,time_did/),           &
         copc_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       copc_vid,                        &
       'units',                         &
       'kg/m^2/h'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      copc_vid,                        &
      'long_name',                     &
      'convective precipitation'               &
       )

 ! u10m 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'u10m',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,time_did/),           &
         u10m_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       u10m_vid,                        &
       'units',                         &
       'm/s'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      u10m_vid,                        &
      'long_name',                     &
      '10 m U wind'               &
       )

 ! v10m 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'v10m',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,time_did/),           &
         v10m_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       v10m_vid,                        &
       'units',                         &
       'm/s'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      v10m_vid,                        &
      'long_name',                     &
      '10 m V wind'               &
       )

 ! weas 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'weas',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,time_did/),           &
         weas_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       weas_vid,                        &
       'units',                         &
       'mm'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      weas_vid,                        &
      'long_name',                     &
      'Water equiv. accum. snow  depth'               &
       )

 ! 3D vars ------------------------------------------------
 ! temp 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'temp',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,lev_did,time_did/),      &
         temp_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       temp_vid,                        &
       'units',                         &
       'K'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      temp_vid,                        &
      'long_name',                     &
      'Temperature'               &
       )

 ! cliq 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'cliq',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,lev_did,time_did/),      &
         cliq_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       cliq_vid,                        &
       'units',                         &
       'kg/kg'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      cliq_vid,                        &
      'long_name',                     &
      'Cloud liquid water mixing ratio'               &
       )

 ! rain 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'rain',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,lev_did,time_did/),      &
         rain_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       rain_vid,                        &
       'units',                         &
       'kg/kg'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      rain_vid,                        &
      'long_name',                     &
      'Rain water mixing ratio'               &
       )

 ! sphu 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'sphu',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,lev_did,time_did/),      &
         sphu_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       sphu_vid,                        &
       'units',                         &
       'kg/kg'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      sphu_vid,                        &
      'long_name',                     &
      'Specific humidity'               &
       )

 ! pres 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'pres',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,lev_did,time_did/),      &
         pres_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       pres_vid,                        &
       'units',                         &
       'Pa'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      pres_vid,                        &
      'long_name',                     &
      'Pressure'               &
       )

 ! alti 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'alti',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,lev_did,time_did/),      &
         alti_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       alti_vid,                        &
       'units',                         &
       'm'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      alti_vid,                        &
      'long_name',                     &
      'Altitude of half-sigma level'               &
       )

 ! winz 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'winz',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,lev_did,time_did/),      &
         winz_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       winz_vid,                        &
       'units',                         &
       'm/s'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      winz_vid,                        &
      'long_name',                     &
      'Zonal wind'               &
       )

 ! winm 
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'winm',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,lev_did,time_did/),      &
         winm_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       winm_vid,                        &
       'units',                         &
       'm/s'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      winm_vid,                        &
      'long_name',                     &
      'Meridional wind'               &
       )

 ! winw
 ncstat=nf90_def_var(                              &
         ncout,                                    &
         'winw',                                   &
         NF90_FLOAT,                               &
         (/nlon_did,nlat_did,lev_did,time_did/),      &
         winw_vid                                  &
         )
 ncstat=nf90_put_att(                   &
       ncout,                           &
       winw_vid,                        &
       'units',                         &
       'm/s'                          &
       )
 ncstat=nf90_put_att(                  &
      ncout,                           &
      winw_vid,                        &
      'long_name',                     &
      'Vertical wind'               &
       )


  ncstat=nf90_put_att(ncout,NF90_GLOBAL,"Title","MOLOCH4CHIMERE SUITE version 2022")
 
  ncstat=nf90_put_att(ncout,NF90_GLOBAL,"Sub-title","METEO fields for CHIMERE v2020")

  ncstat=nf90_put_att(ncout,NF90_GLOBAL,"Meteo info","vars list taken by CHIMERE exdomout file")

  ncstat=nf90_put_att(ncout,NF90_GLOBAL,"Generating_process","Generated by the program m2c.F90")

  ncstat=nf90_put_att(ncout,NF90_GLOBAL,"Domain","IT03M") ! auto

  ncstat=nf90_put_att(ncout,NF90_GLOBAL,"MOLOCH version","v22")
  ncstat=nf90_put_att(ncout,NF90_GLOBAL,"MOLOCH nlon", nlon) ! auto
  ncstat=nf90_put_att(ncout,NF90_GLOBAL,"MOLOCH nlat", nlat)! auto
  ncstat=nf90_put_att(ncout,NF90_GLOBAL,"MOLOCH nlev", nlev)! auto
  ncstat=nf90_put_att(ncout,NF90_GLOBAL,"MOLOCH rotation lon", x0)! auto
  ncstat=nf90_put_att(ncout,NF90_GLOBAL,"MOLOCH rotation lat", y0)! auto

  ncstat=nf90_put_att(ncout,NF90_GLOBAL,"Authors","T. C. Landi, O. Drofa and P. Malguzzi")

  ncstat=nf90_put_att(ncout,NF90_GLOBAL,"Institution","CNR - ISAC")
 
  ! Exit netCDF define mode
  ncstat=nf90_enddef(ncout) 
  ! Write vars
  ncstat=nf90_put_var(ncout,lon_vid,lon)
  ncstat=nf90_put_var(ncout,lat_vid,lat)
 
 ! scrivo le vars che dipendono dal tempo in netcdf

 ! Times
 ncstat=nf90_put_var(                          &
       ncout,                                  &
       time_vid,			       &
       Times(1:iist-1)                         &
       )
       
 ! 2D vars
 ! psfc 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      psfc_vid,                               &
      psfc(:,:,1:iist-1))

 ! swrd 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      swrd_vid,                               &
      swrd(:,:,1:iist-1))

 ! lwrd 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      lwrd_vid,                               &
      lwrd(:,:,1:iist-1))

 ! sshf 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      sshf_vid,                               &
      sshf(:,:,1:iist-1))

 ! slhf 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      slhf_vid,                               &
      slhf(:,:,1:iist-1))

 ! usta 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      usta_vid,                               &
      usta(:,:,1:iist-1))

 ! tem2 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      tem2_vid,                               &
      tem2(:,:,1:iist-1))

 ! soim 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      soim_vid,                               &
      soim(:,:,1:iist-1))

 ! stl1
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      stl1_vid,                               &
      stl1(:,:,1:iist-1))
 
 ! hght 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      hght_vid,                               &
      hght(:,:,1:iist-1))

 ! rh2m 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      rh2m_vid,                               &
      rh2m(:,:,1:iist-1))

 ! lspc 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      lspc_vid,                               &
      lspc(:,:,1:iist-1))

 ! copc 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      copc_vid,                               &
      copc(:,:,1:iist-1))

 ! u10m 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      u10m_vid,                               &
      u10m(:,:,1:iist-1))

 ! v10m 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      v10m_vid,                               &
      v10m(:,:,1:iist-1))

 ! weas 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      weas_vid,                               &
      weas(:,:,1:iist-1))


 ! 3D vars
 ! temp 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      temp_vid,                               &
      temp(:,:,:,1:iist-1))

! cliq 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      cliq_vid,                               &
      cliq(:,:,:,1:iist-1))

! rain 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      rain_vid,                               &
      rain(:,:,:,1:iist-1))

! sphu 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      sphu_vid,                               &
      sphu(:,:,:,1:iist-1))

! pres 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      pres_vid,                               &
      pres(:,:,:,1:iist-1))

! alti 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      alti_vid,                               &
      alti(:,:,:,1:iist-1))

! winz 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      winz_vid,                               &
      winz(:,:,:,1:iist-1))

! winm 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      winm_vid,                               &
      winm(:,:,:,1:iist-1))

! winw 
 ncstat=nf90_put_var(                         &
      ncout,                                  &
      winw_vid,                               &
      winv(:,:,:,1:iist-1))

 time = iist

 print*, "ora iist vale:", iist 
 
 ! Close output netCDF
 ncstat=nf90_close(ncout) 
 
end program m2c 

! ### Fine programma principale ###

! ### Sezione delle subroutines ###
!###############################################################################################################
      subroutine rdmhf_atm(kunit, nlon, nlat, nlev, nfdr, pdr, p, u, v, w, t, q, qcw, qci)

!   Read Model History File with atmospheric variables from kunit

      use mod_moloch, only : ierr_read1

      implicit none

      integer                           :: kunit, nlon, nlat, nlev, nlevg, jlon, jlat, jklev
      integer, dimension(50)            :: nfdr
      real, dimension(100)              :: pdr
      real, dimension(nlon,nlat,nlev)   :: p, u, v, t, q, qcw, qci
      real, dimension(nlon,nlat,nlev+1) :: w

      read(kunit, iostat=ierr_read1) nfdr

      if (ierr_read1 /= 0) then
        close(kunit)
        print*,'EOF reached on file unit ', kunit
        return
      endif

      read(kunit) pdr

      do jklev=1,nlev
      do jlat=1,nlat
        read (kunit) (p(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo
      do jklev=1,nlev
      do jlat=1,nlat
        read(kunit) (u(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo
      do jklev=1,nlev
      do jlat=1,nlat
        read(kunit) (v(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo
      do jklev=1,nlev+1
      do jlat=1,nlat
        read(kunit) (w(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo
      do jklev=1,nlev
      do jlat=1,nlat
        read(kunit) (t(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo
      do jklev=1,nlev
      do jlat=1,nlat
        read(kunit) (q(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo
      do jklev=1,nlev
      do jlat=1,nlat
        read(kunit) (qcw(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo
      do jklev=1,nlev
      do jlat=1,nlat
        read(kunit) (qci(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      return
      end
!###############################################################################################################
      subroutine rdmhf_soil(kunit, nlon, nlat, nlevg, nlev_snow, nfdr, pdr, &
                     rgm, rgq, fice, iceth, albedo, emis1, emis2, &
                     cloudt, prectot, precconvr, precsolidr, &
                     tskin, tgsurf, tg, qskin, qgsurf, qg, fice_soil_surf, fice_soil, &
                     snow, snow_lev, snow_t, snow_fice, snow_age, snow_melt_age, snow_dens, snow_albedo, &
                     cswflr, clwflr, chfluxr, cqfluxr, t2minr, t2maxr, ws10maxr, runoffr, runoff_tot_r, &
                     cwvfluxr, csoilh_bott_fluxr, csoilw_bott_fluxr)

!   Read Model History File with surface/sea/soil variables from kunit

      use mod_moloch, only : ierr_read2

      implicit none

      integer                           :: kunit, nlon, nlat, nlevg, nlev_snow,  jlon, jlat, jklev, ird
      integer, dimension(50)            :: nfdr, nfdr_local
      real, dimension(100)              :: pdr, pdr_local
      real, dimension(nlon,nlat)        :: rgm, rgq, fice, iceth, albedo, emis1, emis2, &
                                           cloudt, prectot, precconvr, precsolidr, &
                                           tskin, tgsurf, qskin, qgsurf, fice_soil_surf, &
                                           snow, snow_albedo, cswflr, clwflr, chfluxr, cqfluxr, t2minr, t2maxr, ws10maxr, &
                                           runoffr, runoff_tot_r, cwvfluxr, csoilh_bott_fluxr, csoilw_bott_fluxr
      real, dimension(nlon,nlat,nlevg)  :: tg, qg, fice_soil
      real, dimension(nlon,nlat,nlev_snow) :: snow_lev, snow_t, snow_fice, snow_age, snow_melt_age, snow_dens
      real, dimension(nlon,nlat)        :: field2d_add

      read(kunit, iostat=ierr_read2) nfdr_local

      if (ierr_read2 /= 0) then
        close(kunit)
        print*,'EOF reached on file unit ', kunit
        return
      endif

      read(kunit) pdr_local

      if (any(nfdr_local(:) /= nfdr(:)).or.any(pdr_local(:) /= pdr(:))) then
          print *, ' Header parameter (nfdr, pdr) read from unit ',kunit, &
 ' not coincide with defined parameters'
          print *,"   STOP"
          stop
      endif

!  Physiographical parameters changing in time

      do jlat=1,nlat
        read (kunit) (field2d_add(jlon,jlat),jlon=1,nlon) ! lai
      enddo

      do jlat=1,nlat
        read (kunit) (field2d_add(jlon,jlat),jlon=1,nlon) ! fveg
      enddo

      do jlat=1,nlat
        read (kunit) (rgm(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read (kunit) (rgq(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read (kunit) (iceth(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read (kunit) (fice(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read (kunit) (albedo(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read (kunit) (emis1(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read (kunit) (emis2(jlon,jlat),jlon=1,nlon)
      enddo

! Prognostic cloud and precipitation variables at the surface

      do jlat=1,nlat
        read(kunit) (cloudt(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (prectot(jlon,jlat),jlon=1,nlon)
      enddo

!      do jlat=1,nlat
!        read(kunit) (precconvr(jlon,jlat),jlon=1,nlon)
!      enddo

      do jlat=1,nlat
        read(kunit) (precsolidr(jlon,jlat),jlon=1,nlon)
      enddo

! Prognostic surface and soil/sea fields

      do jlat=1,nlat
        read(kunit) (tskin(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read (kunit) (tgsurf(jlon,jlat),jlon=1,nlon)
      enddo

      do jklev=1,nlevg
      do jlat=1,nlat
        read(kunit) (tg(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      do jlat=1,nlat
        read(kunit) (qskin(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read (kunit) (qgsurf(jlon,jlat),jlon=1,nlon)
      enddo

      do jklev=1,nlevg
      do jlat=1,nlat
        read(kunit) (qg(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      do jlat=1,nlat
        read (kunit) (fice_soil_surf(jlon,jlat),jlon=1,nlon)
      enddo

      do jklev=1,nlevg
      do jlat=1,nlat
        read(kunit) (fice_soil(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      do jlat=1,nlat
        read(kunit) (snow(jlon,jlat),jlon=1,nlon)
      enddo

      do jklev=1,nlev_snow
      do jlat=1,nlat
        read(kunit) (snow_lev(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      do jklev=1,nlev_snow
      do jlat=1,nlat
        read(kunit) (snow_t(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      do jklev=1,nlev_snow
      do jlat=1,nlat
        read(kunit) (snow_fice(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      do jklev=1,nlev_snow
      do jlat=1,nlat
        read(kunit) (snow_age(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      do jklev=1,nlev_snow
      do jlat=1,nlat
        read(kunit) (snow_melt_age(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      do jklev=1,nlev_snow
      do jlat=1,nlat
        read(kunit) (snow_dens(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      do jlat=1,nlat
        read(kunit) (snow_albedo(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (cswflr(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (clwflr(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (chfluxr(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (cqfluxr(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (t2minr(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (t2maxr(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (ws10maxr(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (runoffr(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (runoff_tot_r(jlon,jlat),jlon=1,nlon)
      enddo

! Reading of additional 2d fields

     do ird=1,3
      do jlat=1,nlat
      read(kunit) (field2d_add(jlon,jlat),jlon=1,nlon)
      enddo
     enddo

!      do jlat=1,nlat
!        read(kunit) (cwvfluxr(jlon,jlat),jlon=1,nlon)
!      enddo
!
!      do jlat=1,nlat
!        read(kunit) (csoilh_bott_fluxr(jlon,jlat),jlon=1,nlon)
!      enddo
!
!      do jlat=1,nlat
!        read(kunit) (csoilw_bott_fluxr(jlon,jlat),jlon=1,nlon)
!      enddo

      return
      end
!###############################################################################################################
      subroutine rd_param_const(nlon, nlat, nlevg, x0, y0, alon0, alat0, dlon, dlat, phig, fmask, qgmax, qgmin)

! Reads from additional input file
! all constant (in time) model physiographical parameters

      implicit none

      integer                           :: nlon, nlat, nlev, nlevg, &
                                           nlon_local, nlat_local, nlev_local, nlevg_local, &
                                           iunit=11, ierr_open, ierr=0, nst_local, nvt_local, &
 jlon, jlat, ird, jklev
      real                              :: x0, y0, alon0, alat0, dlon, dlat, &
                                           x0_local, y0_local, alon0_local, alat0_local, dlon_local, dlat_local
      real, parameter                   :: g=9.807
      real, dimension(nlon,nlat)        :: phig, fmask, zread
      real, dimension(nlon,nlat,nlevg)  :: qgmax, qgmin
      character(len=30) :: filerd="model_param_constant.bin"

     open (iunit,file=trim(filerd),status='old',form='unformatted',iostat=ierr_open)
     if (ierr_open /= 0) then
        print *
        print *,'Not found input ',trim(filerd)
        print *,"   stop"
        stop
      endif

      read (iunit) nlon_local, nlat_local, nlevg_local, dlon_local, dlat_local, x0_local, y0_local, alon0_local, alat0_local, &
 nst_local, nvt_local

      if (nlon_local /= nlon) ierr=ierr+1
      if (nlat_local /= nlat) ierr=ierr+1
      if (nlevg_local /= nlevg) ierr=ierr+1
      if (dlon_local /= dlon) ierr=ierr+1
      if (dlat_local /= dlat) ierr=ierr+1
      if (x0_local /= x0) ierr=ierr+1
      if (y0_local /= y0) ierr=ierr+1
      if (alon0_local /= alon0) ierr=ierr+1
      if (alat0_local /= alat0) ierr=ierr+1

      if (ierr /= 0) then
        print *,"Error in header parameters in input file ,",trim(filerd),", not coincident with defined parameters"
        print *,"Model nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0
        print *,"Read nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nlon_local, nlat_local, nlevg_local, dlon_local, dlat_local, x0_local, y0_local, alon0_local, alat0_local
        print *,"   stop"
        stop
      endif

!  land-sea fraction and orography

      do jlat=1,nlat
        read(iunit) (fmask(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(iunit) (phig(jlon,jlat),jlon=1,nlon)
      enddo
      phig(:,:) = phig(:,:)*g

      do jlat=1,nlat
        read(iunit) (zread(jlon,jlat),jlon=1,nlon) ! htopvar
      enddo
      do ird=1,nst_local+1
        do jlat=1,nlat
          read(iunit) (zread(jlon,jlat),jlon=1,nlon) ! soil_map
        enddo
      enddo
      do ird=1,nvt_local+1
        do jlat=1,nlat
          read(iunit) (zread(jlon,jlat),jlon=1,nlon) ! veg_map
        enddo
      enddo

      do jklev=1,nlevg
        do jlat=1,nlat
          read(iunit) (qgmax(jlon,jlat,jklev),jlon=1,nlon)
        enddo
      enddo

      do jklev=1,nlevg
        do jlat=1,nlat
          read(iunit) (qgmin(jlon,jlat,jklev),jlon=1,nlon)
        enddo
      enddo

      close(iunit)

      return
      end
!#######################################################################
 subroutine rrec2 (kunit, nlon, nlat, vect)
 real vect (nlon,nlat)

 do jlat = 1, nlat
 read (kunit) (vect(jlon,jlat), jlon = 1, nlon)
 enddo

 return
 end
!#######################################################################
!#######################################################################
!#######################################################################
 subroutine rot_grid(x0,y0,x00,y00,dlon,dlat,alon,alat,nlon,nlat)
 
 ! Calculation of standard geographical coordinates (ALON, ALAT) of model
 ! (rotated) grid points (rotation centre: X0, Y0)
 ! X00, Y00 are the SW corner point coordinates in rotated coordinates
 ! Some computations require double precision
 
 implicit none
 
 real :: x0, y0, x00, y00, dlon, dlat
 integer :: nlon, nlat, jlat, jlon
 real, dimension(nlon,nlat) :: alat, alon
 real*8 :: zfac, zx0, zy0, zlon, zlat, zzlat, zaarg, zarg
 
 if (abs(x0)>0.01.or.abs(y0)>0.01) then
 
 ! Case of rotated grid

 zfac = dabs(dacos(-1.d0))/180.d0
 zx0  = dble(x0)*zfac
 zy0  = dble(y0)*zfac
 do jlat=1,nlat
   zlat=(dble(y00)+(jlat-1)*dble(dlat))*zfac
   do jlon=1,nlon
     zlon = (dble(x00)+(jlon-1)*dble(dlon))*zfac
     zzlat= 1.d0/zfac*dasin( dcos(zy0)*dsin(zlat) +  &
                      dsin(zy0)*dcos(zlat)*dcos(zlon) )
     zarg = -dsin(zlat)*dsin(zy0)+dcos(zy0)*dcos(zlat)*dcos(zlon)
     zaarg = zarg/dcos(zfac*zzlat)
     alat(jlon,jlat)=zzlat
     if(zaarg < -1.d0.and.zaarg > -1.00001d0) zaarg = -1.d0
     if(zaarg >  1.d0.and.zaarg <  1.00001d0) zaarg =  1.d0

     if (abs(abs(alat(jlon,jlat))-90.) > 1.e-7) then
       if (zlon < 0.d0) then
         alon(jlon,jlat) = 1.d0/zfac*(zx0-dacos(zaarg))
       else
         alon(jlon,jlat) = 1.d0/zfac*(zx0+dacos(zaarg))
       endif
     else
       alon(jlon,jlat) = 0.
     endif

     if (alon(jlon,jlat) >  180.) alon(jlon,jlat) = alon(jlon,jlat) - 360.
     if (alon(jlon,jlat) < -180.) alon(jlon,jlat) = alon(jlon,jlat) + 360.
   enddo
 enddo

else

! Case of non rotated grid

 do jlon=1,nlon
   alon(jlon,:)=x00+(jlon-1)*dlon
 enddo
 do jlat=1,nlat
   alat(:,jlat)=y00+(jlat-1)*dlat
 enddo
 
endif

return
end subroutine rot_grid
!=======================================================================
subroutine anti_rot_grid(x0,y0,alon,alat,xr,yr,nlon,nlat)

! Calculation of the rotated coordinates (put in XR, YR) with rotation centre X0, Y0
! (given in geographical coordinates) for input grid points given in
! geographical coordinates (defined in ALON, ALAT)
! Some computations require double precision

implicit none

real :: x0, y0
integer :: nlon, nlat, jlat, jlon
real, dimension(nlon,nlat) :: alat, alon, xr, yr
real*8 :: zfac, zx0, zy0, zx, zy, zlon, zlat, zz, pi

if (abs(x0)>0.01.or.abs(y0)>0.01) then

! Case of rotated grid

  pi  = dabs(dacos(-1.d0))
  zfac = pi/180.d0
  zx0  = dble(x0)*zfac
  zy0  = dble(y0)*zfac
  do jlat=1,nlat
  do jlon=1,nlon
    if (jlon==1.and.jlat==1) call sleep (0) ! fix for a problem of ifort 14 (but needs compil. with -O2 or less)
    zx = dble(alon(jlon,jlat))*zfac
    if (zx-zx0 > pi) zx = zx - 2.d0*pi
    if (zx-zx0 < -pi) zx = zx + 2.d0*pi
    zy = dble(alat(jlon,jlat))*zfac
    zlat = dasin( -dcos(zy)*dsin(zy0)*dcos(zx-zx0)+dsin(zy)*dcos(zy0) )
    zz = dsin(zy)-dcos(zy0)*dsin(zlat)
    zz = zz/(dsin(zy0)*dcos(zlat))
    if (zz < -1.d0.and.zz > -1.00001d0) zz = -1.d0
    if (zz > 1.d0.and.zz < 1.00001d0) zz = 1.d0
    if (zx < zx0) then
      zlon = -dacos(zz)
    else
      zlon = dacos(zz)
    endif
    zx = zlon/zfac
    zy = zlat/zfac
    xr(jlon,jlat) = sngl(zx)
    yr(jlon,jlat) = sngl(zy)
  enddo
  enddo

else

! Case of non rotated grid

  xr(:,:) = alon(:,:)
  yr(:,:) = alat(:,:)

endif

return
end subroutine anti_rot_grid

   subroutine h_inter (a, nx, ny, dx, dy, x00, y00, xo, yo, nxo, nyo, f)

!  Horizontal bilinear interpolation between 2-D grids
!  Version to interpolate from a regular grid in input
!  with prescribed values defined in A(NX,NY)
!  Output values defined in vector F(nxo,nyo).

   real a(nx,ny), xo(nxo,nyo), yo(nxo,nyo), f(nxo,nyo)

   dxr = 1./dx
   dyr = 1./dy

   do jo = 1, nyo
   do io = 1, nxo
   i = int((xo(io,jo)-x00)*dxr+1.0)
   j = int((yo(io,jo)-y00)*dyr+1.0)
   xi = x00+(i-1)*dx
   yj = y00+(j-1)*dy
   ip1 = min(i+1,nx)
   jp1 = min(j+1,ny)
   x = (xo(io,jo)-xi)*dxr
   y = (yo(io,jo)-yj)*dyr
   f1 = a(i,j  ) + x*(a(ip1,j  )-a(i,j  ))
   f2 = a(i,jp1) + x*(a(ip1,jp1)-a(i,jp1))
   f(io,jo) = f1 + y*(f2-f1)
   enddo
   enddo

   return
   end subroutine h_inter
!=======================================================================
subroutine interp(alfa, ex1, ex2, npi, xi, g, x, f, nval)

! (Similar to INTERP_SPLINE_1D, but with different input, and does not require
!  call of SUBR. NEAR - moreover, it checks that input coordinates are monotonic).

!  Interpolates with splines with tension in one dimension.
!  The spline is defined imposing that the second derivative is the average
!  of second derivatives computed at the two adjacent points.
!  At interval extremes the second derivative is assumed null.
!  This subroutine also extrapolates out of the interval where the input funtion G is defined

!  INPUT:  function G defined at coordinates XI (CAUTION: can be changed by this subroutine)
!          G(1:NPI) values at irregular but strictly growing coordinates XI(1:NPI)
!  OUTPUT: F(1:NVAL) interpolated values at arbitrary coordinates X(1:NVAL)

!  ALFA: spline tension parameter, comprised between 0 and 1
!  If ALFA=1, pure linear interpolation; if ALFA=0, pure spline

!  EX1: param. determining extrapolation for X < XI(1)
!  EX2: param. determining extrapolation for X > XI(NPI)
!  If EX1=0 OR EX2=0, constant value extrapolation is used at corresponding extreme
!  If EX1=1 OR EX2=1, linear extrapolation is used at corresponding extreme
!  Intermediate values of EX1 AND EX2 give intermediate extrapolation values

  real, dimension(npi )  :: xi, g
  real, dimension(nval)  :: x,  f

  if(alfa.lt..0.or.alfa.gt.1) then
  print*, 'CAUTION: in INTERP, ALFA out of interval 0-1'
  endif
  if(ex1.lt..0.or.ex1.gt.1) then
  print*, 'CAUTION: in INTERP, EX1 out of interval 0-1'
  endif
  if(ex2.lt..0.or.ex2.gt.1) then
  print*, 'CAUTION: in INTERP, EX2 out of interval 0-1'
  endif

! Fix for the case in which coordinates of the input function are not strictly increasing
! Note that this changes the input coordinates in the calling programme

  do k=2,npi
   if(xi(k).le.xi(k-1)) then
   print*, "CAUTION: in INTERP, coordinates of input function changed because not monotonic!"
   exit
   endif
  enddo

  zeps=(xi(npi)-xi(1))*1.e-6   ! Small deviation used to set apart interlaced coordinates
200 do k=2,npi
     if(xi(k).le.xi(k-1)) then
     ximed=0.5*(xi(k)+xi(k-1))
     xi(k-1)=ximed-zeps
     xi(k)=ximed+zeps
     gmed=0.5*(g(k)+g(k-1))
     g(k-1)=gmed
     g(k)=gmed
     endif
    enddo

 do k=2,npi
  if(xi(k).le.xi(k-1)) then
  goto 200
  endif
 enddo

 do 100 jval =1, nval

!  2 cases of extrapolation

 if(x(jval).lt.xi(1)) then
 f(jval) = g(1) + ex1*(g(1)-g(2))/(xi(1)-xi(2)) * (x(jval)-xi(1))
 go to 100
 elseif (x(jval).ge.xi(npi)) then
 f(jval) = g(npi) + ex2*(g(npi)-g(npi-1))/(xi(npi)-xi(npi-1)) * (x(jval)-xi(npi))
 go to 100
 endif

 ir = 0

!  IR is a reference index determining the interpolation interval
!  The interpolation expression is applied also if X=XI(J)

 do j = 1, npi
 if (x(jval).ge.xi(j)) ir = ir + 1
 enddo

 if (ir.eq.1) then
 fmm = 2*g(1) - g(2)
 xmm = 2*xi(1) - xi(2)
 fpp = g(ir+2)
 xpp = xi(ir+2)
 elseif (ir.eq.(npi-1)) then
 fpp = 2*g(npi) - g(npi-1)
 xpp = 2*xi(npi) - xi(npi-1)
 fmm = g(ir-1)
 xmm = xi(ir-1)
 else
 fmm = g(ir-1)
 xmm = xi(ir-1)
 fpp = g(ir+2)
 xpp = xi(ir+2)
 endif

 fm     = g(ir)
 xm     = xi(ir)
 fp     = g(ir+1)
 xp     = xi(ir+1)
 delx   = xp - xm
 delxp  = xpp - xp
 delxm  = xm - xmm
 delx1  = x(jval) - xm
 delx2  = xp - x(jval)
 delxs  = delx**2
 delx1s = delx1**2
 delx2s = delx2**2

!  Spline contribution to interpolation

 spl = fm*(delx2/delx + delx1*delx2s/(delxs*delxm) - delx1s*     &
       delx2/((delx+delxp)*delxs)) + fp*(delx1/delx +            &
       delx1s*delx2/(delxs*delxp) - delx1*delx2s/((delx+delxm)*  &
       delxs)) - fmm * delx1*delx2s/((delx+delxm)*delx*delxm) -  &
       fpp * delx1s*delx2/((delx+delxp)*delx*delxp)

!  Linear interpolation contribution

 clin = (fm*delx2 + fp*delx1)/delx

!  Final interpolation combined using ALFA

 f(jval) = alfa*clin + (1.-alfa)*spl

 100  continue

 return
end subroutine interp
!=======================================================================

      subroutine calendar (nyrin, nmonin, ndayin, nhouin, nminin, iday,        &
                           ihou, imin, nyrc, nmonc, ndayc, nhouc, nminc, ndayr)

!  Defines current calendar date (nyrc, nmonc, ndayc, nhouc, nminc) and
!  day of the year (ndayr) by adding the forecast time (iday, ihou, imin)
!  to the initial date (nyrin, nmonin, ndayin, nhouin, nminin)

      implicit none
      integer imonth(12), nyrin, nmonin, ndayin, nhouin, nminin, iday, ihou, imin, &
               nyrc, nmonc, ndayc, nhouc, nminc, ndayr, j, itday, isum

!  define day of the year ( 1 < ndayr < 366 ) of initial date

      imonth( 1) = 31
      imonth( 2) = 28
      imonth( 3) = 31
      imonth( 4) = 30
      imonth( 5) = 31
      imonth( 6) = 30
      imonth( 7) = 31
      imonth( 8) = 31
      imonth( 9) = 30
      imonth(10) = 31
      imonth(11) = 30
      imonth(12) = 31
      if (mod(nyrin,4).eq.0) imonth(2) = 29

      ndayr = ndayin
      do j = 1, nmonin-1
      ndayr = ndayr + imonth(j)
      enddo

!  update initial hour and minutes

      nminc = nminin + imin
      nhouc = nhouin + ihou
      ndayr = ndayr  + iday
      if (nminc.ge.60) then
      nminc = nminc - 60
      nhouc = nhouc + 1
      endif
      if (nhouc.ge.24) then
      nhouc = nhouc - 24
      ndayr = ndayr + 1
      endif

!  update ndayr and initial day, month, year

      nyrc = nyrin
 1    itday  = 365
      if (mod(nyrc,4).eq.0) itday=366
      if (ndayr.gt.itday) then
      ndayr = ndayr - itday
      nyrc  = nyrc + 1
      else
      go to 2
      endif
      go to 1

 2    imonth(2)=28
      if (mod(nyrc,4).eq.0) imonth(2)=29
      isum=0
      do nmonc = 1, 12
      isum = isum + imonth(nmonc)
      if (ndayr.le.isum) go to 3
      enddo
 3    ndayc = ndayr + imonth(nmonc) - isum

      return
      end subroutine calendar
!#######################################################################
subroutine comp_esk(esat, qsat, t, p, iflag)

! Computes esat from temperature and qsat from absolute temperature and pressure
! IFLAG 1: esat and qsat with respect to water and ice, separately, depending if t>tzer or t<tzer
!       2: esat and qsat with an interpolation at t<tzer between water and ice
!       3: esat and qsat with respect to water also for t<tzer

  tzer = 273.15
  ezer = 611.
  cpv  = 1869.46
  cw   = 4186.8
  rd   = 287.05
  rv   = 461.51
  yliv = 2834170.5
  yliw = 333560.5
  eps  = rd/rv
  ci   = cw/2.
  ylwv = yliv-yliw
  ccw1 = (cpv-cw)/rv
  ccw2 = ylwv/tzer/rv-ccw1
  cci1 = (cpv-ci)/rv
  cci2 = yliv/tzer/rv-cci1

  zt0t = tzer/t
  if (zt0t.le.1.) then
   zesk = ezer*exp(-ccw1*alog(zt0t)+ccw2*(1.-zt0t))  ! Partial pressure over water
   else
    if (iflag.eq.1) then
    zesk = ezer*exp(-cci1*alog(zt0t)+cci2*(1.-zt0t)) ! Partial pressure over ice
    elseif (iflag.eq.2) then
    zratio = 1.04979*(0.5 + 0.5*tanh((t-tzer+9.)/6.))
    zesk = zratio*(ezer*exp(-ccw1*alog(zt0t)+ccw2*(1.-zt0t)))+  &
      (1.-zratio)*(ezer*exp(-cci1*alog(zt0t)+cci2*(1.-zt0t)))
    elseif (iflag.eq.3) then
    zesk = ezer*exp(-ccw1*alog(zt0t)+ccw2*(1.-zt0t))
    else
    print*, "IFLAG out of range in subroutine comp_esk", iflag
    stop
    endif
  endif

  esat=zesk
  qsat = zesk*eps/(p+zesk*(eps-1.))
  return
end subroutine comp_esk
!#######################################################################################################
  subroutine anti_rot_wind (x0, y0, lon, lat, lonr, latr, nlon, nlat, uf, vf, u, v)

! Input: x0, y0, lon, lat, lonr, latr, nlon, nlat, uf, vf
! Output: u, v
! Calculation of horizontal wind components in non rotated geography (u, v)
! for input wind components fields given in rotated geography (uf, vf)
! lon, lat: matrices containing the longitudes and latitudes of lat-lon coordinates
! Computations require double precision, but inaccuracies remain at poles and along the meridian x0

  implicit none
  real :: x0, y0, alon_rot, alat_rot
  integer :: nlon, nlat, jlon, jlat
  real, dimension(nlon,nlat) :: lon, lat, lonr, latr, u, v, uf, vf
  real*8 :: zfac, zx0, zy0, zx, zy, zxx, zyy, zu, zv, zzv, zmod, zphi, zpi

  if (abs(x0) > 0.01.or.abs(y0) > 0.01) then

! Case of rotated grid

  zpi  = dabs(dacos(-1.d0))
  zfac = zpi/180.d0
  zx0  = dble(x0)*zfac
  zy0  = dble(y0)*zfac

  do jlat = 1, nlat
  do jlon = 1, nlon

  alon_rot = lonr(jlon,jlat)
  alat_rot = latr(jlon,jlat)

  zx = dble(lon(jlon,jlat))*zfac
  if (zx-zx0.le.-zpi) zx = zx +2.d0*zpi
  if (zx-zx0.ge. zpi) zx = zx -2.d0*zpi
  zy  = dble(lat(jlon,jlat))*zfac

  zxx = dble(alon_rot)*zfac
  zyy = dble(alat_rot)*zfac
  zu  = dble(uf(jlon,jlat))
  zv  = dble(vf(jlon,jlat))

  if (dabs(dcos(zy))>1.d-2) then  ! more that .5 degrees apart from poles...
    if (dabs(dsin(zx-zx0))>1.d-3) then
    zzv=-dsin(zy0)*dsin(zxx)*zu/dcos(zy)+zv*(dcos(zy0)*dcos(zyy)-dsin(zy0)*dsin(zyy)*dcos(zxx))/dcos(zy)
    zu=dcos(zyy)/dsin(zx-zx0)/dsin(zy0)*(zv-zzv/dcos(zyy)*(dcos(zx-zx0)*dsin(zy0)*dsin(zy)+dcos(zy0)*dcos(zy)))
    zv=zzv
    u(jlon,jlat)=sngl(zu)
    v(jlon,jlat)=sngl(zv)
    else
       if (zx-zx0.lt.-zpi/2.d0.or.zx-zx0.gt.zpi/2.d0) then
       u(jlon,jlat) = -sngl(zu)
       v(jlon,jlat) = -sngl(zv)
       else
       u(jlon,jlat) = sngl(zu)
       v(jlon,jlat) = sngl(zv)
       endif
    endif
  else  ! close to poles...
  zmod = dsqrt(zu**2+zv**2)
  zphi = zpi/2.d0+datan2(zv,zu)
  if(zphi.gt.zpi) zphi = zphi-2.d0*zpi
  if(zphi.lt.zpi) zphi = zphi+2.d0*zpi
  if (zy0.lt.0.d0) zx = -zx             ! south pole
  u(jlon,jlat) = -zmod*dsin(zx-zphi)
  v(jlon,jlat) = -zmod*dcos(zx-zphi)
  endif

  enddo
  enddo

  else

! Case of non rotated grid

  u = uf
  v = vf

  endif

  return
  end subroutine anti_rot_wind
!#######################################################################
    subroutine filt (p, nlon, nlat, nlev)

    real p(nlon,nlat,nlev), p2(nlon,nlat)

!-----------------------
!  Horizontal diffusion
!-----------------------

    do 50 k = 1, nlev

    do jlat = 1, nlat
    do jlon = 2, nlon-1
    p2(jlon,jlat) = .25*(p(jlon-1,jlat,k)+p(jlon+1,jlat,k))+.5*p(jlon,jlat,k)
    enddo
    enddo

    do jlat = 2, nlat-1
    do jlon = 2, nlon-1
    p(jlon,jlat,k) = .25*(p2(jlon,jlat+1)+p2(jlon,jlat-1))+.5*p2(jlon,jlat)
    enddo
    enddo

50  continue

    return
    end subroutine filt
!###############################################################################################################
      subroutine vtsurf (u10, v10, t2, q2, q2rel, hflux, qflux, ustar, nlon, nlat, nlev, fmask, rgm, rgq, &
                         richs, phig, ps, tskin, qskin, u, v, t, p, q, h, dz, b0)

      use mod_moloch, only : gzita, bzita

!  Interpolates wind at 10 m, temperature and spec. humid. at 2 m.

      real, dimension(nlon,nlat,nlev) :: u, v, t, p, q
      real, dimension(nlon,nlat)      :: u10, v10, t2, td2, q2, q2rel, hflux, qflux, richs, &
                                         fmask, phig, rgm, rgq, ps, tskin, qskin, ustar

      real, parameter :: zaholt=1., zbholt=2./3., zcholt=5., zdholt=0.35,                                &
                         ak=.4, zbet=5., zgam=16., yliv=2834170.5, yliw=333560.5, ylwv=yliv-yliw,        &
                         tzer=273.15, pzer=1.e5, ezer= 611., rd=287.05, rv=461.51, eps=rd/rv,            &
                         cpd=1004.6, cvd=cpd-rd, cpv=1869.46, cw=4186.8, ci=cw/2., rdrcp=rd/cpd, g=9.807

!  Constants for computing saturation partial pressure

      real, parameter :: ccw1=(cpv-cw)/rv, ccw2=ylwv/tzer/rv-ccw1, cci1=(cpv-ci)/rv, cci2=yliv/tzer/rv-cci1

!  Businger functions when Ri<0

      psium(zx1,zx2)=alog((1.+zx1)**2*(1.+zx1**2)/((1.+zx2)**2*(1.+zx2**2)))-2.*(atan(zx1)-atan(zx2))
      psiuh(zy1,zy2)=2.*alog((1.+zy1)/(1.+zy2))

!  Holtslag functions when Ri>0

      psism(zx1,zx2)=-zaholt*zx1-zbholt*(zx1-zcholt/zdholt)*exp(-zdholt*zx1) &
                     +zaholt*zx2+zbholt*(zx2-zcholt/zdholt)*exp(-zdholt*zx2)
      psish(zy1,zy2)=-(1.+2./3.*zaholt*zy1)**1.5-zbholt*(zy1-zcholt/zdholt)*exp(-zdholt*zy1) &
                     +(1.+2./3.*zaholt*zy2)**1.5+zbholt*(zy2-zcholt/zdholt)*exp(-zdholt*zy2)

!  Turbulent fluxes at the ground (positive upward)

      zep=1./eps-1.

      do 100 jlat=1,nlat
      jlatp1=min(jlat+1,nlat)
      do 100 jlon=1,nlon
      jlonm1=max(jlon-1,1)

      zphi  = phig(jlon,jlat)*gzita(.5*dz)-g*h*bzita(.5*dz)*log(1.-.5*dz/h)
      za    = (zphi-phig(jlon,jlat))/g
      zua   = 0.5*(u(jlon,jlat,1)+u(jlonm1,jlat,1))
      zva   = 0.5*(v(jlon,jlatp1,1)+v(jlon,jlat,1))
      zmod2  = zua**2 + zva**2 + 0.07
      zmod  = sqrt(zmod2)
      zros  = ps(jlon,jlat)/( rd*tskin(jlon,jlat)*(1.+zep*qskin(jlon,jlat)) )

!  Virtual potential temperature computed with skin temperature

      zconvg = (pzer/ps(jlon,jlat) )**rdrcp*(1.+zep*qskin(jlon,jlat))
      ztevg  = tskin(jlon,jlat)*zconvg

      zconv1 = (pzer/p(jlon,jlat,1))**rdrcp*(1.+zep*q(jlon,jlat,1)  )
      ztetav = t(jlon,jlat,1)*zconv1

!  Richardson number above the ground

      ztebar = .5*(ztevg+ztetav)
      zri    = za*g*(ztetav-ztevg)/(ztebar*zmod2)
      richs(jlon,jlat) = zri ! to be used for the definition of clouds

!  Definition of roughness zrgm, zrgt, zrgq

      if(fmask(jlon,jlat).gt..5) then

!  Computation of Charnok roughness

      zchar = 5.e-4
      zcoch1=.0185*ak**2/g*zmod2
        do jiter = 1, 5
        zcoch2= za/zchar
        zchar = zcoch1/alog(zcoch2)**2
        enddo
      zrgm=zchar

        if(zri.ge.0.25) then
        zrgt = 2.2e-9
        zrgq = zrgt
        elseif(zri.lt.0.) then
        zrgt = 5.e-5
        zrgq = 9.e-5
        else
        zrgt = 3.80227e-6*exp(-zri*29.819387)
        zrgq = zrgt
        endif

        zrgmd=zrgm
        zrgtd=zrgt
        zrgqd=zrgq
      else
      zsea=5.e-5

      zrgm=fmask(jlon,jlat)*zsea+(1.-fmask(jlon,jlat))*rgm(jlon,jlat)
      zrgt=fmask(jlon,jlat)*zsea+(1.-fmask(jlon,jlat))*rgq(jlon,jlat)
      zrgq=zrgt

! Local roughness set to a few cm if larger (for interpolated
! diagnostics at 2 m (and other near surface levels) and 10 m only,
! not for computing fluxes)

      zrgmd=fmask(jlon,jlat)*zsea+(1.-fmask(jlon,jlat))*min(rgm(jlon,jlat), .05)
      zrgtd=fmask(jlon,jlat)*zsea+(1.-fmask(jlon,jlat))*min(rgq(jlon,jlat), .05)
      zrgqd=zrgtd
      endif

!  Computation of za/l

      zalzam=alog((za+zrgm)/zrgm)
      zalzat=alog((za+zrgt)/zrgt)
      zalzaq=alog((za+zrgq)/zrgq)

      zalzamd=alog((za+zrgmd)/zrgmd)
      zalzatd=alog((za+zrgtd)/zrgtd)
      zalzaqd=alog((za+zrgqd)/zrgqd)

      if(zri.ge.0.) then      ! Holtslag functions

      zpsim=0.
      zpsimd=0.
      zpsih=0.
      zpsihd=0.

      do jiter=1,4
      zal=zri*(zalzam-zpsim)**2/(zalzat-zpsih)
      zald=zri*(zalzamd-zpsimd)**2/(zalzatd-zpsihd)
      zx1=zal*(za+zrgm)/za
      zx1d=zald*(za+zrgmd)/za
      zx2=zal*zrgm/za
      zx2d=zald*zrgmd/za
      zy1=zal*(za+zrgt)/za
      zy1d=zald*(za+zrgtd)/za
      zy2=zal*zrgt/za
      zy2d=zald*zrgtd/za
      zpsim=psism(zx1,zx2)
      zpsimd=psism(zx1d,zx2d)
      zpsih=psish(zy1,zy2)
      zpsihd=psish(zy1d,zy2d)
      enddo
      zphim=1.+zaholt*zal+zbholt*zal*(1.+zcholt-zdholt*zal)*exp(-zdholt*zal)
      zphimd=1.+zaholt*zald+zbholt*zald*(1.+zcholt-zdholt*zald)*exp(-zdholt*zald)
      zpsiq=zpsih
      zpsiqd  = zpsihd
      zpsim10= zpsimd*(10./za)
      zpsih2 = zpsihd*(2. /za)
      zpsiq2 = zpsih2

      else                    ! Businger functions

      zpsim=0.
      zpsimd=0.
      zpsih=0.
      zpsihd=0.
      do jiter=1,4
      zal=zri*(zalzam-zpsim)**2/(zalzat-zpsih)
      zald=zri*(zalzamd-zpsimd)**2/(zalzatd-zpsihd)
      zx1=(1.-zgam*zal*(za+zrgm)/za)**.25
      zx1d=(1.-zgam*zald*(za+zrgmd)/za)**.25
      zx2=(1.-zgam*zal*zrgm/za)**.25
      zx2d=(1.-zgam*zald*zrgmd/za)**.25
      zy1=sqrt(1.-zgam*zal*(za+zrgt)/za)
      zy1d=sqrt(1.-zgam*zald*(za+zrgtd)/za)
      zy2=sqrt(1.-zgam*zal*zrgt/za)
      zy2d=sqrt(1.-zgam*zald*zrgtd/za)
      zpsim=psium(zx1,zx2)
      zpsimd=psium(zx1d,zx2d)
      zpsih=psiuh(zy1,zy2)
      zpsihd=psiuh(zy1d,zy2d)
      enddo
      zz1=sqrt(1.-zgam*zal*(za+zrgq)/za)
      zz1d=sqrt(1.-zgam*zald*(za+zrgqd)/za)
      zz2=sqrt(1.-zgam*zal*zrgq/za)
      zz2d=sqrt(1.-zgam*zald*zrgqd/za)
      zpsiq=psiuh(zz1,zz2)
      zpsiqd=psiuh(zz1d,zz2d)
      zx1d=(1.-zgam*(10.+zrgmd)/za*zald)**0.25
      zy1d=sqrt(1.-zgam*zald*(2.+zrgtd)/za)
      zz1d=sqrt(1.-zgam*zald*(2.+zrgqd)/za)
      zpsim10=psium(zx1d,zx2d)
      zpsih2 =psiuh(zy1d,zy2d)
      zpsiq2 =psiuh(zz1d,zz2d)

      endif

      zustar = ak*zmod/(zalzam-zpsim)
      ztstar = ak*(ztetav-ztevg)/(zalzat-zpsih)
      zqstar = ak*(q(jlon,jlat,1)-qskin(jlon,jlat))/(zalzaq-zpsiq)

      ustar(jlon,jlat) = zustar

      zustard = ak*zmod/(zalzamd-zpsimd)
      ztstard = ak*(ztetav-ztevg)/(zalzatd-zpsihd)
      zqstard = ak*(q(jlon,jlat,1)-qskin(jlon,jlat))/(zalzaqd-zpsiqd)

!  Surface fluxes of sensible and latent heat (positive upward)

      zcdt = zustar*ak/(zalzat-zpsih)
      zperflut = -zros*zcdt*cpd/zconvg
      hflux(jlon,jlat)=zperflut*(ztetav-ztevg)

      zcdq = zustar*ak/(zalzaq-zpsiq)
      zlate=ylwv-2360.*(tskin(jlon,jlat)-tzer)
      zperfluq = -zros*zcdq*zlate
      qflux(jlon,jlat)=zperfluq*(q(jlon,jlat,1)-qskin(jlon,jlat))
!-----------------------------------------------------------------------
      zv10=zustard/ak*(alog((10.+zrgmd)/zrgmd)-zpsim10)
      zt2 =ztevg+ztstard/ak*(alog((2. +zrgtd)/zrgtd)-zpsih2)
      zq2 =qskin(jlon,jlat)+zqstard/ak*(alog((2. +zrgqd)/zrgqd)-zpsiq2)
!-----------------------------------------------------------------------

!  Wind components

      u10(jlon,jlat) = zua/zmod*zv10
      v10(jlon,jlat) = zva/zmod*zv10

!  From potential temperature to temperature at 2 m in Celsius

      zconv2 = (pzer/ps(jlon,jlat) )**rdrcp*(1.+zep*zq2)
      zt2=zt2/zconv2
      t2(jlon,jlat) = zt2
      q2(jlon,jlat) = zq2

! Specific humidity near surface must not exceed saturation with respect to water and ice
! Relative humidity computed with respect to water (WMO standard)
! Dew point temperature

      call comp_esk(zeskl,zqsat,t2(jlon,jlat),ps(jlon,jlat),2) ! Computes blended saturation below freezing
      q2(jlon,jlat) = min(zq2,zqsat)

      call comp_esk(zeskl,zqsat,t2(jlon,jlat),ps(jlon,jlat),3) ! Computes saturation to water also below freezing
      q2p=min(zq2,zqsat*1.01) ! assures that Q2REL does not exceed 101% (1% more for smoothing graphics)
      eee=ps(jlon,jlat)*q2p/(0.622*(1.-q2p)+q2p)
      q2rel(jlon,jlat) = eee/zeskl

 100  continue

      return
      end
!###############################################################################################################
      subroutine cpbl (nlon, nlat, nlev, richs, zeta, u, v, t, p, q, qcw, qci, tvirt, pbl)

      real, dimension(nlon,nlat)     :: richs, pbl, zorogr
      real, dimension(nlon,nlat,nlev):: t, p, q, tvirt, u, v, qcw, qci, zeta
      real, dimension(nlon,nlev)     :: zqs, zteta, ztetav, zrich

      real, parameter :: tzer=273.15, pzer=1.e5, ezer= 611., rd=287.05, rv=461.51, eps=rd/rv, &
                         cpd=1004.6, cpv=1869.46, rdrcp=rd/cpd, cw = 4186.8,                  &
                         yliv=2834170.5, yliw=333560.5, ylwv=yliv-yliw,                       &
                         ccw1 = (cpv-cw)/rv, ccw2 = ylwv/tzer/rv-ccw1, g = 9.807
      ntop = nlev - 3
      pbl = 1.

! Computation of dry and moist static stability and of effective static stability
! as a function of relative humidity.
! Computation of the Richardson number.

      do jlat = 2, nlat-1

! Comput. of virtual theta

      do jklev = 1, ntop+1
      do jlon = 2, nlon-1
      zconv = exp(rdrcp*(log(pzer)-log(p(jlon,jlat,jklev))))
      zt0t = tzer/t(jlon,jlat,jklev)
      zesk = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t))     ! partial pressure over water
      zqs(jlon,jklev) = zesk*eps/(p(jlon,jlat,jklev)+zesk*(eps-1.))
      zteta (jlon,jklev) = t    (jlon,jlat,jklev)*zconv
      ztetav(jlon,jklev) = tvirt(jlon,jlat,jklev)*zconv
      enddo
      enddo

!  Computation of the Richardson no. depending on dry and moist (saturated) static stability
!  as Durran & Klemp (1982). (Specific humidity is used in place of mixing ratio).

      do jklev = 2, ntop+1   ! loop on half levels
      do jlon = 2, nlon-1
      dthd   = 2.*(ztetav(jlon,jklev)-ztetav(jlon,jklev-1))/(ztetav(jlon,jklev)+ztetav(jlon,jklev-1)) ! dry
      r_up = zqs(jlon,jklev  )
      r_do = zqs(jlon,jklev-1)
      r_av = 0.5*(r_up + r_do)
      t_av = 0.5*(t(jlon,jlat,jklev-1) + t(jlon,jlat,jklev))
      theta_up = zteta(jlon,jklev)
      theta_do = zteta(jlon,jklev-1)
      zaa = (1. + ylwv*r_av/(rd*t_av))/(1. + eps*ylwv**2*r_av/(cpd*rd*t_av**2))
      dthm = zaa*((theta_up-theta_do)*2./(theta_up + theta_do) + ylwv/(cpd*t_av)*(r_up - r_do)) & ! Durran&Klemp
             - q(jlon,jlat,jklev  ) - qcw(jlon,jlat,jklev-1) -  qci(jlon,jlat,jklev  )          &
             + q(jlon,jlat,jklev-1) + qcw(jlon,jlat,jklev-1) +  qci(jlon,jlat,jklev-1)

! Average relative humidity computed giving some more weight to the layer below

      zrhm = 0.55*q(jlon,jlat,jklev-1)/zqs(jlon,jklev-1) + 0.45*q(jlon,jlat,jklev)/zqs(jlon,jklev)
      zcof = max(-24.5 + 25.*zrhm, 0.)    ! zcof=0. for rh<=0.98, zcof=0.5 for rh=1
      zcof = min(zcof, .85)
      dlogthe = zcof*dthm + (1.-zcof)*dthd ! effective stability near saturation
      zrdz  = 1./(zeta(jlon,jlat,jklev)-zeta(jlon,jlat,jklev-1))
      zdtdz = dlogthe*zrdz

      jlonm1 = max(jlon-1,1)
      jlatp1 = min(jlat+1,nlat)
      zdudz = .5*(u(jlon,jlat,jklev)+u(jlonm1,jlat,jklev)-u(jlon,jlat,jklev-1)-u(jlonm1,jlat,jklev-1))*zrdz
      zdvdz = .5*(v(jlon,jlat,jklev)+v(jlon,jlatp1,jklev)-v(jlon,jlat,jklev-1)-v(jlon,jlatp1,jklev-1))*zrdz
      zshear = zdudz**2 + zdvdz**2
      zbuoy = g*dlogthe*zrdz
      zrich(jlon,jklev) = min (zbuoy/(zshear + 1.e-6), 500.)
      enddo
      enddo   ! end of loop on half levels
      zrich(:,1) = richs(:,jlat)  ! Richardson no. at the surface (computed in vtsurf)

! Compute height (over surface) of top of pbl

      do jlon = 2, nlon-1
      jkpbl = 0
       do jklev = 1, nlev/4
       zriav = .5*(zrich(jlon,jklev+1)+zrich(jlon,jklev))
       if (zriav.gt..25) exit
       jkpbl = jklev
       enddo
      if (jkpbl.eq.0.and.zrich(jlon,1).lt..25) jkpbl = 1
      if (jkpbl.gt.0) pbl(jlon,jlat) = zeta(jlon,jlat,jkpbl)+10.
      pbl(jlon,jlat) = min (pbl(jlon,jlat), 2500.)
      enddo

      enddo ! jlat

      return
      end
!=======================================================================
subroutine outgraph(icentre_code,imodel_code,nx,ny,x0,y0,alon0,alat0,dx,dy,                &
           param_discipl,param_categ,param_ind,lev_type,lev1,lev2,idate0,iperiod,a,zf1,zf2)

  implicit none

  integer :: nx, ny, icentre_code, imodel_code, param_discipl, param_categ, param_ind, &
             lev_type, lev1, lev2, idate0(5), iperiod(3)
  real, dimension(nx,ny) :: a
  real :: x0, y0, alon0, alat0, dx, dy, zf1, zf2

! CNR-ISAC-BO: ICENTRE_CODE=253
! BOLAM: IMODEL_CODE=1, MOLOCH: IMODEL_CODE=2, GLOBO: IMODEL_CODE=3

!              PARAM_DISCIPL PARAM_CATEG PARAM_IND
!    T              0            0          0
!    Q              0            1          0

!                          LEV_TYPE
! Ground or water surface       1
! Hybrid level                105
! Isobaric surface (Pa)       100

  open (11,file='data_output.bin',status='unknown',form='unformatted',position='append')

  write (11) icentre_code,imodel_code
  write (11) nx,ny
  write (11) x0,y0,alon0,alat0,dx,dy
  write (11) param_discipl,param_categ,param_ind
  write (11) lev_type,lev1,lev2
  write (11) idate0(1:5),iperiod(1:3)
  write (11) (a(1:nx,1:ny)*zf1+zf2)

  close (11)

return
end subroutine outgraph
!=======================================================================
