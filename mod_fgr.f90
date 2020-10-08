Module mod_fgr
!! JPCA_117_6196
implicit none
real*8, parameter :: clight=2.99792458D10,av=6.0221367D23,hbar=1.05457266D-34
real*8, parameter :: mass_h=1.007825d0,kb=1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=4184.d0
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter :: au2kg=9.10938291d-31,au2J=4.35974417d-18,au2m=5.2917721092d-11
complex*16,parameter :: iota = (0.d0,1.d0)
real*8 pi,wave_to_J

!! Common Variables
real*8 V_coup,V_exothermicity
real*8 omg1,V_reorg1,g_coup1
real*8 omg2,V_reorg2,g_coup2
real*8 gamma
real*8 temperature,mass,x_cr,beta
real*8 VER_rate

!! integration
real*8 wmax,tmax
integer nclass,n_t

!! diag
integer nold
real*8,allocatable:: work(:)
integer,allocatable:: iwork(:),isuppz(:)

contains
!---------------------------------------------------------- 
!---------------------------------------------------------- 

subroutine setup
  implicit none
  character st_ch
  real*8 c_0,c_e

  pi=dacos(-1.d0)
  wave_to_J=2*pi*clight*hbar
  open(10,file="fgr.inp")
  read(10,*) mass
  read(10,*) V_coup
  read(10,*) V_exothermicity
  read(10,*) omg1
  read(10,*) V_reorg1
  read(10,*) omg2
  read(10,*) gamma
  read(10,*) VER_rate
  read(10,*) nclass
  read(10,*) wmax
  read(10,*) n_t
  read(10,*) tmax
  read(10,*) st_ch
  close(10)
  !----------------------------------------------------------

  if(st_ch.ne.'x') then
    write(6,*) "problem in reading input file"
    stop
  endif

  !---------------------------------------------------------- 

  mass=mass*au2kg
  omg1=omg1*(2*pi*clight)
  omg2=omg2*(2*pi*clight)
  gamma=gamma*(2*pi*clight)
  wmax=wmax*(2*pi*clight)
  V_exothermicity=V_exothermicity*wave_to_J
  V_coup=V_coup*wave_to_J
  V_reorg1=V_reorg1*wave_to_J
  beta=1.d0/(kb*temperature)

  g_coup1=dsqrt(V_reorg1*mass*omg1**2/2.d0)
  call calculate_lambda

end subroutine setup
!---------------------------------------------------------- 

subroutine main
  implicit none
  integer i
  real*8 log_gama
  real*8 k_FGR,k_TST,k_M

  open(20,file="FGR_rate.out")
  write(6,*) "VER_rate (ps-1)","     k_FGR (ps-1)","     k_M (ps-1)"
  write(20,*) "VER_rate (ps-1)","    k_FGR (ps-1)","     k_M (ps-1)"
  do i=1,20
    VER_rate=1.d12+39.d12*(i-1)/19.d0
    call calculate_lambda
    temperature=400.d0!+250.d0*(i-1)/5.d0
    beta=1.d0/(kb*temperature)
    k_TST=omg1/(2*pi)*dexp(-(V_exothermicity-V_reorg1)**2/(4*V_reorg1*kb*temperature))
    k_M=2*pi*V_coup**2/(hbar*dsqrt(4*pi*V_reorg1*kb*temperature))*dexp(-(V_exothermicity-V_reorg1)**2/(4*V_reorg1*kb*temperature))

    call compute_k_FGR(k_FGR)
    write(6,*) VER_rate/1.d12,k_FGR/1.d12,k_M/1.d12
    write(20,*) VER_rate/1.d12,k_FGR/1.d12,k_M/1.d12
  enddo
  close(100)

end subroutine main
!---------------------------------------------------------- 

subroutine compute_k_FGR(k_FGR)
  implicit none
  real*8,intent(out) :: k_FGR
  integer i,j
  real*8 t,tmin,dt
  real*8 omg(nclass),ck(nclass),wt(nclass)
  complex*16 su1,su2,fac(nclass)
  real*8 t1,t2
  real*8 kk

  call cpu_time(t1)

  call convert_spectra(omg,ck)

  tmin=0.d0
  dt=(tmax-tmin)/dfloat(n_t-1)
  
  su2=0.d0
  !call omp_set_num_threads(8)
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,t,wt,fac,su1)
  !$OMP DO SCHEDULE(STATIC) REDUCTION(+:su2)
  do i=1,n_t
    t=tmin+(i-1)*dt
    wt=omg*t
    fac=(1-dcos(wt))/dtanh(0.5d0*beta*hbar*omg)-iota*dsin(wt)
    su1=2/mass*sum(ck**2/omg**3*fac)/hbar
    su2=su2+cdexp(-su1-iota*V_exothermicity*t/hbar)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  su2=su2*dt

  k_FGR=2*V_coup**2*dble(su2)/hbar**2

  call cpu_time(t2)
  !write(6,*) t2-t1

end subroutine compute_k_FGR
!-----------------------------------------------------------------  

function density(w)
  implicit none
  real*8 density,w

  density=0.5*V_reorg2*omg2*omg2 * gamma*w/((w*w-omg2*omg2)**2+(gamma*w)**2)
!  density=mass*omg1*w

end function density
!-----------------------------------------------------------------  

subroutine convert_spectra(omg_n,ck)
  implicit none
  real*8,intent(out):: omg_n(nclass),ck(nclass)
  integer i,j
  real*8 mat(nclass,nclass),en(nclass),vect(nclass,nclass)
  real*8 ck_brown(nclass),omg(nclass)
  real*8 delw

  delw=wmax/real(nclass)

  do i=1,nclass
    omg(i)=i*wmax/real(nclass)
    ck(i)=dsqrt(2*density(omg(i))*mass*omg(i)*delw/pi)
  enddo
  ck(1)=0.d0

  mat=0.d0
  mat(1,1)=0.5*mass*omg1**2+sum(ck**2/(2*mass*omg**2))

  do i=2,nclass
    mat(i,i)=0.5*mass*omg(i)**2
  enddo
  do i=2,nclass
    mat(1,i)=ck(i)/2.d0
    mat(i,1)=ck(i)/2.d0
  enddo

  call diag(mat,nclass,en,vect,nclass)

  omg_n=dsqrt(2*en/mass)

  ck=vect(1,:)*g_coup1

end subroutine convert_spectra
!-----------------------------------------------------------------  

subroutine calculate_lambda
  implicit none
  real*8 w

  w=omg1

  V_reorg2=VER_rate*mass*omg1
  V_reorg2=V_reorg2 * 2 * ((w*w-omg2*omg2)**2+(gamma*w)**2)/(omg2**2*gamma*w)

  g_coup2=dsqrt(V_reorg2*mass*omg2**2/2.d0)

end subroutine calculate_lambda
!-----------------------------------------------------------------  

subroutine diag(mat,n,eigen_value,eigen_vect,m_values)
  !! Diaganalizing matrix using dsyevr. First m_values eigen values and eigenvectors computed.
  !! The module's common variables should contain:

  !! Initialize nold=0 

  !! nold makes sure that everytime value of n changes, work and iwork are re-allocated for optimal performance.
  !! mat is destroyed after use.

  implicit none
  integer,intent(in) :: n,m_values
  real*8,intent(out) :: eigen_value(n),eigen_vect(n,m_values)
  real*8,intent(inout) :: mat(n,n)
  real*8 vl,vu,abstol
  integer il,iu,info,m,AllocateStatus
  integer lwork,liwork

  vl=0.d0;vu=0.d0   !! not referenced
  il=1;iu=m_values
  abstol=0.d0
  info=0

  if(nold.ne.n .or. .not.allocated(work) .or. .not.allocated(iwork) .or. .not.allocated(isuppz)) then
  !if(nold.ne.n) then
    lwork=-1;liwork=-1
    if(allocated(isuppz))deallocate(isuppz)
    if(allocated(work))deallocate(work)
    if(allocated(iwork))deallocate(iwork)
    allocate(isuppz(2*m_values),work(n),iwork(n))
    call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
    lwork=nint(work(1)); liwork=iwork(1)
    deallocate(work,iwork)
    allocate(work(lwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    allocate(iwork(liwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    nold=n
  endif

  lwork=size(work)
  liwork=size(iwork)

  call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
  if(info.ne.0) then
    write(6,*) "problem in diagonalization",info
    stop
  endif

end subroutine diag
!---------------------------------------------------------- 

End Module mod_fgr
