c
c
      subroutine dskout(ll)
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     At the end of the run, this routine at the option of the user 
c     writes out to disk file 'idskf' the namelist input deck 
c     and various computed quantities, including the distn functions.
c     Also, write out to disk file 'idskrf' rf diffusion coefficients,
c     and related quantities.
c..................................................................

      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE

      character*80 line
      
CMPIINSERT_IF_RANK_NE_0_RETURN

      if (lrzmax.le.1) then
        if((idskf.eq. "disabled".and.idskrf.eq."disabled")
     1    . or. n.ne.nstop+1)  return
      else
        if((idskf.eq. "disabled".and.idskrf.eq."disabled")
     1    . or. n.ne.nstop)  return
      endif
      if(ll.gt.1)  go to 4
      if(idskf.ne."disabled")
     1  open(unit=4,file=idskf,delim='apostrophe',status='unknown')
      if(idskrf.ne."disabled")
     1  open(unit=5,file=idskrf,delim='apostrophe',status='unknown')
ccc      close(unit=2) ! YuP: Why here?
      ilen=0

c..................................................................
c     The input namelist file is transcribed onto the beginning
c     of file idskf and/or idskrf
c..................................................................

      open(unit=2,file='cqlinput',delim='apostrophe',status='old')
 1    read(2,1003) line
      if (line(1:3).eq."end") go to 3
      if(idskf.ne."disabled") write(4,1003) line
      if(idskrf.ne."disabled") write(5,1003) line
      go to 1
 3    close(unit=2)
      if(idskf.ne."disabled") 
     +     write(4,1006) '***Begin computed output***'
      if(idskrf.ne."disabled") 
     +     write(5,1006) '***Begin computed output***'
 4    continue 

c..................................................................
c     In the following disk write to file named idskf:
c     This subroutine is called to write data for each FP'd
c          flux surface.
c     ll=  FP flux surface number 
c          (ll=1:lrors, lrors.le.lrzmax, see cqlinput_help))
c          (lrindx(ll) gives flux surface number on the full
c                     radial mesh.(lrindx(ll)=ll if lrors=lrzmax,
c                     and using cql3d mode(cqlpmod="disabled"))).
c     iy,jx= dimensions in theta and u(momentum/mass)
c           (In the case where iy varies with ll, iy(1) will be greatest.)
c     lrors= number of flux surfaces FP'd.
c     lrzmax= number of flux surfaces, including any not FP'd.
c     x = momentum-per-mass(nomalized to maximum 1.)
c         at  each flux surface lrindx(ll)
c     y = theta(radians) mesh at  each flux surface lrindx(ll)
c     rovera= normalized radius (ll)
c             (rho, see Hinton and Haseltine for non-circ).
c             (generally ~sqrt(tor. flux), other coords available.)
c     elecfld = toroidal electric field (volts/cm)
c     bthr(lr_) - the poloidal magnetic field at theta-poloidal = pi/2.
c     btoru(lr_) - the toroidal magnetic field at the same position. 
c     bthr0(lr_), btor0(lr_), are the poloidal and  toroidal
c         magnetic fields at the outer midplane of the flux surface. 
c     reden= electron density at minimum B point on flux surface.
c     temp= initial electron temperature (keV)
c     radmin= plasma minor radius (cms).
c     vnorm= normalization momentum-per-mass (maximum on grid) (cm/sec)
c     vmaxdvt= vnorm/(temp/mass)**0.5
c     eovedd= electric  field, normalized to Driecer field 
c             (calc'd in sub restvty).
c
c     distribution function normalised so that
c         integral( (dx)**3 f) = density at minumum B point.
c
c..................................................................

      if(idskf.ne."disabled") then
        write(4,1004)  ll, iy,jx,lrors,lrzmax,ngen
        write(4,1004)  itl,itu
        write(4,1005)  (x(j),j=1,jx)
        write(4,1005)  (y(i,ll),i=1,iy)
        do 1000 k=1,ngen
          write(4,1005)  bnumb(k),fmass(k)
          write(4,1005)  rovera(lrindx(ll)),elecfld(lrindx(ll)),
     +                   bthr(lrindx(ll)),btoru(lrindx(ll))
          write(4,1005)  bthr0(lrindx(ll)),btor0(lrindx(ll)),
     +                   reden(k,lrindx(ll)),temp(k,lrindx(ll))
          vmaxdvt=vnorm/(4.19e7*sqrt(2.*1000.*temp(k,lrindx(ll))))
          write(4,1005)  radmin,vnorm,vmaxdvt,eovedd
          write(4,1005)  ((f(i,j,k,ll),i=1,iy),j=1,jx)
 1000   continue
      endif


c..................................................................
c     In the following disk write:
c     vptb is lambda=u_parallel_o * tau_bounce.
c     (tau_bounce= portion of bounce time from the outer equatorial
c     plane to the bounce point or the inner equatorial plane).
c     temp1 is the bounce averaged uu-diffusion coefficient (cgs units).
c..................................................................
      if(idskrf.ne."disabled") then
        write(5,1004)  ll,iy,jx,lrors,lrzmax
        write(5,1005)  (y(i,ll),i=1,iy)
        write(5,1005)  (x(j),j=1,jx)
        write(5,1005)  rovera(lrindx(ll)),vnorm
        write(5,1005)  (vptb(i,lrindx(ll)),i=1,iy)
        vn2=vnorm*vnorm
        do 10 k=1,mrfn
          vmaxdvt=vnorm/(4.19e7*sqrt(2.*1000.*temp(k,lrindx(ll))))
          write(5,1005) reden(k,lrindx(ll)),temp(k,lrindx(ll)),vmaxdvt
          do 11 i=1,iy
 11       temp1(i,1)=0.0
          do 12 j=2,jx
            x2=x(j)**2
            do 13 i=1,iy
 13         temp1(i,j)=urfb(i,j,indxlr(lrindx(ll)),k)*vn2/
     /          (vptb(i,lrindx(ll))*x2)
 12       continue
          write(5,1005) ((temp1(i,j),i=1,iy),j=1,jx)
 10     continue
      endif
 1003 format(a80)
 1004 format(16i5)
 1005 format(5e16.8)
 1006 format(a27)
      if(idskf.ne."disabled".and.ll.eq.lrors)  close(unit=4)
      if(idskrf.ne."disabled".and.ll.eq.lrors)  close(unit=5)
      return
      end subroutine dskout

!=======================================================================
!=======================================================================


      subroutine read_data_files(filenm,jtm,ipresent,t_data,kopt)
      !YuP[2021-01-21] read data files.
      !(Initial purpose - coupling with NIMROD. 
      ! Can be extended to coupling with other codes.)
      ! For coupling with NIMROD, 
      ! each file contains data at one time slice.
      !implicit none
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h' !contains enein_t(njenea,ntotala,nbctimea),
               !tein_t(njenea,nbctimea), etc.
CMPIINSERT_INCLUDE

      !INPUT: filenm, jtm, kopt;  Also nstates (from comm.h)
      !OUTPUT: ipresent (=1 if the filenm is present, =0 otherwise)
      !t_data= t[sec], obtained from reading 1st line in given data file.
      !Other output: the data that is read from file is copied into
      ! CQL3D arrays like enein_t(), tein_t(), etc.
  
      character*128,INTENT(IN) :: filenm
      integer,INTENT(IN) :: jtm !Index for data file
      integer,INTENT(IN) :: kopt
      integer,INTENT(OUT):: ipresent
      real*8,INTENT(OUT) :: t_data
      integer kode,istat,iunit ! local
      integer idum !local
      real*8 dum1,dum2,dum3 ! local
      character*8 cdum !local
      integer njene_data !local
      integer irow !local
      character*64 format_data !local, To form a format line
      ! For reading columns :
      integer ix_data
      real*8 r_data,psi_data,b_data,e_data,curr_data,elecd_data
      real*8 dene_data,denD0_data,deni_data,denz_data
      real*8 ti_data,te_data
      real*8,dimension(:),allocatable :: denzz_data ![0:nstates]
      save denzz_data
      
      ipresent=0  !To be changed to 1 if filenm is found
      t_data=0.d0 !To be found from reading data: time slice [sec]
      iunit=14
      open(unit=iunit,file=filenm,status='old',iostat=kode)
      if (kode.ne.0) then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)'subr.read_data_files:',filenm,' cannot be opened'
CMPIINSERT_ENDIF_RANK         
         ipresent=0
      else
         ipresent=1
      endif

 10   format(i6) ! First number in first line: Radial grid size.
 11   format(i6,4x,a6,i6,3x,a6,d12.5) !The whole first line.
      
      if(kopt.eq.0  .and. ipresent.eq.1)then
        !Only check the presence of file; 
        !Also check the radial grid.  Do not read data yet.
        read(iunit,10) njene_data
        njene=njene_data !To be saved into name_decl.h
        !This value is checked in ainsetva against njenea
        goto 200 !-> return/exit
      endif
      
      if(kopt.eq.1 .and. ipresent.eq.1)then
        idum=0
        cdum='none'
        dum1=0.d0

        read(iunit,11,iostat=kode) idum,cdum,idum,cdum,dum1
        if (kode.ne.0) then
CMPIINSERT_IF_RANK_EQ_0
          WRITE(*,*)'subr.read_data_files:',TRIM(filenm)
          WRITE(*,*)'subr.read_data_files: Corrupted 1st line'
          WRITE(*,*)'time[sec]==dum1=',dum1
CMPIINSERT_ENDIF_RANK         
          STOP 'In read_data_files: cannot read 1st line'
        endif
        t_data=dum1 !OUTPUT: to be used to fill bctime() array
        
        !allocate storage for array to read data
        if (.NOT. ALLOCATED(denzz_data)) then
          allocate(denzz_data(0:nstates),STAT=istat)
          !To read data for impurity, all ionization states,
          !starting from 0 (neutral)
        endif
 
        !Allocate time-dependent variable 
        !It is set as pointer in comm.h, 
        !and stored in  common/impur_read_data/dens_imp_t
        !Only needed when read_data='nimrod'.
        !Note that dens_imp_t() is not a namelist var, it is just a storage
        !for density of ionization states at each time slice.
        if(.NOT. ASSOCIATED(dens_imp_t)) then
          allocate(dens_imp_t(0:nstates,1:njene,1:nbctime),STAT=istat)
        endif
    
        !Form format specs.  ONLY FOR NIMROD data files: 
        write(cdum,'(i4)') (nstates+1) ! columns corr to 0:nstates
        format_data='(1x,i4,1x,12d16.9,'//TRIM(cdum)//'(d16.9))' 
        !write(*,*) 'format_data=',format_data
        !write(*,*) 't_data[sec]=',t_data
        
        ! Read row by row
        read(iunit,*) !Skip one row (with headings)
        do irow=1,njene       
           read(iunit,format_data,iostat=kode) 
     &         ix_data,
     &         r_data,psi_data,b_data,e_data,curr_data,elecd_data,
     &         dene_data,denD0_data,deni_data,denz_data,
     &         ti_data,te_data,
     &         (denzz_data(idum),idum=0,nstates) 
           if (kode.ne.0) then
CMPIINSERT_IF_RANK_EQ_0
             WRITE(*,*)'subr.read_data_files: cannot read line=',irow+2
CMPIINSERT_ENDIF_RANK         
             STOP 'In read_data_files: cannot read line'
           endif
     
           !1.Verified:  denz_data ~ sum(denzz_data), up to 7-8 digits.
           ![ denz_data is from the column designated as "nz" ] 
           !We don't need denz_data by itself, so - no adjustment here.
           !Also noticed that sometimes denzz_data(idum) have 
           !small negative values (at initial time slices).
           !---> Resetting to be not lower than 0.0 :
           do kstate=0,nstates
              denzz_data(kstate)= max(denzz_data(kstate),0.d0)
           enddo
           !Also, just in case:
           deni_data= max(deni_data,0.d0) ! For main ion "ni"
           
           !2.TEST for quasineutrality.
           sum_ni_Zi=deni_data !column "ni", which is for the main ion (D+)
           sum_nimp_zstate=0.d0 !To count electrons from impurity, all Zstates
           do kstate=1,nstates !These are additional ions from impur.source.
             sum_nimp_zstate= sum_nimp_zstate
     &                       +denzz_data(kstate)*bnumb_imp(kstate)
           enddo ! kstate
           sum_ni_Zi= sum_ni_Zi +sum_nimp_zstate !total number of free electrons
           err1= abs(dene_data-sum_ni_Zi)/dene_data !Supposed to be ~0
!           if(err1.gt.1.d-7)then
!             write(*,'(a,2e16.9)')
!     &       '  dene_data,sum_ni_Zi=',dene_data,sum_ni_Zi
!           endif
           !RESULTS: At t_data>0, dene_data = sum_ni_Zi up to 7 digit,
           ! but at t_data=0 [the very 1st data file]
           ! the value of dene_data [column "ne"] is corrupted - 
           ! it is set to 0.4e20 at all radial points 
           ! (supposed to be the edge value only)
           !---> Resetting ne to have quasineutrality:
           dene_data= sum_ni_Zi 
           !END TEST
           
           !3. LOWER LIMIT for Te and Ti
           !Noticed that Te and Ti may go down to ~0.3 eV in data tables.
           !---> Resetting to be not lower than temp_min_data :
           te_data= max(te_data,temp_min_data*1d3) 
           ti_data= max(ti_data,temp_min_data*1d3) 
           ! ti_data is in eV, while temp_min_data is in keV in cqlinput namelist
           
           !--------- Now populate CQL3D arrays
           if(irow.eq.1) r0_data=r_data !Assuming 1st data point is at magn.axis
           ryain(irow)= (r_data-r0_data) !Here in [m]
           !The above definition will be normalized by (r_data_max-r_data_min)
           ! after all rows are read.
           
           tein_t(irow,jtm)= te_data*1.d-3 !converted to keV
           tiin_t(irow,jtm)= ti_data*1.d-3 !converted to keV
           !Main Maxwellian species (only one, for now, D+):
           enein_t(irow,kionm(1),jtm)= deni_data*1.d-6 !converted to cm^-3
           !Electron species:
           if(kelecg.ne.0)then !General e species
           enein_t(irow,kelecg,jtm)= dene_data*1.d-6 !converted to cm^-3
           endif
           if(kelecm.ne.0)then !Maxwellian e
           enein_t(irow,kelecm,jtm)= dene_data*1.d-6 !converted to cm^-3
           endif
           
           !Electric field [V/cm]
           elecin_t(irow,jtm)= e_data*1.d-2 !converted V/m to V/cm
           !Note: in NIMROD file, e_data is labeled as "E.B/B",
           !so - is it parallel? Should we reset efflag="parallel" ?
           !Also the proper sign is of some concern.
           !Maybe need to use bsign here? 
           !write(*,'(a,e12.6)')'elecin_t()=',elecin_t(irow,jtm)
           
           !Impurity - densities of each ionization state
           do kstate=0,nstates
              dens_imp_t(kstate,irow,jtm)=denzz_data(kstate)*1.d-6 !to cm^-3
           enddo
           
        enddo ! irow
        
        ryain(:)= ryain(:)/(r_data-r0_data) !For radcoord='rminmax' only!
        !assuming last data point is at rho=1 (psi=psilim)
        !write(*,*)'ryain=',ryain(1:njene) ! checked: [0;1] range
        !write(*,*)kelecg,kelecm
        
      endif ! kopt.eq.1 .and. ipresent.eq.1
      
      
 200  continue
      close(iunit)

      return
      end subroutine read_data_files
