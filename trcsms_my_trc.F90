MODULE trcsms_my_trc
   !!======================================================================
   !!                         ***  MODULE trcsms_my_trc  ***
   !! TOP :   Main module of the MY_TRC tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_my_trc
   !!----------------------------------------------------------------------
   !!   'key_my_trc'                                               CFC tracers
   !!----------------------------------------------------------------------
   !! trc_sms_my_trc       : MY_TRC model main routine
   !! trc_sms_my_trc_alloc : allocate arrays specific to MY_TRC sms
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trdmod_oce
   USE trdmod_trc
   
   use sbc_oce

   USE sms_lobster
   USE lbclnk
   USE trdmod_oce
   USE trdmod_trc
   USE iom
   USE prtctl_trc 
   USE fldread
 
  
# include "top_substitute.h90"

   IMPLICIT NONE
   PUBLIC

   PUBLIC   trc_sms_my_trc       ! called by trcsms.F90 module
   PUBLIC   trc_sms_my_trc_alloc ! called by trcini_my_trc.F90 module

   ! parameters
      REAL(wp) ::  q_myt1 !production rate of tracer1
      REAL(wp) ::  q_myt2 !production rate of tracer2
      REAL(wp) ::  coeff 
  
      REAL(wp), POINTER, DIMENSION(:,:,:,:)  ::   s_pa      !sinking rate
      REAL(wp),POINTER, DIMENSION(:,:,:,:)   ::  s_th      !sinking rate
      REAL(wp), POINTER, DIMENSION(:,:,:,:)  ::   k_th      !adsorption rate
      REAL(wp),POINTER, DIMENSION(:,:,:,:)   ::  k_pa      !adsorption rate
      REAL(wp) ::  j_pa      !desorption rate
      REAL(wp), POINTER, DIMENSION(:,:,:)  ::   rfse3t

   ! Defined HERE the arrays specific to MY_TRC sms and ALLOCATE them in trc_sms_my_trc_alloc
      CHARACTER(len=100), PUBLIC ::   cn_dir       = './'    !: Root directory for location of river file
      TYPE(FLD_N)                ::   sn_Pa_diss ! information about the diss. Pa file to be read
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: Pa_diss   ! Array that receives values from netCDF
      TYPE(FLD), ALLOCATABLE, DIMENSION(:)    :: sf_Pa_diss ! structure variable (PUBLIC for TAM)

      TYPE(FLD_N)                             :: sn_Pa_part
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: Pa_part
      TYPE(FLD), ALLOCATABLE, DIMENSION(:)    :: sf_Pa_part


      TYPE(FLD_N)                             :: sn_Th_diss
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: Th_diss
      TYPE(FLD), ALLOCATABLE, DIMENSION(:)    :: sf_Th_diss 


      TYPE(FLD_N)                             :: sn_Th_part
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: Th_part 
      TYPE(FLD), ALLOCATABLE, DIMENSION(:)    :: sf_Th_part 

      TYPE(FLD_N)                             :: sn_j_th 
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: j_th 
      TYPE(FLD), ALLOCATABLE, DIMENSION(:)    :: sf_j_th 

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcsms_my_trc.F90 3294 2012-01-28 16:44:18Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_my_trc( kt )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_my_trc  ***
      !!
      !! ** Purpose :   main routine of MY_TRC model
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean e-step index
      INTEGER ::   jn,jk,jj,ji, ibot,ibotm1, ibotm2  ! dummy loop index 
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztrmyt


!!----------------------------------------------------------------------

      IF( nn_timing == 1 )  CALL timing_start('trc_sms_my_trc')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_my_trc:  MY_TRC model really??'
      IF( l_trdtrc )  CALL wrk_alloc( jpi, jpj, jpk, ztrmyt)
  

       coeff=1.0/(365.0*86400.0)
       j_pa=0.31*coeff

       DO jk=1,jpk
       rfse3t(:,:,jk)=1.0/fse3t(:,:,jk)
       END DO



!add variable ice_birth

       !tracer ice for Pa
       DO jk=1,16
         trn(:,:,jk,jpmyt5)=fr_i
       END DO

       !tracer ice for Th
       DO jk=1,16
         trn(:,:,jk,jpmyt6)=fr_i
       END DO

       !set sinking/adsorption rates for Pa and Th
       DO jk=1,jpk
         s_pa(:,:,jk,jpmyt1)=tmask(:,:,jk)*(601*trb(:,:,jk,jpmyt5)+500)*coeff
         s_th(:,:,jk,jpmyt2)=tmask(:,:,jk)*(226*trb(:,:,jk,jpmyt6)+500)*coeff
       END DO
       DO jk=1,jpk
         k_pa(:,:,jk,jpmyt1)=tmask(:,:,jk)*(-0.05*trb(:,:,jk,jpmyt5)+0.06)*coeff
         k_th(:,:,jk,jpmyt2)=tmask(:,:,jk)*(-0.67*trb(:,:,jk,jpmyt6)+0.75)*coeff
       END DO

       !add sinking term for the 'ice'
       DO jk=2,jpk
         tra(:,:,jk,jpmyt5)=tmask(:,:,jk)*(tra(:,:,jk,jpmyt5)+&
                            (s_pa(:,:,jk-1,jpmyt1)*trb(:,:,jk-1,jpmyt5)&
                            -s_pa(:,:,jk  ,jpmyt1)*trb(:,:,jk  ,jpmyt5))*rfse3t(:,:,jk-1))
       END DO
       DO jk=2,jpk
         tra(:,:,jk,jpmyt6)=tmask(:,:,jk)*(tra(:,:,jk,jpmyt6)+&
                            (s_th(:,:,jk-1,jpmyt2)*trb(:,:,jk-1,jpmyt6)&
                            -s_th(:,:,jk  ,jpmyt2)*trb(:,:,jk,jpmyt6))*rfse3t(:,:,jk-1))
       END DO

!add bc for diss.
       DO jk=1,jpk
         WHERE ((gphit < 69))
                trn(:,:,jk,jpmyt1)= Pa_diss(:,:,jk)
                trn(:,:,jk,jpmyt2)= Th_diss(:,:,jk)
         END WHERE
       END DO

!add bc for part.
       DO jk=1,jpk
         WHERE ((gphit < 69))
                trn(:,:,jk,jpmyt3)= Pa_part(:,:,jk)
                trn(:,:,jk,jpmyt4)= Th_part(:,:,jk)
         END WHERE
       END DO




 !calculate diss Pa equations
       DO jk=1,jpk
         WHERE (gphit >= 69)
         tra(:,:,jk,jpmyt1) = tmask(:,:,jk)*( q_myt1  &
                             - k_pa(:,:,jk,jpmyt1)*trb(:,:,jk,jpmyt1) &
                             + j_pa*trb(:,:,jk,jpmyt3))
         END WHERE
       END DO

 !calculate part Pa equations
       DO jk=1,jpk
         WHERE (gphit >= 69)
         tra(:,:,jk,jpmyt3) = tmask(:,:,jk)*(k_pa(:,:,jk,jpmyt1)*trb(:,:,jk,jpmyt1) &
                            - j_pa*trb(:,:,jk,jpmyt3))     
         END WHERE
       END DO

       DO jk=2,jpk
            WHERE (gphit >= 69)
             tra(:,:,jk,jpmyt3) = tmask(:,:,jk)*(tra(:,:,jk,jpmyt3)+ &
                              (s_pa(:,:,jk-1,jpmyt1)*trb(:,:,jk-1,jpmyt3)&
                              -s_pa(:,:,jk  ,jpmyt1)*trb(:,:,jk  ,jpmyt3))*rfse3t(:,:,jk-1))
            END WHERE
       END DO

!calculate diss Th equations
       DO jk=1,jpk
         WHERE (gphit >= 69)
         tra(:,:,jk,jpmyt2) =  tmask(:,:,jk)*(q_myt2  &
                             - k_th(:,:,jk,jpmyt2)*trb(:,:,jk,jpmyt2) &
                             + j_th(:,:)*coeff*trb(:,:,jk,jpmyt4))
         END WHERE
       END DO

!calculate part Th equations
       DO jk=1,jpk
         WHERE (gphit >= 69)
         tra(:,:,jk,jpmyt4) = tmask(:,:,jk)*(k_th(:,:,jk,jpmyt2)*trb(:,:,jk,jpmyt2) &
                            - j_th(:,:)*coeff*trb(:,:,jk,jpmyt4))
         END WHERE
       END DO

       DO jk=2,jpk
            WHERE (gphit >= 69) 
             tra(:,:,jk,jpmyt4) = tmask(:,:,jk)*(tra(:,:,jk,jpmyt4)+ &
                                (s_th(:,:,jk-1,jpmyt2)*trb(:,:,jk-1,jpmyt4)&
                                -s_th(:,:,jk  ,jpmyt2)*trb(:,:,jk  ,jpmyt4))*rfse3t(:,:,jk-1))
            END WHERE
       END DO




      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_sms_my_trc')
      !
   END SUBROUTINE trc_sms_my_trc


   INTEGER FUNCTION trc_sms_my_trc_alloc()
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_sms_my_trc_alloc  ***
      !!----------------------------------------------------------------------
      !
      INTEGER           ::   ierror ! temporary integer
      ! ALLOCATE here the arrays specific to MY_TRC
      ! ALLOCATE( tab(...) , STAT=trc_sms_my_trc_alloc )
      ALLOCATE( rfse3t(jpi,jpj,jpk),k_th(jpi,jpj,jpk,jptra), s_pa(jpi,jpj,jpk,jptra),k_pa(jpi,jpj,jpk,jptra), s_th(jpi,jpj,jpk,jptra))


      ALLOCATE( Pa_diss(jpi,jpj,jpk), STAT=trc_sms_my_trc_alloc )
      ALLOCATE( sf_Pa_diss(1), STAT=ierror )         ! Create sf_Pa_diss structure (Pa_diss inflow)

      ALLOCATE( Pa_part(jpi,jpj,jpk), STAT=trc_sms_my_trc_alloc )
      ALLOCATE( sf_Pa_part(1), STAT=ierror )

      ALLOCATE( Th_diss(jpi,jpj,jpk), STAT=trc_sms_my_trc_alloc )
      ALLOCATE( sf_Th_diss(1), STAT=ierror )

      ALLOCATE( Th_part(jpi,jpj,jpk), STAT=trc_sms_my_trc_alloc )
      ALLOCATE( sf_Th_part(1), STAT=ierror )

      ALLOCATE( j_th(jpi,jpj), STAT=trc_sms_my_trc_alloc )
      ALLOCATE( sf_j_th(1), STAT=ierror )

!      IF(lwp) WRITE(numout,*)
!      IF(lwp) WRITE(numout,*) '          Pa_diss inflow read in a file'
      IF( ierror > 0 ) THEN
            CALL ctl_stop( 'trc_sms_my_trc_alloc: unable to allocate sf_Pa_diss structure' )   ;   RETURN
      ENDIF
      ALLOCATE( sf_Pa_diss(1)%fnow(jpi, jpj, jpk)   )
      ALLOCATE( sf_Th_diss(1)%fnow(jpi, jpj, jpk)   )
      ALLOCATE( sf_Th_part(1)%fnow(jpi, jpj, jpk)   )
      ALLOCATE( sf_Pa_part(1)%fnow(jpi, jpj, jpk)   )

      ALLOCATE( sf_j_th(1)%fnow(jpi, jpj, 1)   )

      trc_sms_my_trc_alloc = 0      ! set to zero if no array to be allocated
      !
      IF( trc_sms_my_trc_alloc /= 0 ) CALL ctl_warn('trc_sms_my_trc_alloc : failed to allocate arrays')
      !
   END FUNCTION trc_sms_my_trc_alloc


#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No MY_TRC model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sms_my_trc( kt )             ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms_my_trc: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_my_trc
#endif

   !!======================================================================
END MODULE trcsms_my_trc
