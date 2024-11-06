
subroutine user_create_perturbation

  use ModGITM
  use ModInputs

  NDensityS(nLons/2,nLats/2,2,1:3,1) = NDensityS(nLons/2,nLats/2,2,1:3,1)*50.0
  temperature(nLons/2,nLats/2,2,1) = temperature(nLons/2,nLats/2,2,1)*100.0

end subroutine user_create_perturbation


subroutine user_perturbation

  use ModGITM
  use ModInputs
  use ModNumConst
  use ModKind, ONLY: Real8_
  use ModTime
  use ModSources, only: UserHeatingRate

  implicit none

  real         :: latcenter, loncenter, amp
  real         :: latwidth, lonwidth
  real         :: f(-1:nLons+2,-1:nLats+2)
  integer      :: iBlock, iSpecies, iLat, iLon, iAlt
  real(Real8_) :: PerturbTimeStart, PerturbTimeEnd, MidTime
  real         :: tsave = 0.0
  real         :: dla, dlo, lac, loc
  
  UserHeatingRate = 0.0

  ! Start 1 hour1 after start of the simulation
  PerturbTimeStart = StartTime + 1.0*3600.0+60.0
  ! End 1 minute after that
  PerturbTimeEnd   = PerturbTimeStart + 60.0

  dla = latitude(2,1) - latitude(1,1)
  lac = (LatEnd+LatStart)/2.0
  dlo = longitude(2,1) - longitude(1,1)
  loc = (LonEnd+LonStart)/2.0

  latcenter = lac + dla/2.0
  loncenter = loc + dlo/2.0
  latwidth  = dla/2.0
  lonwidth  = dlo/2.0

  if  (CurrentTime < PerturbTimeStart) then
     tsave = sum(temperature(1:nLons,1:nLats,1,1:nBlocks)) / &
          (nLons*nLats*nBlocks)
  endif

  DuringPerturb = .false.
  
  if  (CurrentTime >= PerturbTimeStart .and. &
       CurrentTime <  PerturbTimeEnd) then

     MidTime = (PerturbTimeStart + PerturbTimeEnd)/2.0

     DuringPerturb = .true.

     amp = exp(-((CurrentTime-MidTime)/(PerturbTimeEnd-PerturbTimeStart)*5)**2)
!     write(*,*) "Perturbing!", tsave, amp,latcenter/3.1415*180.0,&
!          loncenter/3.1415*180.0, latwidth, lonwidth

     do iBlock = 1, nBlocks
        do iLon = 1, nLons
           do iLat = 1, nLats
              if ((abs(latitude(iLat,iBlock)-latcenter) < 4*latwidth) .and. &
                  ( abs(longitude(iLon,iBlock)-loncenter) < 4*lonwidth)) then
                 f(iLon,iLat) = & !amp*&
                      exp(-((latitude(iLat,iBlock)-latcenter)/latwidth)**2) * &
                      exp(-((longitude(iLon,iBlock)-loncenter)/lonwidth)**2)
              else
                 f(iLon,iLat) = 0.0
              endif

              iAlt = 1
              UserHeatingRate(iLon,iLat,iAlt,iBlock) =  &
                   1000.0 * 4.184e6 * f(iLon,iLat) / &
                   cellvolume(iLon,iLat,iAlt,iBlock) / &
                   TempUnit(iLon,iLat,iAlt) / &
                   cp(iLon,iLat,iAlt,iBlock) / &
                   rho(iLon,iLat,iAlt,iBlock)
           enddo
        enddo

!        temperature(:,:,-1,iBlock) = tsave + 5.0 * f * tsave
!        temperature(:,:, 0,iBlock) = tsave + 5.0 * f * tsave
!        temperature(:,:, 1,iBlock) = tsave + 5.0 * f * tsave
!        do iSpecies = 1, nSpecies
!           VerticalVelocity(:, :, -1, iSpecies, iBlock) = 2*f
!           VerticalVelocity(:, :,  0, iSpecies, iBlock) = 2*f
!           VerticalVelocity(:, :,  1, iSpecies, iBlock) = 2*f
!        enddo
!        Velocity(:,:,-1, iUp_, iBlock) = 2*f
!        Velocity(:,:, 0, iUp_, iBlock) = 2*f
!        Velocity(:,:, 1, iUp_, iBlock) = 2*f
     enddo

  endif

!  NDensityS(nLons/2,nLats/2,2,1:3,1) = NDensityS(nLons/2,nLats/2,2,1:3,1)*50.0
!  temperature(nLons/2,nLats/2,2,1) = temperature(nLons/2,nLats/2,2,1)*100.0

end subroutine user_perturbation


! ----------------------------------------------------------------
! If you want to output some specific variables, then do that here.
! In ModUserGITM, there are two variables defined, UserData2D and UserData3D.
! To output variables:
! 1. Figure out which variable you want to output.
! 2. Go into the code where the variable is set and copy it into
!    UserData3D or UserData2D.
! 3. Do this for each variable you want to output.
! 4. Edit output_header_user below, making a list of all the variables
!    that you have added.  Make sure you leave longitude, latitude, and
!    altitude in the list of variables.
! 5. Count the number of variables that you output (including 
!    Lon, Lat, and Alt). Change nVarsUser3d or nVarsUser2d in the 
!    subroutines towards the top of this file.
! 6. If you add more than 40 variables, you probably should check 
!    nUserOutputs in ModUserGITM.f90 and make sure that this number is
!    larger than the number of variables that you added.
! 7. Recompile and run. Debug. Repeat 7.
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!
! ----------------------------------------------------------------

subroutine set_nVarsUser3d

  use ModUserGITM

  ! Make sure to include Lat, Lon, and Alt

  nVarsUser3d = 6 + 3 + 3 ! add 3D winds & Ion Drift

  if (nVarsUser3d-3 > nUserOutputs) &
       call stop_gitm("Too many user outputs!! Increase nUserOutputs!!")

end subroutine set_nVarsUser3d

! ----------------------------------------------------------------
!
! ----------------------------------------------------------------

subroutine set_nVarsUser2d

  use ModUserGITM
  use ModSources, only:ED_N_Energies

  ! Make sure to include Lat, Lon, and Alt

  !nVarsUser2d = 10+ED_N_Energies
  nVarsUser2d = 3 + 4 + 4 + 3 & ! Neutral wind
       + 3 & !! follow Cheng's code, yuhong on 020524, add j
       + 9 & !! follow Cheng's code, yuhong on 062024, add j-comp
       + 9 & !! add winds at 110km, 270km and Vi at 400km
       + 9 & !! follow Cheng, add dB, 06/065/2024
       + 1   !! add magnetospheric source FAC, 08/28/2024

  if (nVarsUser2d-3 > nUserOutputs) &
       call stop_gitm("Too many user outputs!! Increase nUserOutputs!!")

end subroutine set_nVarsUser2d

! ----------------------------------------------------------------
!
! ----------------------------------------------------------------

subroutine set_nVarsUser1d

  use ModUserGITM

  ! Make sure to include Lat, Lon, and Alt

  nVarsUser1d = 4

  if (nVarsUser2d-3 > nUserOutputs) &
       call stop_gitm("Too many user outputs!! Increase nUserOutputs!!")

end subroutine set_nVarsUser1d

! ----------------------------------------------------------------
!
! ----------------------------------------------------------------

subroutine output_header_user(cType, iOutputUnit_)

  use ModUserGITM
  use ModSources, only:ED_Energies, ED_N_Energies

  implicit none

  character (len=5), intent(in) :: cType
  integer, intent(in)           :: iOutputUnit_
  integer :: n

  ! ------------------------------------------
  ! 3D Output Header
  ! ------------------------------------------

  if (cType(1:2) == '3D') then 

     write(iOutputUnit_,*) "NUMERICAL VALUES"
     write(iOutputUnit_,"(I7,6A)") nVarsUser3d, " nvars"
     write(iOutputUnit_,"(I7,7A)") nAlts+4, " nAltitudes"
     write(iOutputUnit_,"(I7,7A)") nLats+4, " nLatitudes"
     write(iOutputUnit_,"(I7,7A)") nLons+4, " nLongitudes"

     write(iOutputUnit_,*) ""

     write(iOutputUnit_,*) "VARIABLE LIST"
     write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
     write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
     write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"
     write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Eden"
     write(iOutputUnit_,"(I7,A1,a)")  5, " ", "Wi(up)" 
     write(iOutputUnit_,"(I7,A1,a)")  6, " ", "ExB(up)"
     !! add 3D winds & Ion Drift, yuhong, 06/20/2024
     write(iOutputUnit_,"(I7,A1,a)") 7, " ", "V!Di!N (east)"
     write(iOutputUnit_,"(I7,A1,a)") 8, " ", "V!Di!N (north)"
     write(iOutputUnit_,"(I7,A1,a)") 9, " ", "V!Di!N (up)"
     write(iOutputUnit_,"(I7,A1,a)") 10, " ", "V!Dn!N (east)"
     write(iOutputUnit_,"(I7,A1,a)") 11, " ", "V!Dn!N (north)"
     write(iOutputUnit_,"(I7,A1,a)") 12, " ", "V!Dn!N (up)"

  endif

  ! ------------------------------------------
  ! 2D Output Header
  ! ------------------------------------------

  if (cType(1:2) == '2D') then 

     write(iOutputUnit_,*) "NUMERICAL VALUES"
     write(iOutputUnit_,"(I7,6A)") nVarsUser2d, " nvars"
     write(iOutputUnit_,"(I7,7A)")     1, " nAltitudes"
     write(iOutputUnit_,"(I7,7A)") nLats, " nLatitudes"
     write(iOutputUnit_,"(I7,7A)") nLons, " nLongitudes"

     write(iOutputUnit_,*) ""
     write(iOutputUnit_,*) "NO GHOSTCELLS"
     write(iOutputUnit_,*) ""

     write(iOutputUnit_,*) "VARIABLE LIST"
     write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
     write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
     write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"

     write(iOutputUnit_,"(I7,A1,a)")  4, " ", "EFlux"
     write(iOutputUnit_,"(I7,A1,a)")  5, " ", "Potential"
     write(iOutputUnit_,"(I7,A1,a)")  6, " ", "TEC"
     write(iOutputUnit_,"(I7,A1,a)")  7, " ", "intJH"

     write(iOutputUnit_,"(I7,A1,a)")  8, " ", "Vi(East)"
     write(iOutputUnit_,"(I7,A1,a)")  9, " ", "Vi(North)"
     write(iOutputUnit_,"(I7,A1,a)")  10, " ", "Vi(Up)"  
     write(iOutputUnit_,"(I7,A1,a)")  11, " ", "Vn(East)"
     write(iOutputUnit_,"(I7,A1,a)")  12, " ", "Vn(North)"
     
     write(iOutputUnit_,"(I7,A1,a)")  13, " ", "O_N2"
     write(iOutputUnit_,"(I7,A1,a)")  14, " ", "Vi(Up_all)"

     !------------------------ wind-dynamo current ------------------
     ! follow Cheng's code, yuhong on 020524
     !---------------------------------------------------------------
     write(iOutputUnit_,"(I7,A1,a)") 15, " ", "Jp1"
     write(iOutputUnit_,"(I7,A1,a)") 16, " ", "Jp1_UcrossB"
     write(iOutputUnit_,"(I7,A1,a)") 17, " ", "Jp1_EField"
     !
     write(iOutputUnit_,"(I7,A1,a)") 18, " ", "J2D (east)"
     write(iOutputUnit_,"(I7,A1,a)") 19, " ", "J2D (north)"
     write(iOutputUnit_,"(I7,A1,a)") 20, " ", "J2D (up)"
     !
     write(iOutputUnit_,"(I7,A1,a)") 21, " ", "J2D_UcrossB (east)"
     write(iOutputUnit_,"(I7,A1,a)") 22, " ", "J2D_UcrossB (north)"
     write(iOutputUnit_,"(I7,A1,a)") 23, " ", "J2D_UcrossB (up)"
     !
     write(iOutputUnit_,"(I7,A1,a)") 24, " ", "J2D_EField (east)"
     write(iOutputUnit_,"(I7,A1,a)") 25, " ", "J2D_EField (north)"
     write(iOutputUnit_,"(I7,A1,a)") 26, " ", "J2D_EField (up)"
     !
     write(iOutputUnit_,"(I7,A1,a)") 27, " ", "dBTotal (east)"
     write(iOutputUnit_,"(I7,A1,a)") 28, " ", "dBTotal (north)"
     write(iOutputUnit_,"(I7,A1,a)") 29, " ", "dBTotal (up)"
     !
     write(iOutputUnit_,"(I7,A1,a)") 30, " ", "dBUcB (east)"
     write(iOutputUnit_,"(I7,A1,a)") 31, " ", "dBUcB (north)"
     write(iOutputUnit_,"(I7,A1,a)") 32, " ", "dBUcB (up)"
     !
     write(iOutputUnit_,"(I7,A1,a)") 33, " ", "dBEField (east)"
     write(iOutputUnit_,"(I7,A1,a)") 34, " ", "dBEField (north)"
     write(iOutputUnit_,"(I7,A1,a)") 35, " ", "dBEField (up)"
     !---------------------------------------------------------------
     write(iOutputUnit_,"(I7,A1,a)") 36, " ", "U!DE!N@110km"
     write(iOutputUnit_,"(I7,A1,a)") 37, " ", "U!DN!N@110km"
     write(iOutputUnit_,"(I7,A1,a)") 38, " ", "U!DU!N@110km"
     !
     write(iOutputUnit_,"(I7,A1,a)") 39, " ", "U!DE!N@270km"
     write(iOutputUnit_,"(I7,A1,a)") 40, " ", "U!DN!N@270km"
     write(iOutputUnit_,"(I7,A1,a)") 41, " ", "U!DU!N@270km"
     !
     write(iOutputUnit_,"(I7,A1,a)") 42, " ", "V!DE!N@400km"
     write(iOutputUnit_,"(I7,A1,a)") 43, " ", "V!DN!N@400km"
     write(iOutputUnit_,"(I7,A1,a)") 44, " ", "V!DU!N@400km"
     !
     write(iOutputUnit_,"(I7,A1,a)") 45, " ", "FAC"

     
     !write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Potential (kV)"
     !write(iOutputUnit_,"(I7,A1,a)")  5, " ", "Average Energy (keV)"
     !write(iOutputUnit_,"(I7,A1,a)")  6, " ", "Total Energy (ergs)"
     !write(iOutputUnit_,"(I7,A1,a)")  7, " ", "Discrete Average Energy (keV)"
     !write(iOutputUnit_,"(I7,A1,a)")  8, " ", "Discrete Total Energy (ergs)"
     !write(iOutputUnit_,"(I7,A1,a)")  9, " ", "Wave Average Energy (keV)"
     !write(iOutputUnit_,"(I7,A1,a)") 10, " ", "Wave Total Energy (ergs)"
     !do n=1,ED_N_Energies
     !   write(iOutputUnit_,"(I7,A6,1P,E9.3,A11)") 10+n, " Flux@",ED_energies(n), "eV (/cm2/s)"
     !enddo
  endif

  write(iOutputUnit_,*) ""

end subroutine output_header_user

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_3dUser(iBlock, iOutputUnit_)

  use ModGITM
  use ModUserGITM

  implicit none

  integer, intent(in) :: iBlock, iOutputUnit_
  integer :: iAlt, iLat, iLon, i ! yuhong, 06/20/2024

  do iAlt=-1,nAlts+2
     do iLat=-1,nLats+2
        do iLon=-1,nLons+2
           write(iOutputUnit_)       &
                Longitude(iLon,iBlock), &
                Latitude(iLat,iBlock), &
                Altitude_GB(iLon, iLat, iAlt, iBlock),&
                !UserData3D(iLon,iLat,iAlt,1:nVarsUser3d-3,iBlock)
                IDensityS(iLon,iLat,iAlt,10,iBlock), &
                Ivelocity(iLon,iLat,iAlt,3,iBlock), &
                ExB(iLon,iLat,iAlt,3), &
                (Ivelocity(iLon,iLat,iAlt,i,iBlock),i=1,3), &
                (Velocity(iLon,iLat,iAlt,i,iBlock),i=1,3)
        enddo
     enddo
  enddo

end subroutine output_3dUser

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_2dUser(iBlock, iOutputUnit_)

  use ModGITM
  use ModUserGITM

  use ModCurrents ! follow Cheng's code, yuhong on 020524
  use ModdB ! follow Cheng, yuhong 06/05/2024

  implicit none

  integer, intent(in) :: iBlock, iOutputUnit_
  integer :: iAlt, iLat, iLon, i
  integer :: iiAlt, iiLat, iiLon, iAlt1, iAlt2, iAlt3  ! follow Cheng, add dB, yuhong, 05/06/2024
  
  real :: alts_in(nAlts)

  do iLat=1,nLats
  !do iLat=-1,nLats+2 ! why using this
     iiLat = min(max(iLat,1),nLats)
     do iLon=1,nLons
     !do iLon=-1,nLons+2 ! why using this
        iiLon = min(max(iLon,1),nLons)
        !!!
        !iLonAll=LongitudeIndex(iiLon,iBlock)
        !iLatAll=LatitudeIndex(iiLat,iBlock)
        !!!
        alts_in = Altitude_GB(iLon,iLat,1:nAlts,iBlock)
        iAlt = minloc(abs(alts_in-250000.),dim=1) ! Change to 250 km
        iAlt1 = minloc(abs(alts_in-110000.),dim=1) ! Change to 110 km
        iAlt2 = minloc(abs(alts_in-270000.),dim=1) ! Change to 270 km
        iAlt3 = minloc(abs(alts_in-400000.),dim=1) ! Change to 400 km  

        write(iOutputUnit_)       &
             Longitude(iLon,iBlock), &
             Latitude(iLat,iBlock), &
             Altitude_GB(iLon, iLat, iAlt, iBlock),&
             ElectronEnergyFlux(iLon,iLat),& ! Energy flux
             Gedy_pot(iLon, iLat, iAlt),& ! potential
             scTEC(iLon,iLat), & ! TEC
             HeightIntegratedJH(iLon,iLat), & ! Height-integrated Joule heating
             (Ivelocity(iLon,iLat,iAlt,i,iBlock),i=1,2), & ! Vi
             ExB(iLon,iLat,iAlt,3), & ! ExB up
             (Velocity(iLon,iLat,iAlt,i,iBlock),i=1,2), & ! Vn
             on2ratio(iLon,iLat), & ! O/N2 ratio
             Ivelocity(iLon,iLat,iAlt,3,iBlock), & ! Vi up
             !! follow Cheng's code, yuhong on 020524
             Jp1(iLon,iLat), &
             Jp1UcB(iLon,iLat), &
             Jp1EF(iLon,iLat), &
             (J2D(iiLon,iiLat,i),i=1,3), &
             (J2DUcB(iiLon,iiLat,i),i=1,3), &
             (J2DEF(iiLon,iiLat,i),i=1,3), &
             (dBTotal(iiLon,iiLat,i,iBlock)*1.e9,i=1,3), &  !dB in nT
             (dBUcB(iiLon,iiLat,i,iBlock)*1.0e9,i=1,3), &
             (dBEF(iiLon,iiLat,i,iBlock)*1.0e9,i=1,3), &
             (Velocity(iLon,iLat,iAlt1,i,iBlock),i=1,3), &
             (Velocity(iLon,iLat,iAlt2,i,iBlock),i=1,3), &
             (Ivelocity(iLon,iLat,iAlt3,i,iBlock),i=1,3), &
             Gedy_fac(iLon, iLat, iAlt)
        
     enddo
  enddo

end subroutine output_2dUser

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_1dUser(iBlock, iOutputUnit_)

  use ModGITM
  use ModUserGITM

  implicit none

  integer, intent(in) :: iBlock, iOutputUnit_
  integer :: iAlt, iLat, iLon

  iAlt = 1
  do iLat=1,nLats
     do iLon=1,nLons
        write(iOutputUnit_)       &
             Longitude(iLon,iBlock), &
             Latitude(iLat,iBlock), &
             Altitude_GB(iLon, iLat, iAlt, iBlock),&
             UserData2D(iLon,iLat,iAlt,1:nVarsUser2d-3,iBlock)
             
     enddo
  enddo

end subroutine output_1dUser

