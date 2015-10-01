! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------
    
Program sibveg

! This code creates CCAM vegie data using the 1km SiB dataset

Implicit None

Character*80, dimension(:,:), allocatable :: options
Character*80, dimension(1:10) :: fname
Character*80 topofile
Character*80 landtypeout
Character*80 newtopofile
Integer binlimit, nopts, month
real zmin
Logical fastsib,siblsmask,ozlaipatch,usedean

usedean = .true.

Namelist/vegnml/ topofile,fastsib,                  &
                 landtypeout,siblsmask,newtopofile, &
                 binlimit,month,                    &
                 zmin,ozlaipatch,usedean

#ifndef stacklimit
! For linux only - removes stacklimit on all processors
call setstacklimit(-1)
#endif 

Write(6,*) 'SIBVEG - SiB/DG 1km to CC grid (SEP-15)'

! Read switches
nopts=1
Allocate (options(nopts,2))
options(:,1) = (/ '-s' /)
options(:,2) = ''

Call readswitch(options,nopts)
Call defaults(options,nopts)

! Read namelist
Write(6,*) 'Input &vegnml namelist'
Read(5,NML=vegnml)
Write(6,*) 'Namelist accepted'

! Generate veg data
fname(1)=topofile
fname(8)=landtypeout
fname(10)=newtopofile

Call createveg(options,nopts,fname,fastsib,siblsmask,ozlaipatch,month,binlimit,zmin,usedean)

Deallocate(options)

Stop
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine displays the help message
!

Subroutine help()

Implicit None

Write(6,*)
Write(6,*) "Usage:"
Write(6,*) "  sibveg -s size < sibveg.nml"
Write(6,*)
Write(6,*) "Options:"
Write(6,*) "  -s size      size of array used for reading SiB data"
Write(6,*) "               (typically =500).  The larger the array, the"
Write(6,*) "               faster and more accurate the output."
Write(6,*)
Write(6,*) "Namelist:"
Write(6,*) "  The namelist specifies what data to store and the filenames"
Write(6,*) "  to use.  For example:"
Write(6,*)
Write(6,*) '  &vegnml'
Write(6,*) '    month=0'
Write(6,*) '    topofile="topout"'
Write(6,*) '    newtopofile="topoutb"'
Write(6,*) '    landtypeout="veg"'
Write(6,*) '    fastsib=t'
Write(6,*) '    siblsmask=t'
Write(6,*) '    ozlaipatch=f'
Write(6,*) '    binlimit=2'
Write(6,*) '    zmin=40.'
Write(6,*) '    usedean=t'
Write(6,*) '  &end'
Write(6,*)
Write(6,*) '  where:'
Write(6,*) '    month         = the month to process (1-12, 0=all)'
Write(6,*) '    topofile      = topography (input) file'
Write(6,*) '    newtopofile   = Output topography file name'
Write(6,*) '                    (if siblsmask=t)'
Write(6,*) '    landtypeout   = Land-use classification filename'
Write(6,*) '    fastsib       = Turn on fastsib mode (see notes below)'
Write(6,*) '    siblsmask     = Define land/sea mask from SiB dataset'
Write(6,*) '    ozlaipatch    = Replace Australian LAI with Lee, S. data'
Write(6,*) '    binlimit      = The minimum ratio between the grid'
Write(6,*) '                    length scale and the length scale of'
Write(6,*) '                    the aggregated land-use data (see notes'
Write(6,*) '                    below).'
Write(6,*) '    zmin          = Reference height for blending roughness lengths'
Write(6,*) '                    (usually set to the first model level)'
Write(6,*) '    usedean       = Use 6km Australian land-use data'
Write(6,*)
Write(6,*) 'NOTES: Fastsib mode will speed up the code by aggregating'
Write(6,*) '       land-use data at a coarser resolution before'
Write(6,*) '       processing.  The degree of aggregation is determined'
Write(6,*) '       by the avaliable memory (i.e., -s switch).   Usually,'
Write(6,*) '       fastsib is used to test the output and then the'
Write(6,*) '       dataset is subsequently regenerated with fastsib=f.'
Write(6,*)
Write(6,*) '       During the binning of land-use data, the length scale'
Write(6,*) '       eventually becomes sufficently small so that binlimit'
Write(6,*) '       can no longer be satisfied.  Under these circumstances'
Write(6,*) '       the code will use the minimum length scale of the'
Write(6,*) '       SiB dataset (e.g., 1km) for all data that is'
Write(6,*) '       subsequently binned.  In the case where the grid scale'
Write(6,*) '       is less than the minimum length scale of the SiB'
Write(6,*) '       dataset, the code will interpolate instead of binning.'
Write(6,*)
Stop

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine determins the default values for the switches
!

Subroutine defaults(options,nopts)

Implicit None

Integer nopts
Character(len=*), dimension(nopts,2), intent(inout) :: options
Integer siz
Integer locate

siz=locate('-s',options(:,1),nopts)

If (options(siz,2).EQ.'') then
  options(siz,2)='500'
End if

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine processes the sib data
!

Subroutine createveg(options,nopts,fname,fastsib,siblsmask,ozlaipatch,month,binlimit,zmin,usedean)

Use ccinterp

Implicit None

Logical, intent(in) :: fastsib,siblsmask,ozlaipatch,usedean
Integer, intent(in) :: nopts,binlimit,month
real, intent(in) :: zmin
Character(len=*), dimension(1:nopts,1:2), intent(in) :: options
Character(len=*), dimension(1:10), intent(in) :: fname
Character*80, dimension(1:3) :: outputdesc
character*90 filename
Character*80 returnoption,csize,filedesc
Character*47 header
Character*9 formout
Character*2 monthout
real, dimension(:,:,:), allocatable :: laidata
Real, dimension(:,:,:), allocatable :: landdata,rawlanddata,soildata,rlld,albvisdata,albnirdata
Real, dimension(:,:), allocatable :: gridout,lsdata,urbandata,oceandata,rdata
Real, dimension(1:3,1:2) :: alonlat
Real, dimension(1:2) :: lonlat
Real, dimension(1:12) :: atime
Real, dimension(1) :: alvl
Real schmidt,dsx,ds
Integer, dimension(:,:), allocatable :: idata
Integer, dimension(1:2) :: sibdim
Integer, dimension(1:4) :: dimnum,dimid,dimcount
Integer, dimension(0:4) :: ncidarr
Integer, dimension(1:6) :: adate
Integer, dimension(1:9) :: varid
Integer sibsize,tunit,i,j,ierr,sibmax(1),mthrng
integer tt

mthrng=1
if (month==0) then
  mthrng=12
end if
if ((month<0).or.(month>12)) then
  write(6,*) "ERROR: Invalid month ",month
  stop
end if

csize=returnoption('-s',options,nopts)
Read(csize,FMT=*,IOSTAT=ierr) sibsize
If (ierr.NE.0) then
  Write(6,*) 'ERROR: Invalid array size.  Must be an integer.'
  Stop
End if

! Read topography file
tunit=1
Call readtopography(tunit,fname(1),sibdim,lonlat,schmidt,dsx,header)
Write(6,*) "Dimension : ",sibdim
Write(6,*) "lon0,lat0 : ",lonlat
Write(6,*) "Schmidt   : ",schmidt

Allocate(gridout(1:sibdim(1),1:sibdim(2)),rlld(1:sibdim(1),1:sibdim(2),1:2))
Allocate(rawlanddata(1:sibdim(1),1:sibdim(2),0:81),albvisdata(1:sibdim(1),1:sibdim(2),0:mthrng-1))
Allocate(laidata(1:sibdim(1),1:sibdim(2),0:mthrng-1),albnirdata(1:sibdim(1),1:sibdim(2),0:mthrng-1))
Allocate(soildata(1:sibdim(1),1:sibdim(2),0:8),lsdata(1:sibdim(1),1:sibdim(2)))
Allocate(urbandata(1:sibdim(1),1:sibdim(2)),landdata(1:sibdim(1),1:sibdim(2),0:42))
Allocate(oceandata(1:sibdim(1),1:sibdim(2)))

! Determine lat/lon to CC mapping
Call ccgetgrid(rlld,gridout,sibdim,lonlat,schmidt,ds)

! Read sib data
Call getdata(rawlanddata,lonlat,gridout,rlld,sibdim,81,      sibsize,'land',  fastsib,ozlaipatch,binlimit,month,usedean)
Call getdata(soildata,   lonlat,gridout,rlld,sibdim,8,       sibsize,'soil',  fastsib,ozlaipatch,binlimit,month,usedean)
Call getdata(laidata,    lonlat,gridout,rlld,sibdim,mthrng-1,sibsize,'lai',   fastsib,ozlaipatch,binlimit,month,usedean)
Call getdata(albvisdata, lonlat,gridout,rlld,sibdim,mthrng-1,sibsize,'albvis',fastsib,ozlaipatch,binlimit,month,usedean)
Call getdata(albnirdata, lonlat,gridout,rlld,sibdim,mthrng-1,sibsize,'albnir',fastsib,ozlaipatch,binlimit,month,usedean)

deallocate(gridout)

write(6,*) "Preparing data..."
! extract appended urban data
urbandata(:,:)=sum(rawlanddata(:,:,30:50),3) !+rawlanddata(:,:,81) ! neglect Dean's urban

! remove SiB classes 14-50 that contain urban
Call sibfix(landdata,rawlanddata,rlld,sibdim)

deallocate(rawlanddata)

if (siblsmask) then
  write(6,*) "Using SiB land/sea mask"
  ! extract ocean data
  where ((landdata(:,:,0)+landdata(:,:,42))>0.)
    oceandata=landdata(:,:,0)/(landdata(:,:,0)+landdata(:,:,42))
  elsewhere
    oceandata=0.
  end where
  lsdata=real(nint(landdata(:,:,0)+landdata(:,:,42)))
  call cleantopo(tunit,fname(1),fname(10),lsdata,oceandata,sibdim)
else
  write(6,*) "Using topography land/sea mask"
  call gettopols(tunit,fname(1),lsdata,sibdim)
end if

urbandata=min(urbandata,(1.-lsdata))

! Clean-up soil, lai, veg, albedo and urban data
Call cleansib(landdata,lsdata,rlld,sibdim)
Call cleanreal(soildata,8,lsdata,rlld,sibdim)
Call cleanreal(laidata,mthrng-1,lsdata,rlld,sibdim)
Call cleanreal(albvisdata,mthrng-1,lsdata,rlld,sibdim)
Call cleanreal(albnirdata,mthrng-1,lsdata,rlld,sibdim)
do i=1,sibdim(1)
  do j=1,sibdim(2)
    if (lsdata(i,j)>=0.5) then
      albvisdata(i,j,:)=0.08 ! 0.07 in Masson (2003)
      albnirdata(i,j,:)=0.08 ! 0.20 in Masson (2003)
    end if
  end do
end do

Deallocate(rlld,oceandata)
Allocate(idata(1:sibdim(1),1:sibdim(2)),rdata(1:sibdim(1),1:sibdim(2)))

! Prep nc output
dimnum(1:2)=sibdim(1:2) ! CC grid dimensions
dimnum(3)=1 ! Turn off level
dimnum(4)=1 ! Number of months in a year
adate=0 ! Turn off date
adate(2)=1 ! time units=months

do tt=1,mthrng
  write(6,*) "Writing month ",tt,"/",mthrng
  
  if (mthrng==1) then
    filename=fname(8)
  else
    write(filename,"(A,'.',I2.2)") trim(fname(8)),tt
  end if

  Call ncinitcc(ncidarr,filename,dimnum(1:3),dimid,adate)
  outputdesc(1)='soil'
  outputdesc(2)='Soil classification'
  outputdesc(3)='none'
  Call ncaddvargen(ncidarr,outputdesc,5,2,varid(2),1.,0.)
  outputdesc(1)='albvis'
  outputdesc(2)='Soil albedo (VIS)'
  outputdesc(3)=''
  Call ncaddvargen(ncidarr,outputdesc,5,3,varid(3),1.,0.)
  outputdesc(1)='albnir'
  outputdesc(2)='Soil albedo (NIR)'
  outputdesc(3)=''
  Call ncaddvargen(ncidarr,outputdesc,5,3,varid(4),1.,0.)
  outputdesc(1)='rsmin'
  outputdesc(2)='RSmin'
  outputdesc(3)=''
  Call ncaddvargen(ncidarr,outputdesc,5,3,varid(5),1.,0.)
  outputdesc(1)='rough'
  outputdesc(2)='Roughness length'
  outputdesc(3)='m'
  Call ncaddvargen(ncidarr,outputdesc,5,3,varid(6),1.,0.)
  outputdesc(1)='lai'
  outputdesc(2)='Leaf Area Index'
  outputdesc(3)=''
  Call ncaddvargen(ncidarr,outputdesc,5,3,varid(7),1.,0.)
  outputdesc(1)='landtype'
  outputdesc(2)='Land-use classification'
  outputdesc(3)='none'
  Call ncaddvargen(ncidarr,outputdesc,5,2,varid(8),1.,0.)
  outputdesc(1)='urban'
  outputdesc(2)='Urban fraction'
  outputdesc(3)='none'
  Call ncaddvargen(ncidarr,outputdesc,5,2,varid(9),1.,0.)

  call ncatt(ncidarr,'lon0',lonlat(1))
  call ncatt(ncidarr,'lat0',lonlat(2))
  call ncatt(ncidarr,'schmidt',schmidt)
  call ncatt(ncidarr,'sibvegversion',2015.)

  Call ncenddef(ncidarr)
  alonlat(:,1)=(/ 1., real(sibdim(1)), 1. /)
  alonlat(:,2)=(/ 1., real(sibdim(2)), 1. /)
  alvl=1.
  if (mthrng==12) then
    atime(1)=Real(tt) ! Define Months
  else
    atime(1)=real(month)
  end if
  Call nclonlatgen(ncidarr,dimid,alonlat,alvl,atime,dimnum)


  ! Write soil type
  Write(6,*) 'Write soil type.'
  Call calsoilnear(landdata,soildata,lsdata,sibdim,idata)
  dimcount=(/ sibdim(1), sibdim(2), 1, 1 /)
  rdata(:,:)=real(idata)
  Call ncwritedatgen(ncidarr,rdata(:,:),dimcount,varid(2))

  ! Write albedo file
  Write(6,*) 'Write albedo.'
  dimcount=(/ sibdim(1), sibdim(2), 1, 1 /)
  Call ncwritedatgen(ncidarr,albvisdata(:,:,tt-1),dimcount,varid(3))
  dimcount=(/ sibdim(1), sibdim(2), 1, 1 /)
  Call ncwritedatgen(ncidarr,albnirdata(:,:,tt-1),dimcount,varid(4))

  ! Write rsmin file
  Write(6,*) 'Write rsmin.'
  Call calrsmin(landdata,laidata(:,:,tt-1),sibdim,1,rdata)
  Where(rdata.GT.995.)
    rdata=995.
  End Where
  dimcount=(/ sibdim(1), sibdim(2), 1, 1 /)
  Call ncwritedatgen(ncidarr,rdata,dimcount,varid(5))

  ! Write roughness file
  Write(6,*) 'Write roughness files.'
  Call calrough(landdata,laidata(:,:,tt-1),sibdim,1,rdata,zmin)
  Where(rdata.LT.0.01)
    rdata=0.01
  End Where
  dimcount=(/ sibdim(1), sibdim(2), 1, 1 /)
  Call ncwritedatgen(ncidarr,rdata,dimcount,varid(6))

  ! Write lai file
  Write(6,*) 'Write lai files.'
  dimcount=(/ sibdim(1), sibdim(2), 1, 1 /)
  Call ncwritedatgen(ncidarr,laidata(:,:,tt-1),dimcount,varid(7))

  ! SiB type
  Write(6,*) 'Write land-use type'
  Do i=1,sibdim(1)
    Do j=1,sibdim(2)
      if (lsdata(i,j)>=0.5) then
        idata(i,j)=0
      else
        sibmax=Maxloc(landdata(i,j,1:41))
        idata(i,j)=sibmax(1)
      end if
    End do
  End do
  where (idata>13)
    idata=idata+101-13
  end where
  where (idata>0.and.idata<=13)
    idata=idata+31
  end where
  where (idata>100)
    idata=idata-101
  end where
  dimcount=(/ sibdim(1), sibdim(2), 1, 1 /)
  rdata(:,:)=real(idata)
  Call ncwritedatgen(ncidarr,rdata(:,:),dimcount,varid(8))

  ! Urban
  Write(6,*) 'Write urban fraction'
  dimcount=(/ sibdim(1), sibdim(2), 1, 1 /)
  Call ncwritedatgen(ncidarr,urbandata,dimcount,varid(9))

  Call ncclose(ncidarr)
  
end do

Deallocate(landdata,soildata,albvisdata,albnirdata,idata,rdata,urbandata,laidata,lsdata)

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fix SiB data
!

subroutine sibfix(landdata,rawdata,rlld,sibdim)

implicit none

integer, dimension(1:2), intent(in) :: sibdim
real, dimension(1:sibdim(1),1:sibdim(2),0:81), intent(in) :: rawdata
real, dimension(sibdim(1),sibdim(2),1:2), intent(in) :: rlld
real, dimension(1:sibdim(1),1:sibdim(2),0:42), intent(out) :: landdata
real, dimension(1:sibdim(1),1:sibdim(2),0:42) :: landtemp
logical, dimension(1:sibdim(1),1:sibdim(2)) :: sermsk
integer i,ilon,ilat,pxy(2)
real nsum,wsum

landtemp(:,:,0:22)=rawdata(:,:,0:22)
do i=31,50 ! remove urban attached to SiB
  sermsk=rlld(:,:,1)>112..and.rlld(:,:,1)<154.5.and.rlld(:,:,2)>-45..and.rlld(:,:,2)<-10.
  where (.not.sermsk)
    landtemp(:,:,i-30)=landtemp(:,:,i-30)+rawdata(:,:,i)
  end where
end do

landdata(:,:,1:12) =landtemp(:,:,1:12)                 ! SIB's
landdata(:,:,0)    =landtemp(:,:,19)                   ! redefine ocean
landdata(:,:,13)   =landtemp(:,:,20)                   ! redefine ice
landdata(:,:,14:41)=rawdata(:,:,51:78)                 ! Dean's
landdata(:,:,42)   =rawdata(:,:,79)+rawdata(:,:,80)    ! redefine lake

landtemp(:,:,0:42) =landdata(:,:,0:42)
sermsk=sum(landdata(:,:,1:41),3)>0.
if (.not.any(sermsk)) return

call fill_cc_a(landtemp(:,:,1:41),sibdim(1),41,sermsk)

do ilat=1,sibdim(2)
  do ilon=1,sibdim(1)
    wsum=landdata(ilon,ilat,0)+landdata(ilon,ilat,42) ! water
    if (wsum<1.) then
      if (.not.sermsk(ilon,ilat)) then
        landdata(ilon,ilat,1:41)=landtemp(ilon,ilat,1:41)
      end if
      nsum=sum(landdata(ilon,ilat,1:41)) ! land
      landdata(ilon,ilat,1:41)=landdata(ilon,ilat,1:41)*max(1.-wsum,0.)/nsum
    end if
  end do
  if ( mod(ilat,100)==0 .or. ilat==sibdim(2) ) then
    write(6,*) "Searching ",ilat,"/",sibdim(2)
  end if
end do

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! clean land

subroutine cleansib(dataout,lsdata,rlld,sibdim)

implicit none

integer, dimension(1:2), intent(in) :: sibdim
real, dimension(1:sibdim(1),1:sibdim(2),0:42), intent(inout) :: dataout
real, dimension(1:sibdim(1),1:sibdim(2)), intent(in) :: lsdata
real, dimension(1:sibdim(1),1:sibdim(2),1:2), intent(in) :: rlld
real, dimension(1:sibdim(1),1:sibdim(2),0:42) :: datain
logical, dimension(1:sibdim(1),1:sibdim(2)) :: sermsk,ocnmsk
integer ilon,ilat,pxy(2)
real nsum,wsum

datain=dataout
sermsk=sum(datain(:,:,1:41),3).gt.0.
ocnmsk=(datain(:,:,0)+datain(:,:,42))>0.
if (.not.any(sermsk)) then
  dataout(:,:,0)=1.
  dataout(:,:,1:)=0.
  return
end if

do ilat=1,sibdim(2)
  do ilon=1,sibdim(1)
    if (lsdata(ilon,ilat)<0.5) then
      if (.not.sermsk(ilon,ilat)) then
        call findnear(pxy,ilon,ilat,sermsk,rlld,sibdim)
        dataout(ilon,ilat,1:41)=datain(pxy(1),pxy(2),1:41)
      end if
      nsum=sum(dataout(ilon,ilat,1:41))
      dataout(ilon,ilat,1:41)=dataout(ilon,ilat,1:41)*(1.-lsdata(ilon,ilat))/nsum
    else
      dataout(ilon,ilat,1:41)=0.
    end if
    if (lsdata(ilon,ilat)>=0.5) then
      if (.not.ocnmsk(ilon,ilat)) then
        call findnear(pxy,ilon,ilat,ocnmsk,rlld,sibdim)
        dataout(ilon,ilat,0)=datain(pxy(1),pxy(2),0)
        dataout(ilon,ilat,42)=datain(pxy(1),pxy(2),42)
      end if
      wsum=dataout(ilon,ilat,0)+dataout(ilon,ilat,42)
      dataout(ilon,ilat,0)=dataout(ilon,ilat,0)*lsdata(ilon,ilat)/wsum
      dataout(ilon,ilat,42)=dataout(ilon,ilat,42)*lsdata(ilon,ilat)/wsum
    else
      dataout(ilon,ilat,0)=0.
      dataout(ilon,ilat,42)=0.
    end if
  end do
end do

return
end

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! clean soil data
!

subroutine cleandata(dataout,num,lsdata,rlld,sibdim)

implicit none

integer, intent(in) :: num
integer, dimension(1:2), intent(in) :: sibdim
real, dimension(1:sibdim(1),1:sibdim(2),0:num), intent(inout) :: dataout
real, dimension(1:sibdim(1),1:sibdim(2)), intent(in) :: lsdata
real, dimension(1:sibdim(1),1:sibdim(2),1:2), intent(in) :: rlld
real, dimension(1:sibdim(1),1:sibdim(2),0:num) :: datain
logical, dimension(1:sibdim(1),1:sibdim(2)) :: sermsk
integer ilon,ilat,pxy(2)
real nsum

datain=dataout
sermsk=sum(datain(:,:,1:num),3).gt.0.
if (.not.any(sermsk)) return

do ilon=1,sibdim(1)
  do ilat=1,sibdim(2)
    if (1-nint(lsdata(ilon,ilat)).eq.1) then
      if (.not.sermsk(ilon,ilat)) then
        call findnear(pxy,ilon,ilat,sermsk,rlld,sibdim)
        dataout(ilon,ilat,:)=datain(pxy(1),pxy(2),:)
      end if
      nsum=sum(dataout(ilon,ilat,1:num))
      dataout(ilon,ilat,1:num)=dataout(ilon,ilat,1:num)*(1.-lsdata(ilon,ilat))/nsum
      dataout(ilon,ilat,0)=lsdata(ilon,ilat)
    else
      dataout(ilon,ilat,1:num)=0.
      dataout(ilon,ilat,0)=1.
    end if
  end do
end do

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! clean real data
!

subroutine cleanreal(dataout,num,lsdata,rlld,sibdim)

implicit none

integer, intent(in) :: num
integer, dimension(1:2), intent(in) :: sibdim
real, dimension(1:sibdim(1),1:sibdim(2),0:num), intent(inout) :: dataout
real, dimension(1:sibdim(1),1:sibdim(2)), intent(in) :: lsdata
real, dimension(1:sibdim(1),1:sibdim(2),1:2), intent(in) :: rlld
real, dimension(1:sibdim(1),1:sibdim(2),0:num) :: datain
logical, dimension(1:sibdim(1),1:sibdim(2)) :: sermsk
integer ilon,ilat,pxy(2)
real nsum

datain=dataout

sermsk=.true.
do ilon=0,num
  sermsk=sermsk.and.(datain(:,:,ilon)>0.)
end do
if (.not.any(sermsk)) return

do ilon=1,sibdim(1)
  do ilat=1,sibdim(2)
    if ((1-nint(lsdata(ilon,ilat))==1).and.(.not.sermsk(ilon,ilat))) then
      call findnear(pxy,ilon,ilat,sermsk,rlld,sibdim)
      dataout(ilon,ilat,:)=datain(pxy(1),pxy(2),:)
    end if
  end do
end do

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine rewrites the land/sea mask in the topography
! data file.
!

Subroutine cleantopo(topounit,toponame,topoout,lsmskin,oceanin,sibdim)

use netcdf_m

Implicit None

Integer, intent(in) :: topounit
Integer, dimension(2), intent(in) :: sibdim
integer, dimension(3) :: spos,npos
integer, dimension(3) :: dimid
Integer ilout,ierr,ia,ib,i
integer ncid,lnctopo,varid
Character(len=*), intent(in) :: toponame,topoout
Character*80 formout
Character*47 dc
Real, dimension(sibdim(1),sibdim(2)), intent(in) :: lsmskin,oceanin
Real, dimension(sibdim(1),sibdim(2)) :: topo,sd,lsmsk
real, dimension(sibdim(2)) :: dum
Real, dimension(1) :: ra,rb,rc,rd
ilout=Min(sibdim(1),30) ! To be compatiable with terread

Write(6,*) "Adjust topography data for consistancy with land-sea mask"

ierr=nf_open(toponame,nf_nowrite,ncid)
if (ierr==0) then
  lnctopo=1
  spos=1
  npos(1)=sibdim(1)
  npos(2)=sibdim(2)
  npos(3)=1
  ierr=nf_get_att_real(ncid,nf_global,'lon0',ra(1))
  ierr=nf_get_att_real(ncid,nf_global,'lat0',rb(1))
  ierr=nf_get_att_real(ncid,nf_global,'schmidt',rc(1))
  ierr=nf_inq_varid(ncid,'zs',varid)
  ierr=nf_get_vara_real(ncid,varid,spos,npos,topo)
  ierr=nf_inq_varid(ncid,'lsm',varid)
  ierr=nf_get_vara_real(ncid,varid,spos,npos,lsmsk)
  ierr=nf_inq_varid(ncid,'tsd',varid)
  ierr=nf_get_vara_real(ncid,varid,spos,npos,sd)
  ierr=nf_close(ncid)
else
  lnctopo=0
  Open(topounit,FILE=toponame,FORM='formatted',STATUS='old',IOSTAT=ierr)
  Read(topounit,*,IOSTAT=ierr) ia,ib,ra,rb,rc,rd,dc
  Read(topounit,*,IOSTAT=ierr) topo ! Topography data
  Read(topounit,*,IOSTAT=ierr) lsmsk ! land/sea mask (to be replaced)
  Read(topounit,*,IOSTAT=ierr) sd ! Topography standard deviation
  Close(topounit)
end if

If (ierr.NE.0) then
  Write(6,*) "ERROR: Cannot read file ",trim(toponame)
  Stop
End if

! MJT notes - SiB merges both lakes and ocean together, so we
! cannot split them here anymore.

lsmsk=Real(1-nint(lsmskin))
where ((nint(oceanin)==1).and.(nint(lsmskin)==1))
!  topo(:,:)=0.
  sd(:,:)=0.
end where

if (lnctopo==1) then
  ierr=nf_create(topoout,nf_clobber,ncid)
  if (ierr/=0) then
    write(6,*) "ERROR creating output topography file ",ierr
    stop
  end if
  ierr=nf_def_dim(ncid,'longitude',sibdim(1),dimid(1))
  ierr=nf_def_dim(ncid,'latitude',sibdim(2),dimid(2))
  ierr=nf_def_dim(ncid,'time',nf_unlimited,dimid(3))
  ierr=nf_def_var(ncid,'longitude',nf_float,1,dimid(1),varid)
  ierr=nf_def_var(ncid,'latitude',nf_float,1,dimid(2),varid)
  ierr=nf_def_var(ncid,'time',nf_float,1,dimid(3),varid)
  ierr=nf_def_var(ncid,'zs',nf_float,3,dimid(1:3),varid)
  ierr=nf_def_var(ncid,'lsm',nf_float,3,dimid(1:3),varid)
  ierr=nf_def_var(ncid,'tsd',nf_float,3,dimid(1:3),varid)
  ierr=nf_put_att_real(ncid,nf_global,'lon0',nf_real,1,ra)
  ierr=nf_put_att_real(ncid,nf_global,'lat0',nf_real,1,rb)
  ierr=nf_put_att_real(ncid,nf_global,'schmidt',nf_real,1,rc)
  ierr=nf_enddef(ncid)
  do i=1,sibdim(2)
    dum(i)=real(i)
  end do
  ierr=nf_inq_varid(ncid,'longitude',varid)
  ierr=nf_put_vara_real(ncid,varid,spos(1:1),npos(1:1),dum(1:sibdim(1)))
  ierr=nf_inq_varid(ncid,'latitude',varid)
  ierr=nf_put_vara_real(ncid,varid,spos(2:2),npos(2:2),dum(1:sibdim(2)))
  ierr=nf_inq_varid(ncid,'time',varid)
  dum(1)=0.
  ierr=nf_put_vara_real(ncid,varid,spos(3:3),npos(3:3),dum(1:1))
  ierr=nf_inq_varid(ncid,'zs',varid)
  ierr=nf_put_vara_real(ncid,varid,spos,npos,topo)
  ierr=nf_inq_varid(ncid,'lsm',varid)
  ierr=nf_put_vara_real(ncid,varid,spos,npos,lsmsk)
  ierr=nf_inq_varid(ncid,'tsd',varid)
  ierr=nf_put_vara_real(ncid,varid,spos,npos,sd)
  ierr=nf_close(ncid)
else
  Open(topounit,FILE=topoout,FORM='formatted',STATUS='replace',IOSTAT=ierr)
  Write(topounit,'(i4,i6,2f10.3,f6.3,f10.0," ",a39)',IOSTAT=ierr) ia,ib,ra,rb,rc,rd,dc
  Write(formout,'("(",i3,"f7.0)")',IOSTAT=ierr) ilout
  Write(topounit,formout,IOSTAT=ierr) topo ! Topography data
  Write(formout,'("(",i3,"f4.1)")',IOSTAT=ierr) ilout
  Write(topounit,formout,IOSTAT=ierr) lsmsk ! land/sea mask
  Write(formout,'("(",i3,"f6.0)")',IOSTAT=ierr) ilout
  Write(topounit,formout,IOSTAT=ierr) sd ! Topography standard deviation
  Close(topounit)
end if

If (ierr.NE.0) then
  Write(6,*) "ERROR: Cannot write file ",trim(topoout)
  Stop
End if

Return
End

