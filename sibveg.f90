Program sibveg

! This code creates CCAM vegie data using the 1km SiB dataset

Implicit None

Character*80, dimension(:,:), allocatable :: options
Character*80, dimension(1:11) :: fname
Character*80 topofile,albvisout,albnirout,rsminout
Character*80 soilout,landtypeout,laiout,roughout
Character*80 newtopofile,urbanout,soilmethod
Integer binlimit, nopts, month
real zmin
Logical fastsib,siblsmask,ozlaipatch

Namelist/vegnml/ topofile,albvisout,albnirout,rsminout, &
                 soilout,fastsib,laiout,roughout, &
                 landtypeout,siblsmask,newtopofile, &
                 binlimit,urbanout,soilmethod,month, &
		 zmin,ozlaipatch

Write(6,*) 'SIBVEG - SiB 1km to CC grid (AUG-08)'

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
fname(2)=soilout
fname(3)=albvisout
fname(4)=albnirout
fname(5)=rsminout
fname(6)=roughout
fname(7)=laiout
fname(8)=landtypeout
fname(9)=urbanout
fname(10)=newtopofile
fname(11)=soilmethod

Call createveg(options,nopts,fname,fastsib,siblsmask,ozlaipatch,month,binlimit,zmin)

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
Write(6,*) '    soilout="soil"'
Write(6,*) '    albvisout="albvis"'
Write(6,*) '    albnirout="albnir"'
Write(6,*) '    rsminout="rsmin"'
Write(6,*) '    roughout="rough"'
Write(6,*) '    laiout="lai"'
Write(6,*) '    urbanout="urban"'
Write(6,*) '    landtypeout="veg"'
Write(6,*) '    fastsib=t'
Write(6,*) '    siblsmask=t'
Write(6,*) '    ozlaipatch=f'
Write(6,*) '    soilmethod="near"'
Write(6,*) '    binlimit=2'
Write(6,*) '    zmin=40.'
Write(6,*) '  &end'
Write(6,*)
Write(6,*) '  where:'
Write(6,*) '    month         = the month to process (1-12, 0=all)'
Write(6,*) '    topofile      = topography (input) file'
Write(6,*) '    newtopofile   = Output topography file name'
Write(6,*) '                    (if siblsmask=t)'
Write(6,*) '    soilout       = Soil filename (Zobler)'
Write(6,*) '    albvisout     = Albedo (VIS) filename'
Write(6,*) '    albnirout     = Albedo (NIR) filename'
Write(6,*) '    rsminout      = RSmin filename'
Write(6,*) '    roughout      = Roughness filename'
Write(6,*) '    laiout        = Leaf Area Index filename'
Write(6,*) '    urbanout      = Urban cover fraction filename'
Write(6,*) '    landtypeout   = Land-use classification filename'
Write(6,*) '    fastsib       = Turn on fastsib mode (see notes below)'
Write(6,*) '    siblsmask     = Define land/sea mask from SiB dataset'
Write(6,*) '    ozlaipatch    = Replace Australian LAI with Lee, S. data'
Write(6,*) '    soilmethod    = Method to use for soil classification'
Write(6,*) '                    Valid methods are:'
Write(6,*) '                    usda   - Use Zobler soil classification'
Write(6,*) '                             via USDA types (8 types, 9=Ice)'
Write(6,*) '                    near   - Use Zobler soil classification'
Write(6,*) '                             via nearest sand/clay fraction'
Write(6,*) '                             (8 types, 9=Ice)'
Write(6,*) '    binlimit      = The minimum ratio between the grid'
Write(6,*) '                    length scale and the length scale of'
Write(6,*) '                    the aggregated land-use data (see notes'
Write(6,*) '                    below).'
Write(6,*) '    zmin          = Reference height for blending roughness lengths'
Write(6,*) '                    (usually set to the first model level)'
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

Subroutine createveg(options,nopts,fname,fastsib,siblsmask,ozlaipatch,month,binlimit,zmin)

Use ccinterp

Implicit None

Logical, intent(in) :: fastsib,siblsmask,ozlaipatch
Integer, intent(in) :: nopts,binlimit,month
real, intent(in) :: zmin
Character(len=*), dimension(1:nopts,1:2), intent(in) :: options
Character(len=*), dimension(1:11), intent(in) :: fname
Character*80, dimension(1:3) :: outputdesc
Character*80 returnoption,csize,filedesc
Character*45 header
Character*9 formout
Character*2 monthout
real, dimension(:,:,:), allocatable :: laidata
Real, dimension(:,:,:), allocatable :: landdata,rawlanddata,soildata,rlld,albvisdata,albnirdata,rdata
Real, dimension(:,:), allocatable :: gridout,lsdata,urbandata,oceandata
Real, dimension(1:3,1:2) :: alonlat
Real, dimension(1:2) :: lonlat
Real, dimension(1:12) :: atime
Real, dimension(1) :: alvl
Real schmidt,dsx,ds,urbanfrac
Integer, dimension(:,:), allocatable :: idata
Integer, dimension(1:2) :: sibdim
Integer, dimension(1:4) :: dimnum,dimid,dimcount
Integer, dimension(0:4) :: ncidarr
Integer, dimension(1:6) :: adate
Integer, dimension(1:9) :: varid
Integer sibsize,tunit,i,j,ierr,sibmax(1),mthrng

mthrng=1
if (month.eq.0) then
  mthrng=12
end if
if ((month.lt.0).or.(month.gt.12)) then
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
Allocate(rawlanddata(1:sibdim(1),1:sibdim(2),0:50),albvisdata(1:sibdim(1),1:sibdim(2),0:mthrng-1))
Allocate(laidata(1:sibdim(1),1:sibdim(2),0:mthrng-1),albnirdata(1:sibdim(1),1:sibdim(2),0:mthrng-1))
Allocate(soildata(1:sibdim(1),1:sibdim(2),0:1),lsdata(1:sibdim(1),1:sibdim(2)))
Allocate(urbandata(1:sibdim(1),1:sibdim(2)),landdata(1:sibdim(1),1:sibdim(2),0:13))
Allocate(oceandata(1:sibdim(1),1:sibdim(2)))

! Determine lat/lon to CC mapping
Call ccgetgrid(rlld,gridout,sibdim,lonlat,schmidt,ds)

! Read sib data
Call getdata(rawlanddata,lonlat,gridout,rlld,sibdim,50,sibsize,'land',fastsib,ozlaipatch,binlimit,month)
Call getdata(soildata,lonlat,gridout,rlld,sibdim,1,sibsize,'soil',fastsib,ozlaipatch,binlimit,month)
Call getdata(laidata,lonlat,gridout,rlld,sibdim,mthrng-1,sibsize,'lai',fastsib,ozlaipatch,binlimit,month)
Call getdata(albvisdata,lonlat,gridout,rlld,sibdim,mthrng-1,sibsize,'albvis',fastsib,ozlaipatch,binlimit,month)
Call getdata(albnirdata,lonlat,gridout,rlld,sibdim,mthrng-1,sibsize,'albnir',fastsib,ozlaipatch,binlimit,month)

write(6,*) "Preparing data..."
! extract appended urban data
urbandata(:,:)=sum(rawlanddata(:,:,30:50),3)
! extract ocean data
oceandata(:,:)=rawlanddata(:,:,19)
! remove SiB classes >13
Call sibfix(landdata,rawlanddata,rlld,sibdim)

if (siblsmask) then
  write(6,*) "Using SiB land/sea mask"
  lsdata=landdata(:,:,0)
  call cleantopo(tunit,fname(1),fname(10),lsdata,oceandata,sibdim)
else
  write(6,*) "Using topography land/sea mask"
  call gettopols(tunit,fname(1),lsdata,sibdim)
end if

urbandata=min(urbandata,(1.-lsdata))

! Clean-up soil, lai, veg, albedo and urban data
Call cleandata(landdata,13,lsdata,rlld,sibdim)
Call cleanreal(soildata,1,lsdata,rlld,sibdim)
Call cleanreal(laidata,mthrng-1,lsdata,rlld,sibdim)
Call cleanreal(albvisdata,mthrng-1,lsdata,rlld,sibdim)
Call cleanreal(albnirdata,mthrng-1,lsdata,rlld,sibdim)
do i=1,sibdim(1)
  do j=1,sibdim(2)
    if (lsdata(i,j).ge.0.5) then
      albvisdata(i,j,:)=0.08 ! 0.07 in Masson (2003)
      albnirdata(i,j,:)=0.08 ! 0.20 in Masson (2003)
    end if
  end do
end do

Deallocate(gridout,rlld,rawlanddata,oceandata)
Allocate(idata(1:sibdim(1),1:sibdim(2)),rdata(1:sibdim(1),1:sibdim(2),1:mthrng))

! Prep nc output
dimnum(1:2)=sibdim(1:2) ! CC grid dimensions
dimnum(3)=1 ! Turn off level
dimnum(4)=mthrng ! Number of months in a year
adate=0 ! Turn off date
adate(2)=1 ! time units=months
Call ncinitcc(ncidarr,'veg.nc',dimnum(1:3),dimid,adate)
outputdesc=(/ 'soil', 'Soil classification', 'none' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(2),1.,0.)
outputdesc=(/ 'albvis', 'Albedo (VIS)', '' /)
Call ncaddvargen(ncidarr,outputdesc,5,3,varid(3),1.,0.)
outputdesc=(/ 'albnir', 'Albedo (NIR)', '' /)
Call ncaddvargen(ncidarr,outputdesc,5,3,varid(4),1.,0.)
outputdesc=(/ 'rsmin', 'RSmin', '' /)
Call ncaddvargen(ncidarr,outputdesc,5,3,varid(5),1.,0.)
outputdesc=(/ 'rough', 'Roughness length', 'm' /)
Call ncaddvargen(ncidarr,outputdesc,5,3,varid(6),1.,0.)
outputdesc=(/ 'lai', 'Leaf Area Index', '' /)
Call ncaddvargen(ncidarr,outputdesc,5,3,varid(7),1.,0.)
outputdesc=(/ 'landtype', 'Land-use classification', 'none' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(8),1.,0.)
outputdesc=(/ 'urban', 'Urban fraction', 'none' /)
Call ncaddvargen(ncidarr,outputdesc,5,2,varid(9),1.,0.)
Call ncenddef(ncidarr)
alonlat(:,1)=(/ 1., real(sibdim(1)), 1. /)
alonlat(:,2)=(/ 1., real(sibdim(2)), 1. /)
alvl=1.
if (mthrng.eq.12) then
  Do i=1,12
    atime(i)=Real(i) ! Define Months
  End do
else
  atime(1)=real(month)
end if
Call nclonlatgen(ncidarr,dimid,alonlat,alvl,atime,dimnum)


! Write soil type
Write(6,*) 'Write soil type file.'
select case(fname(11))
  case('near')
    Call calsoilnear(landdata,soildata,lsdata,sibdim,idata)
  case('usda')
    Call calsoilusda(landdata,soildata,lsdata,sibdim,idata)
  case DEFAULT
    write(6,*) "ERROR: Unknown soil method ",trim(fname(11))
    stop
end select
Write(formout,'(1h(,i3,2hi3,1h))') sibdim(1)
Open(1,File=fname(2))
Write(1,'(i3,i4,2f8.3,f6.3,f8.0," ",a39)') sibdim(1),sibdim(2),lonlat(1),lonlat(2),schmidt,ds,'soil'
Write(1,formout) idata
Close(1)
dimcount=(/ sibdim(1), sibdim(2), 1, 1 /)
Call ncwritedatgen(ncidarr,Real(idata),dimcount,varid(2))

! Write albedo file
Write(6,*) 'Write albedo files.'
Write(formout,'("(",i3,"f4.0)" )') sibdim(1)
if (mthrng.eq.12) then
  Do i=1,12
    Write(monthout,'(I2.2)') i
    filedesc=trim(fname(3))//'.'//trim(monthout)
    Open(1,File=filedesc)
    Write(1,'(i3,i4,2f8.3,f6.3,f8.0," ",a39)') sibdim(1),sibdim(2),lonlat(1),lonlat(2),schmidt,ds,'albedo'
    Write(1,formout) albvisdata(:,:,i-1)*100.
    Close(1)
  End Do
else
  Open(1,File=fname(3))
  Write(1,'(i3,i4,2f8.3,f6.3,f8.0," ",a39)') sibdim(1),sibdim(2),lonlat(1),lonlat(2),schmidt,ds,'albedo'
  Write(1,formout) albvisdata(:,:,0)*100.
  Close(1)
end if
dimcount=(/ sibdim(1), sibdim(2), mthrng, 1 /)
Call ncwritedatgen(ncidarr,albvisdata,dimcount,varid(3))

if (mthrng.eq.12) then
  Do i=1,12
    Write(monthout,'(I2.2)') i
    filedesc=trim(fname(4))//'.'//trim(monthout)
    Open(1,File=filedesc)
    Write(1,'(i3,i4,2f8.3,f6.3,f8.0," ",a39)') sibdim(1),sibdim(2),lonlat(1),lonlat(2),schmidt,ds,'albnir'
    Write(1,formout) albnirdata(:,:,i-1)*100.
    Close(1)
  End Do
else
  Open(1,File=fname(4))
  Write(1,'(i3,i4,2f8.3,f6.3,f8.0," ",a39)') sibdim(1),sibdim(2),lonlat(1),lonlat(2),schmidt,ds,'albnir'
  Write(1,formout) albnirdata(:,:,0)*100.
  Close(1)
end if
dimcount=(/ sibdim(1), sibdim(2), mthrng, 1 /)
Call ncwritedatgen(ncidarr,albnirdata,dimcount,varid(4))

! Write rsmin file
Write(6,*) 'Write rsmin files.'
Call calrsmin(landdata,laidata,sibdim,mthrng,rdata)
Where(rdata.GT.995.)
  rdata=995.
End Where
Write (formout,'("(",i3,"f5.0)" )') sibdim(1)
if (mthrng.eq.12) then
  Do i=1,12
    Write(monthout,'(I2.2)') i
    filedesc=trim(fname(5))//'.'//trim(monthout)
    Open(1,File=filedesc)
    Write(1,'(i3,i4,2f8.3,f6.3,f8.0," ",a39)') sibdim(1),sibdim(2),lonlat(1),lonlat(2),schmidt,ds,'rsmin'
    Write(1,formout) rdata(:,:,i)
    Close(1)
  End Do
else
  Open(1,File=fname(5))
  Write(1,'(i3,i4,2f8.3,f6.3,f8.0," ",a39)') sibdim(1),sibdim(2),lonlat(1),lonlat(2),schmidt,ds,'rsmin'
  Write(1,formout) rdata(:,:,1)
  Close(1)
end if
dimcount=(/ sibdim(1), sibdim(2), mthrng, 1 /)
Call ncwritedatgen(ncidarr,rdata,dimcount,varid(5))

! Write roughness file
Write(6,*) 'Write roughness files.'
Call calrough(landdata,laidata,sibdim,mthrng,rdata,zmin)
Where(rdata.LT.0.01)
  rdata=0.01
End Where
Write (formout,'("(",i3,"f6.0)" )') sibdim(1)
if (mthrng.eq.12) then
  Do i=1,12
    Write(monthout,'(I2.2)') i
    filedesc=trim(fname(6))//'.'//trim(monthout)
    Open(1,File=filedesc)
    Write(1,'(i3,i4,2f8.3,f6.3,f8.0," ",a39)') sibdim(1),sibdim(2),lonlat(1),lonlat(2),schmidt,ds,'rough'
    Write(1,formout) rdata(:,:,i)*100.
    Close(1)
  End Do  
else
    Open(1,File=fname(6))
    Write(1,'(i3,i4,2f8.3,f6.3,f8.0," ",a39)') sibdim(1),sibdim(2),lonlat(1),lonlat(2),schmidt,ds,'rough'
    Write(1,formout) rdata(:,:,1)*100.
    Close(1)
end if
dimcount=(/ sibdim(1), sibdim(2), mthrng, 1 /)
Call ncwritedatgen(ncidarr,rdata,dimcount,varid(6))

! Write veg frac file
!Write(6,*) 'Write vegetation fraction files.'
!Call calgreen(laidata,sibdim,mthrng,rdata)
!Write (formout,'("(",i3,"f5.0)" )') sibdim(1)
!if (mthrng.eq.12) then
!  Do i=1,12
!    Write(monthout,'(I2.2)') i
!    filedesc=trim(fname(7))//'.'//trim(monthout)
!    Open(1,File=filedesc)
!    Write(1,'(i3,i4,2f8.3,f6.3,f8.0," ",a39)') sibdim(1),sibdim(2),lonlat(1),lonlat(2),schmidt,ds,'vegfrac'
!    Write(1,formout) rdata(:,:,i)*100.
!    Close(1)
!  End Do
!else
!  Open(1,File=fname(7))
!  Write(1,'(i3,i4,2f8.3,f6.3,f8.0," ",a39)') sibdim(1),sibdim(2),lonlat(1),lonlat(2),schmidt,ds,'vegfrac'
!  Write(1,formout) rdata(:,:,1)*100.
!  Close(1)
!end if
!dimcount=(/ sibdim(1), sibdim(2), mthrng, 1 /)
!Call ncwritedatgen(ncidarr,rdata,dimcount,varid(7))


! Write lai file
Write(6,*) 'Write lai files.'
Write (formout,'("(",i3,"f5.0)" )') sibdim(1)
if (mthrng.eq.12) then
  Do i=1,12
    Write(monthout,'(I2.2)') i
    filedesc=trim(fname(7))//'.'//trim(monthout)
    Open(1,File=filedesc)
    Write(1,'(i3,i4,2f8.3,f6.3,f8.0," ",a39)') sibdim(1),sibdim(2),lonlat(1),lonlat(2),schmidt,ds,'lai'
    Write(1,formout) laidata(:,:,i-1)*100.
    Close(1)
  End Do
else
  Open(1,File=fname(7))
  Write(1,'(i3,i4,2f8.3,f6.3,f8.0," ",a39)') sibdim(1),sibdim(2),lonlat(1),lonlat(2),schmidt,ds,'lai'
  Write(1,formout) laidata(:,:,0)*100.
  Close(1)
end if
dimcount=(/ sibdim(1), sibdim(2), mthrng, 1 /)
Call ncwritedatgen(ncidarr,laidata,dimcount,varid(7))

! SiB type
Write(6,*) 'Write land-use type'
Do i=1,sibdim(1)
  Do j=1,sibdim(2)
    if (lsdata(i,j).ge.0.5) then
      idata(i,j)=0
    else
      sibmax=Maxloc(landdata(i,j,1:13))
      idata(i,j)=sibmax(1)
    end if
  End do
End do
Write(formout,'(1h(,i3,2hi3,1h))') sibdim(1)
Open(1,File=fname(8))
Write(1,'(i3,i4,2f8.3,f6.3,f8.0," ",a39)') sibdim(1),sibdim(2),lonlat(1),lonlat(2),schmidt,ds,'sib'
Write(1,formout) idata
Close(1)
dimcount=(/ sibdim(1), sibdim(2), 1, 1 /)
Call ncwritedatgen(ncidarr,Real(idata),dimcount,varid(8))

! Urban
Write(6,*) 'Write urban fraction'
urbanfrac=0.6
Write(formout,'("(",i3,"f4.0)" )') sibdim(1)
Open(1,File=fname(9))
Write(1,'(i3,i4,2f8.3,f6.3,f8.0," ",a39)') sibdim(1),sibdim(2),lonlat(1),lonlat(2),schmidt,ds,'urban'
Write(1,formout) urbandata(:,:)*urbanfrac*100.
Close(1)
dimcount=(/ sibdim(1), sibdim(2), 1, 1 /)
Call ncwritedatgen(ncidarr,urbandata*urbanfrac,dimcount,varid(9))

Call ncclose(ncidarr)

Deallocate(landdata,soildata,albvisdata,albnirdata,idata,rdata,urbandata,laidata,lsdata)

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fix SiB data
!

subroutine sibfix(landdata,rawdata,rlld,sibdim)

implicit none

integer, dimension(1:2), intent(in) :: sibdim
real, dimension(1:sibdim(1),1:sibdim(2),0:50), intent(in) :: rawdata
real, dimension(sibdim(1),sibdim(2),1:2), intent(in) :: rlld
real, dimension(1:sibdim(1),1:sibdim(2),0:13), intent(out) :: landdata
real, dimension(1:sibdim(1),1:sibdim(2),0:22) :: landtemp
logical, dimension(1:sibdim(1),1:sibdim(2)) :: sermsk
integer i,ilon,ilat,pxy(2)
real nsum,wsum

landtemp(:,:,0:22)=rawdata(:,:,0:22)
do i=31,50
  landtemp(:,:,i-30)=landtemp(:,:,i-30)+rawdata(:,:,i)
end do
landdata(:,:,1:12)=landtemp(:,:,1:12)
landdata(:,:,0)=landtemp(:,:,19)+landtemp(:,:,22) ! redefine water
landdata(:,:,13)=landtemp(:,:,20) ! redefine ice
landtemp(:,:,0:13)=landdata(:,:,0:13)

sermsk=sum(landdata(:,:,1:13),3).gt.0.
if (.not.any(sermsk)) return

do ilon=1,sibdim(1)
  do ilat=1,sibdim(2)
    wsum=landdata(ilon,ilat,0) ! water
    if (wsum.lt.1.) then
      nsum=sum(landdata(ilon,ilat,1:13)) ! land
      if (nsum.eq.0.) then
        call findnear(pxy,ilon,ilat,sermsk,rlld,sibdim)
        landdata(ilon,ilat,1:13)=landtemp(pxy(1),pxy(2),1:13)  
        nsum=sum(landdata(ilon,ilat,1:13))
      end if
      landdata(ilon,ilat,1:13)=landdata(ilon,ilat,1:13)*(1.-wsum)/nsum
    end if
  end do
end do

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! clean land and soil data
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
  sermsk=sermsk.and.(datain(:,:,ilon).gt.0.)
end do
if (.not.any(sermsk)) return

do ilon=1,sibdim(1)
  do ilat=1,sibdim(2)
    if ((1-nint(lsdata(ilon,ilat)).eq.1).and.(.not.sermsk(ilon,ilat))) then
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

Implicit None

Integer, intent(in) :: topounit
Integer, dimension(1:2), intent(in) :: sibdim
Integer ilout,ierr,ia,ib
Character(len=*), intent(in) :: toponame,topoout
Character*80 formout
Character*47 dc
Real, dimension(1:sibdim(1),1:sibdim(2)), intent(in) :: lsmskin,oceanin
Real, dimension(1:sibdim(1),1:sibdim(2)) :: topo,sd,lsmsk
Real ra,rb,rc,rd
ilout=Min(sibdim(1),30) ! To be compatiable with terread

Write(6,*) "Adjust topography data for consistancy with land-sea mask"

Open(topounit,FILE=toponame,FORM='formatted',STATUS='old',IOSTAT=ierr)
Read(topounit,*,IOSTAT=ierr) ia,ib,ra,rb,rc,rd,dc
Read(topounit,*,IOSTAT=ierr) topo ! Topography data
Read(topounit,*,IOSTAT=ierr) lsmsk ! land/sea mask (to be replaced)
Read(topounit,*,IOSTAT=ierr) sd ! Topography standard deviation
Close(topounit)

If (ierr.NE.0) then
  Write(6,*) "ERROR: Cannot read file ",trim(toponame)
  Stop
End if

lsmsk=Real(1-nint(lsmskin))
where (((1-nint(oceanin)).eq.0).and.((1-nint(lsmskin)).eq.0))
  topo(:,:)=0.
  sd(:,:)=0.
end where

Open(topounit,FILE=topoout,FORM='formatted',STATUS='replace',IOSTAT=ierr)
Write(topounit,'(i3,i4,2f8.3,f6.3,f8.0," ",a39)',IOSTAT=ierr) ia,ib,ra,rb,rc,rd,dc
Write(formout,'("(",i3,"f7.0)")',IOSTAT=ierr) ilout
Write(topounit,formout,IOSTAT=ierr) topo ! Topography data
Write(formout,'("(",i3,"f4.1)")',IOSTAT=ierr) ilout
Write(topounit,formout,IOSTAT=ierr) lsmsk ! land/sea mask
Write(formout,'("(",i3,"f6.0)")',IOSTAT=ierr) ilout
Write(topounit,formout,IOSTAT=ierr) sd ! Topography standard deviation
Close(topounit)

If (ierr.NE.0) then
  Write(6,*) "ERROR: Cannot write file ",trim(toponame)
  Stop
End if

Return
End
