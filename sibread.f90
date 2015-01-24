! This subroutine is to extract (in memory) data from the SiB dataset.
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads the sib, soil, albedo or urban data and maps
! it to a CC grid.  Data is specified by datatype (e.g., datatype=land
! or datatype=soil).
!

Subroutine getdata(dataout,glonlat,grid,tlld,sibdim,num,sibsize,datatype,fastsib,ozlaipatch,binlimit,month)

Use ccinterp

Implicit None

Integer, intent(in) :: sibsize,num,binlimit,month
Integer, dimension(1:2), intent(in) :: sibdim
Integer, dimension(1:sibdim(1),1:sibdim(2)) :: countn
Integer, dimension(1:2) :: lldim,lldim_x,llstore,pxy
Integer nscale,nscale_x,nface,subsec,mode,tmp
Integer i,j,k,lci,lcj,nx,ny,ni,nj,ci,cj,in,ix,jn,jx,is
Integer basesize,scalelimit,minscale,imth
Character(len=*), intent(in) :: datatype
Real, dimension(1:sibdim(1),1:sibdim(2),0:num), intent(out) :: dataout
Real, dimension(1:sibdim(1),1:sibdim(2)), intent(in) :: grid
Real, dimension(1:sibdim(1),1:sibdim(2),1:2), intent(in) :: tlld
Real, dimension(1:sibdim(1),1:sibdim(2),1:2) :: rlld
Real, dimension(1:sibdim(1),1:sibdim(2)) :: zsum
Real, dimension(1:2), intent(in) :: glonlat
Real, dimension(:,:,:), allocatable :: coverout
Real, dimension(1:2) :: latlon
Real, dimension(1:2,1:2) :: sll
Real, dimension(1:2,1:2,0:num) :: covertemp
Real aglon,aglat,alci,alcj,serlon,serlat,slonn,slatx,elon,elat,tscale,baselon
Real ipol,callon,callat,indexlon,indexlat
Logical, intent(in) :: fastsib,ozlaipatch
Logical, dimension(:,:), allocatable :: sermask

dataout=0.
countn=0

! Determine scale limits
nscale=999

baselon=real(int(glonlat(1)-180.))
rlld=tlld
Do While (Any(rlld(:,:,1).LT.baselon))
  Where (rlld(:,:,1).LT.baselon)
    rlld(:,:,1)=rlld(:,:,1)+360.
  End where
End do
Do While (Any(rlld(:,:,1).GT.(baselon+360.)))
  Where (rlld(:,:,1).GT.(baselon+360.))
    rlld(:,:,1)=rlld(:,:,1)-360.
  End where
End do

Select Case(datatype)
  Case('land')
    Write(6,*) 'Process USGS land-use dataset.'
    scalelimit=1
  Case('soil')
    Write(6,*) 'Process HWSD soil dataset.'
    scalelimit=4
  Case('lai')
    Write(6,*) 'Process mod15_BU LAI dataset.'
    scalelimit=4
  Case('albvis')
    Write(6,*) 'Process mod15_BU albedo (VIS) dataset.'
    scalelimit=6
  Case('albnir')
    Write(6,*) 'Process mod15_BU albedo (NIR) dataset.'
    scalelimit=6          
  Case DEFAULT
    Write(6,*) 'ERROR: Cannot find data ',trim(datatype)
    Stop
End Select

If (fastsib) then

  ! Step over scales
  mode=0
  Do While (Any(countn.EQ.0).AND.(nscale.GT.scalelimit))

    latlon=(/ baselon, 90. /)
    Call findsmallscale(nscale,scalelimit,latlon,llstore,grid,(countn.EQ.0),rlld,subsec,sll,sibsize,sibdim)

    slonn=sll(1,1)
    slatx=sll(2,2)

    minscale=nscale*binlimit

    Write(6,*) 'Bin'
    Write(6,*) 'nscale       = ',nscale
    Write(6,*) 'subsec       = ',subsec
    Write(6,*) 'sll          = ',sll
    Write(6,*) 'llstore      = ',llstore

    If (subsec.NE.0) then

      Do nx=1,subsec
        Do ny=1,subsec

          Write(6,*) 'nx,ny,subsec = ',nx,ny,subsec
      
          lldim=llstore
          ! Determine top corner lat/lon
          Call latlonconvert(nscale,latlon,lldim,slonn,slatx,nx,ny)
      
          Write(6,*) 'orig latlon  = ',latlon
          Write(6,*) 'orig lldim   = ',lldim

          ! Check if there are any points of interest on this tile
          Call searchdim(mode,sll,nscale,real(nscale),latlon,lldim,grid,(countn.EQ.0),rlld,sibdim)
          Call scaleconvert(nscale,tmp,lldim,sll,sibsize)
          mode=2
      
          latlon(1)=sll(1,1)
          latlon(2)=sll(2,2)

          Write(6,*) 'mod latlon   = ',latlon
          Write(6,*) 'mod lldim    = ',lldim

          ! Bin
          If (All(lldim.GT.0)) then

            Allocate(coverout(lldim(1),lldim(2),0:num))
	  
            Select Case(datatype)
	          Case('land')
                Call sibread(latlon,nscale,lldim,coverout)
              Case('soil')
	            Call kmconvert(nscale,nscale_x,lldim,lldim_x,4)
                Call soilread(latlon,nscale_x,lldim_x,coverout)
              Case('lai')
	            Call kmconvert(nscale,nscale_x,lldim,lldim_x,4)
	            if (num.eq.11) then
	              do imth=1,12
                    Call lairead(latlon,nscale_x,lldim_x,imth,coverout(:,:,imth-1))
                  end do
                else
                  Call lairead(latlon,nscale_x,lldim_x,month,coverout(:,:,0))
                end if
              Case('albvis','albnir')
	            Call kmconvert(nscale,nscale_x,lldim,lldim_x,6)
	            if (num.eq.11) then
	              do imth=1,12
                    Call albedoread(latlon,nscale_x,lldim_x,imth,coverout(:,:,imth-1),datatype)
                  end do
                else
                  Call albedoread(latlon,nscale_x,lldim_x,month,coverout(:,:,0),datatype)
                end if
              Case DEFAULT
                Write(6,*) 'ERROR: Cannot find data ',trim(datatype)
                Stop
            End Select

            Write(6,*) 'Start bin'
            Do i=1,lldim(1)
              Do j=1,lldim(2)
                aglon=callon(latlon(1),i,nscale)
                aglat=callat(latlon(2),j,nscale)
                Call lltoijmod(aglon,aglat,alci,alcj,nface)
                lci = nint(alci)
                lcj = nint(alcj)
                lcj = lcj+nface*sibdim(1)
                If (grid(lci,lcj).GE.real(minscale)) then
                  If (sum(coverout(i,j,:)).eq.0.) then
	                If (countn(lci,lcj).EQ.0) Then
                      dataout(lci,lcj,:)=-1. ! Missing value?
                      countn(lci,lcj)=1
                    End if
                  Else
                    If (dataout(lci,lcj,0).LT.0.) Then
                      dataout(lci,lcj,:)=0. ! reset missing point after finding non-trival data
                      countn(lci,lcj)=0
                    End If
                    dataout(lci,lcj,:)=dataout(lci,lcj,:)+coverout(i,j,:)
                    countn(lci,lcj)=countn(lci,lcj)+1
                  End If
                End if
              End Do
            End Do
            Write(6,*) 'Bin complete'

            Deallocate(coverout)

          Else
            Write(6,*) 'No points in valid range'
          End If
      
        End Do
      End Do

    Else
      Write(6,*) 'Skip'
    End If
  
  End Do

Else

  Select Case(datatype)
    Case('land')
      Call sibstream(sibdim,dataout,countn)
    Case('soil')
      Call soilstream(sibdim,dataout,countn)
    Case('lai')
      if (num.eq.11) then
        do imth=1,12
          Call laistream(sibdim,imth,dataout(:,:,imth-1),countn)
        end do
      else
        Call laistream(sibdim,month,dataout(:,:,0),countn)
      end if
    Case('albvis','albnir')
      if (num.eq.11) then
        do imth=1,12
          Call albedostream(sibdim,imth,dataout(:,:,imth-1),countn,datatype)
        end do
      else
        Call albedostream(sibdim,month,dataout(:,:,0),countn,datatype)
      end if
    Case DEFAULT
      Write(6,*) 'ERROR: Cannot find data ',trim(datatype)
      Stop
  End Select

End If

! Interpolate
Write(6,*) 'Interpolate'
Allocate(sermask(1:2,1:2))
nscale=scalelimit

latlon=(/ baselon, 90. /)
llstore=(/ 43200/nscale , 21600/nscale /)
Call searchdim(4,sll,nscale,0.,latlon,llstore,grid,(countn.EQ.0),rlld,sibdim)
Call scaleconvert(nscale,subsec,llstore,sll,sibsize)
slonn=sll(1,1)
slatx=sll(2,2)

Write(6,*) 'nscale       = ',nscale
Write(6,*) 'subsec       = ',subsec
Write(6,*) 'sll          = ',sll
Write(6,*) 'llstore      = ',llstore

If (subsec.NE.0) then
  Do nx=1,subsec
    Do ny=1,subsec

      Write(6,*) 'nx,ny,subsec = ',nx,ny,subsec

      lldim=llstore
      Call latlonconvert(nscale,latlon,lldim,slonn,slatx,nx,ny)

      Write(6,*) 'orig latlon  = ',latlon
      Write(6,*) 'orig lldim   = ',lldim

      ! overlap tiles for interpolation
      If (nx.NE.subsec) lldim(1)=lldim(1)+1
      If (ny.NE.subsec) lldim(2)=lldim(2)+1
    
      ! Check if there are any points of interest on this tile
      Call searchdim(4,sll,nscale,0.,latlon,lldim,grid,(countn.EQ.0),rlld,sibdim)
      Call scaleconvert(nscale,tmp,lldim,sll,sibsize)
      If (Any(lldim(:).EQ.1)) lldim=0
      
      latlon(1)=sll(1,1)
      latlon(2)=sll(2,2)

      Write(6,*) 'mod latlon   = ',latlon
      Write(6,*) 'mod lldim    = ',lldim

      If ((lldim(1).GT.0).AND.(lldim(2).GT.0)) then

        Allocate(coverout(lldim(1),lldim(2),0:num))
	
        Select Case(datatype)
          Case('land')
            Call sibread(latlon,nscale,lldim,coverout)
          Case('soil')
            Call kmconvert(nscale,nscale_x,lldim,lldim_x,4)
            Call soilread(latlon,nscale_x,lldim_x,coverout)
          Case('lai')
            Call kmconvert(nscale,nscale_x,lldim,lldim_x,4)
            if (num.eq.11) then
              do imth=1,12
                Call lairead(latlon,nscale_x,lldim_x,imth,coverout(:,:,imth-1))
              end do
            else
              Call lairead(latlon,nscale_x,lldim_x,month,coverout(:,:,0))
            end if
          Case('albvis','albnir')
            Call kmconvert(nscale,nscale_x,lldim,lldim_x,6)
            if (num.eq.11) then
              do imth=1,12
                Call albedoread(latlon,nscale_x,lldim_x,imth,coverout(:,:,imth-1),datatype)
              end do
            else
              Call albedoread(latlon,nscale_x,lldim_x,month,coverout(:,:,0),datatype)
            end if
         Case DEFAULT
            Write(6,*) 'ERROR: Cannot find data ',trim(datatype)
            Stop
        End Select

        Do lci=1,sibdim(1)
          Do lcj=1,sibdim(2)
            If (countn(lci,lcj).EQ.0) then
  
              aglon=rlld(lci,lcj,1)
              aglat=rlld(lci,lcj,2)
	      
              serlon=indexlon(aglon,latlon(1),nscale)
              serlat=indexlat(aglat,latlon(2),nscale)

              i=int(serlon)
              j=int(serlat)

              If ((i.GE.1).AND.(i.LT.lldim(1)).AND.(j.GE.1).AND.(j.LT.lldim(2))) Then
	      
                ni=i+1 ! 4 point interpolation
                nj=j+1
	      
                serlon=serlon-real(i)
                serlat=serlat-real(j)
	
                If (datatype.eq.'land') then		
                  covertemp(1:2,1:2,:)=coverout(i:ni,j:nj,:)
                else
                  Call realfill(covertemp,coverout,lldim,i,ni,j,nj,num)
                End if	    
                Do k=0,num
                  dataout(lci,lcj,k)=ipol(covertemp(:,:,k),serlon,serlat)
                End Do
                countn(lci,lcj)=1
              End If
                
            End If
          End Do
        End Do
        Deallocate(coverout)

      Else
        Write(6,*) 'No points in valid range'
      End If
    End Do
  End Do

Else
  Write(6,*) 'Skip'
End If

Deallocate(sermask)


If (Any(countn.LT.1)) then
  ! This will happen near the poles
  Write(6,*) "Interpolation failed - Using near nbr"
  Allocate(sermask(1:sibdim(1),1:sibdim(2)))
  sermask(:,:)=countn(:,:).GT.0
  If (Any(sermask)) then
    Do lci=1,sibdim(1)
      Do lcj=1,sibdim(2)
        If (countn(lci,lcj).EQ.0) then
          call findnear(pxy,lci,lcj,sermask,rlld,sibdim)
          dataout(lci,lcj,:)=dataout(pxy(1),pxy(2),:)
     	  countn(lci,lcj)=1
        End if
      End do
    End Do
  Else
    Write(6,*) 'WARN: Cannot find any non-trivial near nbrs'
    Write(6,*) '      Assume data is trivial'
    dataout=0.
    countn=1
  End if	
  Deallocate(sermask)
End If

where(dataout.lt.0.)
  dataout=0.
end where

Do k=0,num
  dataout(:,:,k)=dataout(:,:,k)/Real(countn)
End Do

if (ozlaipatch.and.(datatype.eq.'lai')) call ozlai(sibdim,num,dataout,rlld,month)

Write(6,*) "Task complete"

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads sib data down to nscale=1km resolution
!

Subroutine sibread(latlon,nscale,lldim,coverout)

Implicit None

Integer, intent(in) :: nscale
Real, dimension(1:2), intent(in) :: latlon
Integer, dimension(1:2), intent(in) :: lldim
Real, dimension(lldim(1),lldim(2),0:81), intent(out) :: coverout
Integer*1, dimension(1:43200,1:nscale) :: databuffer
Integer*1, dimension(1:43200) :: datatemp
Integer, dimension(1:2,1:2) :: jin,jout
Integer ilat,ilon,jlat,recpos
Integer, dimension(1:2) :: llint
integer, dimension(860,700) :: idean
integer tix,tiy,jj,ii
integer ni1,nj1,ni2,nj2,ni3,nj3,ni4,nj4

open(11,file='dean31_int',form='formatted',status='old')
read(11,*) ni1,nj1
read(11,*) ni2,nj2
read(11,*) ni3,nj3
read(11,*) ni4,nj4
do jj=1,700
  read(11,*) (idean(ii,jj),ii=1,860)
end do
close(11)

! Must be compiled using 1 byte record lengths
Open(10,FILE='gsib2u.img',ACCESS='DIRECT',FORM='UNFORMATTED',RECL=43200)

! To speed-up the code, 43200x(nscale) blocks of the sib file are read one row
! at a time.  The data is then averaged in memory.  This system speeds-up the
! code considerably.  However, there are limitations as to the size of data
! that can be stored in memory.

Call solvejshift(latlon(1),jin,jout,120)

Do ilat=1,lldim(2)

  if ((mod(ilat,10)==0).or.(ilat==lldim(2))) then
    Write(6,*) 'USGS - ',ilat,'/',lldim(2)
  end if
  
  ! Read data
  llint(2)=nint((90.-latlon(2))*120.)+(ilat-1)*nscale
  Do jlat=1,nscale
    recpos=llint(2)+jlat
    Read(10,REC=recpos) datatemp
    ! Dean G's data
    tiy=nint((45.+90.-(real(llint(2)+jlat)-0.5)/120.)/0.05-0.5)
    if (tiy>=1.and.tiy<=700) then
      do ilon=1,43200
        tix=nint(((real(ilon)-0.5)/120.-180.-112.)/0.05+0.5)
        if (tix>=1.and.tix<=860) then
          datatemp(ilon)=idean(tix,tiy)+50
          if (datatemp(ilon)==50) datatemp(ilon)=19 ! water
        end if
      end do
    end if
    
    ! Shift lon to zero
    databuffer(jin(1,1):jin(1,2),jlat)=datatemp(jout(1,1):jout(1,2))
    databuffer(jin(2,1):jin(2,2),jlat)=datatemp(jout(2,1):jout(2,2))
  End Do
  
  Do ilon=1,lldim(1)
    llint(1)=(ilon-1)*nscale
    Call dataconvert(databuffer(llint(1)+1:llint(1)+nscale,1:nscale),coverout(ilon,ilat,:),nscale,81)
  End Do
End Do

Close(10)

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads soil data down to nscale_4=1
! (i.e., 4km) resolution.
!

Subroutine soilread(latlon,nscale_4,lldim_4,coverout)

Implicit None

Integer, intent(in) :: nscale_4
Real, dimension(1:2), intent(in) :: latlon
Integer, dimension(1:2), intent(in) :: lldim_4
Real, dimension(lldim_4(1),lldim_4(2),0:8), intent(out) :: coverout
real, dimension(0:13) :: faosoil
Integer*1, dimension(1:10800,1:nscale_4) :: databuffer
Integer*1, dimension(1:10800) :: datatemp
Integer, dimension(1:2,1:2) :: jin,jout
Integer ilat,ilon,jlat,recpos,i
Integer, dimension(1:2) :: llint_4
real nsum
integer, dimension(0:13), parameter :: masmap=(/ 0, 1, 1, 4, 2, 4, 7, 2, 2, 5, 6, 3, 8, 9 /)

! Must be compiled using 1 byte record lengths
Open(20,FILE='usda4.img',ACCESS='DIRECT',FORM='UNFORMATTED',RECL=10800)

! To speed-up the code, 43200x(nscale) blocks of the sib file are read
! at a time.  The data is then averaged in memory.  This system speeds-up the
! code considerably.  However, there are limitations as to the size of data
! that can be stored in memory.

Call solvejshift(latlon(1),jin,jout,30)

Do ilat=1,lldim_4(2)

  if ((mod(ilat,10).eq.0).or.(ilat.eq.lldim_4(2))) then
    Write(6,*) 'HWSD - ',ilat,'/',lldim_4(2)
  end if
  
  ! Read data
  llint_4(2)=nint((90.-latlon(2))*30.)+(ilat-1)*nscale_4
  Do jlat=1,nscale_4
    recpos=llint_4(2)+jlat
    Read(20,REC=recpos) datatemp
    ! Shift lon to zero
    databuffer(jin(1,1):jin(1,2),jlat)=datatemp(jout(1,1):jout(1,2))
    databuffer(jin(2,1):jin(2,2),jlat)=datatemp(jout(2,1):jout(2,2))
  End Do
  
  Do ilon=1,lldim_4(1)
    llint_4(1)=(ilon-1)*nscale_4
    Call dataconvert(databuffer(llint_4(1)+1:llint_4(1)+nscale_4,1:nscale_4),faosoil,nscale_4,13)
    nsum=sum(faosoil(1:13))
    if (nsum.gt.0.) then
      coverout(ilon,ilat,:)=0.
      do i=1,13
        coverout(ilon,ilat,masmap(i)-1)=coverout(ilon,ilat,masmap(i)-1)+faosoil(i)/nsum
      end do
    else
      coverout(ilon,ilat,:)=0.
    end if
  End Do
End Do

Close(20)

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads LAI data down to nscale_4=1
! (i.e., 4km) resolution.
!

Subroutine lairead(latlon,nscale_4,lldim_4,imth,dataout)

Implicit None

Integer, intent(in) :: nscale_4,imth
Real, dimension(1:2), intent(in) :: latlon
Integer, dimension(1:2), intent(in) :: lldim_4
Real, dimension(lldim_4(1),lldim_4(2)), intent(out) :: dataout
Integer, dimension(1:10800,1:nscale_4) :: databuffer
Integer*1, dimension(1:10800) :: datatemp
Integer, dimension(1:2) :: llint_4
Integer ilat,ilon,jlat,recpos
Integer, dimension(1:2,1:2) :: jin,jout
Character*10 fname
Logical, dimension(1:nscale_4,1:nscale_4) :: sermask

Call solvejshift(latlon(1),jin,jout,30)

write(fname,'("slai",I2.2,".img")') imth
write(6,*) 'Reading ',trim(fname)

! Must be compiled using 1 byte record lengths
Open(30,FILE=fname,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=10800)

! To speed-up the code, 43200x(nscale) blocks of the sib file are read
! at a time.  The data is then averaged in memory.  This system speeds-up the
! code considerably.  However, there are limitations as to the size of data
! that can be stored in memory.

Do ilat=1,lldim_4(2)

  if ((mod(ilat,10).eq.0).or.(ilat.eq.lldim_4(2))) then
    Write(6,*) 'LAI - ',ilat,'/',lldim_4(2)
  end if
  
  ! Read data
  llint_4(2)=nint((90.-latlon(2))*30.)+(ilat-1)*nscale_4
  Do jlat=1,nscale_4
    recpos=llint_4(2)+jlat
    Read(30,REC=recpos) datatemp
    ! Shift lon to zero
    databuffer(jin(1,1):jin(1,2),jlat)=datatemp(jout(1,1):jout(1,2))
    databuffer(jin(2,1):jin(2,2),jlat)=datatemp(jout(2,1):jout(2,2))
  End Do
  
  where (databuffer.eq.0)
    databuffer=1 ! min LAI
  end where
  where (databuffer.lt.0)
    databuffer=databuffer+256
  end where
  
  Do ilon=1,lldim_4(1)
    llint_4(1)=(ilon-1)*nscale_4
    sermask=(databuffer(llint_4(1)+1:llint_4(1)+nscale_4,1:nscale_4).gt.0)
    sermask=sermask.and.(databuffer(llint_4(1)+1:llint_4(1)+nscale_4,1:nscale_4).lt.100)
    if (Any(sermask)) then
      dataout(ilon,ilat)=real(sum(databuffer(llint_4(1)+1:llint_4(1)+nscale_4,1:nscale_4),sermask)) &
                        /(real(count(sermask))*10.)
    else
      dataout(ilon,ilat)=0. ! missing value flag
    end if
  End Do
End Do

Close(30)

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads albedo data down to nscale_6=1
! (i.e., 6km) resolution.
!

Subroutine albedoread(latlon,nscale_6,lldim_6,imth,dataout,datatype)

Implicit None

Integer, intent(in) :: nscale_6,imth
Real, dimension(1:2), intent(in) :: latlon
Integer, dimension(1:2), intent(in) :: lldim_6
Real, dimension(lldim_6(1),lldim_6(2)), intent(out) :: dataout
Character(len=*), intent(in) :: datatype
Integer, dimension(1:7200,1:nscale_6) :: databuffer
Integer*1, dimension(1:7200) :: datatemp
Integer, dimension(1:2) :: llint_6
Integer ilat,ilon,jlat,recpos
Integer, dimension(1:2,1:2) :: jin,jout
Character*12 fname
Logical, dimension(1:nscale_6,1:nscale_6) :: sermask

Call solvejshift(latlon(1),jin,jout,20)

select case(datatype)
  case ('albvis')
    write(fname,'("albvis",I2.2,".img")') imth
  case ('albnir')
    write(fname,'("albnir",I2.2,".img")') imth
end select
write(6,*) 'Reading ',trim(fname)

! Must be compiled using 1 byte record lengths
Open(40,FILE=fname,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=7200)

! To speed-up the code, 43200x(nscale) blocks of the sib file are read
! at a time.  The data is then averaged in memory.  This system speeds-up the
! code considerably.  However, there are limitations as to the size of data
! that can be stored in memory.

Do ilat=1,lldim_6(2)

  if ((mod(ilat,10).eq.0).or.(ilat.eq.lldim_6(2))) then
    Write(6,*) 'Albedo - ',ilat,'/',lldim_6(2)
  end if
  
  ! Read data
  llint_6(2)=nint((90.-latlon(2))*20.)+(ilat-1)*nscale_6
  Do jlat=1,nscale_6
    recpos=llint_6(2)+jlat
    Read(40,REC=recpos) datatemp
    ! Shift lon to zero
    databuffer(jin(1,1):jin(1,2),jlat)=datatemp(jout(1,1):jout(1,2))
    databuffer(jin(2,1):jin(2,2),jlat)=datatemp(jout(2,1):jout(2,2))
  End Do
  
  where (databuffer.lt.0)
    databuffer=databuffer+256
  end where
  
  Do ilon=1,lldim_6(1)
    llint_6(1)=(ilon-1)*nscale_6
    sermask=(databuffer(llint_6(1)+1:llint_6(1)+nscale_6,1:nscale_6).gt.0)
    sermask=sermask.and.(databuffer(llint_6(1)+1:llint_6(1)+nscale_6,1:nscale_6).le.100)
    if (Any(sermask)) then
      dataout(ilon,ilat)=real(sum(databuffer(llint_6(1)+1:llint_6(1)+nscale_6,1:nscale_6),sermask)) &
                        /(real(count(sermask))*100.)
    else
      dataout(ilon,ilat)=0. ! missing value flag
    end if
  End Do
End Do

Close(40)

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads sibsystems data at nscale=1km resolution
! (i.e., no storage, simply read and bin)
!

Subroutine sibstream(sibdim,coverout,countn)

Use ccinterp

Implicit None

Integer, dimension(2), intent(in) :: sibdim
Real, dimension(1:sibdim(1),1:sibdim(2),0:81), intent(out) :: coverout
Real aglon,aglat,alci,alcj
Real callon,callat
Integer, dimension(1:sibdim(1),1:sibdim(2)), intent(out) :: countn
integer, dimension(860,700) :: idean
Integer*1, dimension(1:43200) :: databuffer
Integer ilat,ilon,lci,lcj,nface,cpos
integer ii,jj,tix,tiy
integer ni1,nj1,ni2,nj2,ni3,nj3,ni4,nj4

coverout=0
countn=0

Write(6,*) "Read USGS data (stream)"

open(11,file='dean31_int',form='formatted',status='old')
read(11,*) ni1,nj1
read(11,*) ni2,nj2
read(11,*) ni3,nj3
read(11,*) ni4,nj4
do jj=1,700
  read(11,*) (idean(ii,jj),ii=1,860)
end do
close(11)

! Must be compiled using 1 byte rsibrd lengths
Open(10,FILE='gsib2u.img',ACCESS='DIRECT',FORM='UNFORMATTED',RECL=43200)

Do ilat=1,21600

  if (mod(ilat,10).eq.0) then
    Write(6,*) 'USGS - ',ilat,'/ 21600'
  end if
  
  ! Read data
  Read(10,REC=ilat) databuffer
  aglat=callat(90.,ilat,1)
  
  Do ilon=1,43200
    aglon=callon(-180.,ilon,1)
    
    if (.not.(aglat>-45..and.aglat<-10..and.aglon>112..and.aglon<154.5)) then
      Call lltoijmod(aglon,aglat,alci,alcj,nface)
      lci = nint(alci)
      lcj = nint(alcj)
      lcj = lcj+nface*sibdim(1)
    
      cpos=databuffer(ilon)
      coverout(lci,lcj,cpos)=coverout(lci,lcj,cpos)+1.
      countn(lci,lcj)=countn(lci,lcj)+1
    end if
    
  End Do
End Do

do tiy=1,700
  aglat=-44.975+real(tiy-1)*0.05
  do tix=1,850
    aglon=112.25+real(tix-1)*0.05

    Call lltoijmod(aglon,aglat,alci,alcj,nface)
    lci = nint(alci)
    lcj = nint(alcj)
    lcj = lcj+nface*sibdim(1)
    
    cpos=databuffer(ilon)+50
    if (cpos==50) cpos=19
    coverout(lci,lcj,cpos)=coverout(lci,lcj,cpos)+1.
    countn(lci,lcj)=countn(lci,lcj)+1
  end do
end do

Close(10)

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads soil at nscale_4=1
! (i.e., 4km) resolution.
! (i.e., no storage, simply read and bin)
!

Subroutine soilstream(sibdim,coverout,countn)

Use ccinterp

Implicit None

Integer, dimension(2), intent(in) :: sibdim
Real, dimension(1:sibdim(1),1:sibdim(2),0:8), intent(out) :: coverout
Real aglon,aglat,alci,alcj
Real callon,callat
Integer, dimension(1:sibdim(1),1:sibdim(2)), intent(out) :: countn
Integer*1, dimension(1:10800) :: databuffer
Integer ilat,ilon,lci,lcj,nface,cpos
integer, dimension(0:13), parameter :: masmap=(/ 0, 1, 1, 4, 2, 4, 7, 2, 2, 5, 6, 3, 8, 9 /)

coverout=0
countn=0

Write(6,*) "Read HWSD data (stream)"

! Must be compiled using 1 byte rsibrd lengths
Open(20,FILE='usda4.img',ACCESS='DIRECT',FORM='UNFORMATTED',RECL=10800)

Do ilat=1,5400

  if (mod(ilat,10).eq.0) then
    Write(6,*) 'HWSD - ',ilat,'/ 5400'
  end if
  
  ! Read data
  Read(20,REC=ilat) databuffer
  aglat=callat(90.,ilat,1)
  
  Do ilon=1,10800
    
    aglon=callon(-180.,ilon,1)
    
    Call lltoijmod(aglon,aglat,alci,alcj,nface)
    lci = nint(alci)
    lcj = nint(alcj)
    lcj = lcj+nface*sibdim(1)
    
    cpos=databuffer(ilon)
    if ((cpos.ge.1).and.(cpos.le.13)) then
      If (coverout(lci,lcj,0).LT.0.) then
        coverout(lci,lcj,:)=0.
        countn(lci,lcj)=0
      End If
      coverout(lci,lcj,masmap(cpos)-1)=coverout(lci,lcj,masmap(cpos)-1)+1.
      countn(lci,lcj)=countn(lci,lcj)+1
    else
      If (countn(lci,lcj).EQ.0) then
        coverout(lci,lcj,:)=-1.
        countn(lci,lcj)=1
      End if    
    end if
    
  End Do
End Do

Close(20)

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads LAI at nscale_4=1
! (i.e., 4km) resolution.
! (i.e., no storage, simply read and bin)
!

Subroutine laistream(sibdim,imth,coverout,countn)

Use ccinterp

Implicit None

integer, intent(in) :: imth
Integer, dimension(2), intent(in) :: sibdim
Real, dimension(1:sibdim(1),1:sibdim(2)), intent(out) :: coverout
Real aglon,aglat,alci,alcj
Real callon,callat
Integer, dimension(1:sibdim(1),1:sibdim(2)), intent(out) :: countn
Integer, dimension(1:10800) :: databuffer
Integer*1, dimension(1:10800) :: datatemp
Integer ilat,ilon,lci,lcj,nface
Character*10 fname

coverout=0
countn=0

write(fname,'("slai",I2.2,".img")') imth
write(6,*) 'Reading (stream) ',trim(fname)

! Must be compiled using 1 byte record lengths
Open(30,FILE=fname,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=10800)

Do ilat=1,5400
  aglat=callat(90.,ilat,1)

  if (mod(ilat,10).eq.0) then
    Write(6,*) 'LAI - ',ilat,'/ 5400'
  end if
  
  ! Read data
  Read(30,REC=ilat) datatemp
  
  databuffer=datatemp
  where (databuffer.eq.0)
    databuffer=1
  end where
  where (databuffer.lt.0)
    databuffer=databuffer+256
  end where
  
  Do ilon=1,10800
    aglon=callon(-180.,ilon,1)
    
    Call lltoijmod(aglon,aglat,alci,alcj,nface)
    lci = nint(alci)
    lcj = nint(alcj)
    lcj = lcj+nface*sibdim(1)
    
    If (databuffer(ilon).LT.100) then
      If (coverout(lci,lcj).LT.0.) then
        coverout(lci,lcj)=0.
        countn(lci,lcj)=0
      End If
      coverout(lci,lcj)=coverout(lci,lcj)+real(databuffer(ilon))/10.
      countn(lci,lcj)=countn(lci,lcj)+1
    Else
      If (countn(lci,lcj).EQ.0) then
        coverout(lci,lcj)=-1.
        countn(lci,lcj)=1
      End if
    End if    
  End Do
End Do

Close(30)

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads albedo at nscale_6=1
! (i.e., 6km) resolution.
! (i.e., no storage, simply read and bin)
!

Subroutine albedostream(sibdim,imth,coverout,countn,datatype)

Use ccinterp

Implicit None

integer, intent(in) :: imth
Integer, dimension(2), intent(in) :: sibdim
Real, dimension(1:sibdim(1),1:sibdim(2)), intent(out) :: coverout
character(len=*), intent(in) :: datatype
Real aglon,aglat,alci,alcj
Real callon,callat
Integer, dimension(1:sibdim(1),1:sibdim(2)), intent(out) :: countn
Integer, dimension(1:7200) :: databuffer
Integer*1, dimension(1:7200) :: datatemp
Integer ilat,ilon,lci,lcj,nface
Character*12 fname

coverout=0
countn=0

select case(datatype)
  case ('albvis')
    write(fname,'("albvis",I2.2,".img")') imth
  case ('albnir')
    write(fname,'("albnir",I2.2,".img")') imth
end select
write(6,*) 'Reading (stream) ',trim(fname)

! Must be compiled using 1 byte record lengths
Open(40,FILE=fname,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=7200)

Do ilat=1,3600
  aglat=callat(90.,ilat,1)

  if (mod(ilat,10).eq.0) then
    Write(6,*) 'Albedo - ',ilat,'/ 3600'
  end if
  
  ! Read data
  Read(40,REC=ilat) datatemp
  
  databuffer=datatemp
  where (databuffer.lt.0)
    databuffer=databuffer+256
  end where
  
  Do ilon=1,7200
    aglon=callon(-180.,ilon,1)
    
    Call lltoijmod(aglon,aglat,alci,alcj,nface)
    lci = nint(alci)
    lcj = nint(alcj)
    lcj = lcj+nface*sibdim(1)
    
    If (databuffer(ilon).LE.100) then
      If (coverout(lci,lcj).LT.0.) then
        coverout(lci,lcj)=0.
        countn(lci,lcj)=0
      End If
      coverout(lci,lcj)=coverout(lci,lcj)+real(databuffer(ilon))/100.
      countn(lci,lcj)=countn(lci,lcj)+1
    Else
      If (countn(lci,lcj).EQ.0) then
        coverout(lci,lcj)=-1.
        countn(lci,lcj)=1
      End if
    End if    
  End Do
End Do

Close(40)

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine aligns the data with the requested lat/lon
!

Subroutine solvejshift(lonin,jin,jout,nscale)

Implicit None

Real, intent(in) :: lonin
Integer, intent(in) :: nscale
Integer, dimension(1:2,1:2), intent(out) :: jin,jout
Integer jshift

jshift=Mod(nint(lonin*real(nscale))+180*nscale,360*nscale)
If (jshift.LT.0) jshift=jshift+360*nscale
jin(1,1)=1
jout(1,1)=Mod(jin(1,1)+jshift,360*nscale)
If (jout(1,1).EQ.0) jout(1,1)=360*nscale
jout(1,2)=360*nscale
jin(1,2)=jout(1,2)-jout(1,1)+jin(1,1)
jin(2,1)=jin(1,2)+1
jin(2,2)=360*nscale
jout(2,1)=1
jout(2,2)=jout(1,1)-1
If (jin(2,1).GT.jin(2,2)) then
  jin(2,:)=jin(1,:)
  jout(2,:)=jout(1,:)
End if

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine averages sib data
!

Subroutine dataconvert(datain,raw,nscale,num)

Implicit None

Integer, intent(in) :: nscale,num
Integer*1, dimension(1:nscale,1:nscale), intent(in) :: datain
Real, dimension(0:num), intent(out) :: raw
Integer i,j,datatmp

! Aggregate land use
! Faster to step over grid once
raw=0.
Do i=1,nscale
  Do j=1,nscale
    datatmp=datain(i,j)
    if (datatmp.LT.0) datatmp=datatmp+256
    raw(datatmp)=raw(datatmp)+1.
  End Do
End Do

raw=raw/real(nscale**2)

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine converts from 1km nscale to adj km nscale
!

Subroutine kmconvert(nscale,nscale_x,lldim,lldim_x,adj)

Implicit None

Integer, intent(in) :: nscale,adj
Integer, intent(out) :: nscale_x
Integer, dimension(1:2), intent(in) :: lldim
Integer, dimension(1:2), intent(out) :: lldim_x
Integer i

nscale_x=Int(nscale/adj)
If (nscale_x.LT.1) nscale_x=1

Do i=1,2
  lldim_x(i)=Int(Real(lldim(i))*Real(nscale)/(Real(nscale_x)*real(adj)))
  If (lldim_x(i).LT.1) lldim_x(i)=1
End Do

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine replaces zero values with the average of the
! surrounding points (ignoring ocean points).  This helps the interpolation
! scheme
!

Subroutine realfill(coverout,coverin,lldim,i,ni,j,nj,num)

Implicit None

Integer, intent(in) :: i,ni,j,nj,num
Integer, dimension(2), intent(in) :: lldim
Integer ti,tj,tk,ia,ib,ja,jb
Real, dimension(1:2,1:2,0:num), intent(out) :: coverout
Real, dimension(1:lldim(1),1:lldim(2),0:num), intent(in) :: coverin
Logical sermask(3,3)

coverout=0.

If (Any(coverin(i:ni,j:nj,:).gt.0.)) then
  Do ti=i,ni
    Do tj=j,nj
      If (Any(coverin(ti,tj,:).le.0.)) then
        ia=Max(ti-1,1)
        ib=Min(ti+1,lldim(1))
        ja=Max(tj-1,1)
        jb=Min(tj+1,lldim(2))
    
        sermask(1:ib-ia+1,1:jb-ja+1)=.TRUE.
        Do tk=0,num
          sermask(1:ib-ia+1,1:jb-ja+1)=sermask(1:ib-ia+1,1:jb-ja+1).AND.(coverin(ia:ib,ja:jb,tk).gt.0.)
        End Do
        Do tk=0,num
          coverout(ti-i+1,tj-j+1,tk)=Sum(coverin(ia:ib,ja:jb,tk),sermask(1:ib-ia+1,1:jb-ja+1))/Real(Count(sermask(1:ib-ia+1,1:jb-ja+1)))
        End do
      Else
        coverout(ti-i+1,tj-j+1,:)=coverin(ti,tj,:)
      End if
    End Do
  End Do
End if

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine determines how many subsections (blocks) the data
! needs to be broken into so as to fit into memory (i.e., sibsize)
!

Subroutine scaleconvert(nscale,subsec,lldim,sll,sibsize)

Implicit None

Integer, intent(out) :: subsec
Integer, intent(in) :: nscale,sibsize
Integer, dimension(1:2), intent(out) :: lldim
Real, dimension(1:2,1:2), intent(in) :: sll
Integer i,j

i=nint((sll(1,2)-sll(1,1))*120./Real(nscale))
j=nint((sll(2,2)-sll(2,1))*120./Real(nscale))

subsec=int(sqrt(real(i)*real(j)/(real(sibsize)**2)))+1
If (subsec.LT.1) subsec=1
  
lldim(1)=nint(real(i)/real(subsec))
lldim(2)=nint(real(j)/real(subsec))

If ((real(lldim(1)*nscale*subsec)).LT.((sll(1,2)-sll(1,1))*120.)) lldim(1)=lldim(1)+1
If ((real(lldim(2)*nscale*subsec)).LT.((sll(2,2)-sll(2,1))*120.)) lldim(2)=lldim(2)+1
If ((nint((90.-sll(2,2))*120.)+lldim(2)*nscale).GT.21600) lldim(2)=(21600-nint((90.-sll(2,2))*120.))/nscale
If ((lldim(1)*nscale).GT.43200) lldim(1)=43200/nscale

If ((lldim(1).LT.1).OR.(lldim(2).LT.1)) Then
  lldim=(/ 0, 0 /)
  subsec=0
End If

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine adjusts the latlon array for a specified subsection
!

Subroutine latlonconvert(nscale,latlon,lldim,slonn,slatx,nx,ny)

Implicit None

Integer, intent(in) :: nscale
Real, dimension(1:2), intent(out) :: latlon
Integer, dimension(1:2), intent(in) :: lldim
Real, intent(in) :: slonn,slatx
Integer, intent(in) :: nx,ny

latlon=(/ slonn+Real((nx-1)*lldim(1)*nscale)/120., slatx-Real((ny-1)*lldim(2)*nscale)/120. /)
if (latlon(2).LT.-90.) latlon(2)=-90.

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This function calculates lon from an array index
!

Real function callon(latlon,i,nscale)

Implicit None

Real, intent(in) :: latlon
Integer, intent(in) :: i,nscale

callon=(Real(i-1)+0.5)*real(nscale)/120.+latlon

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This function rounds up
!

Integer Function rndup(x)

Implicit None

Real, intent(in) :: x

rndup=int(x)
if (x.GT.real(rndup)) rndup=rndup+1

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This function calculates lat from an array index
!

Real function callat(latlon,i,nscale)

Implicit None

Real, intent(in) :: latlon
Integer, intent(in) :: i,nscale

callat=latlon-(Real(i-1)+0.5)*real(nscale)/120.

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This function calculates the array index for a specified lon
!

Real function indexlon(aglon,latlon,nscale)

Implicit None

Real, intent(in) :: aglon,latlon
Integer, intent(in) :: nscale

indexlon=(aglon-latlon)*120./real(nscale)+0.5
	    
Return
End	    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This function calculates the array index for a specified lat
!

Real function indexlat(aglat,latlon,nscale)

Implicit None

Real, intent(in) :: aglat,latlon
Integer, intent(in) :: nscale
   
indexlat=(-aglat+latlon)*120./real(nscale)+0.5

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine determines the boundary of a subsection for a
! specified scale.
!

Subroutine searchdim(mode,sll,nscale,scalelimit,latlon,lldim,grid,maskn,rlld,sibdim)

Implicit None

Integer, intent(in) :: mode,nscale
Integer, dimension(1:2), intent(in) :: lldim,sibdim
Real, intent(in) :: scalelimit
Real, dimension(1:2,1:2), intent(out) :: sll
Real, dimension(1:2), intent(in) :: latlon
Real, dimension(1:sibdim(1),1:sibdim(2)), intent(in) :: grid
Real, dimension(1:sibdim(1),1:sibdim(2),1:2), intent(in) :: rlld
Real, dimension(1:sibdim(1),1:sibdim(2),1:2) :: tlld
Real, dimension(1:2,1:2) :: templl
Integer, dimension(1:2,1:2,1:2) :: posll
Integer i,j
Integer rndup
Logical, dimension(1:sibdim(1),1:sibdim(2)) :: sermask
Integer, dimension(1:sibdim(1),1:sibdim(2)), intent(in) :: maskn

tlld=rlld

templl(1,1)=latlon(1)
templl(1,2)=latlon(1)+real(lldim(1)*nscale)/120.
templl(2,2)=latlon(2)
templl(2,1)=latlon(2)-real(lldim(2)*nscale)/120.

Do i=1,2
  If (templl(2,i).LT.-90.) templl(2,i)=-90.
  If (templl(2,i).GT.90.) templl(2,i)=90.
End Do

sermask=(tlld(:,:,1).GE.templl(1,1)).AND.(tlld(:,:,1).LE.templl(1,2)) &
        .AND.(tlld(:,:,2).GE.templl(2,1)).AND.(tlld(:,:,2).LE.templl(2,2))
sermask=sermask.AND.maskn

Select Case(mode)
  Case(0)
    ! Use all grid points
    sll=templl
    Return
  Case(1)
    sermask=sermask.AND.(grid.LE.scalelimit)
  Case(2)
    sermask=sermask.AND.(grid.GE.scalelimit)
  Case(3)
    sermask=sermask.AND.(grid.EQ.scalelimit)
  Case(4)
    ! Do nothing
  Case Default
    Write(6,*) 'ERROR: Internal error.  Unsupported mode in searchdim'
    Stop
End Select

If (.NOT.Any(sermask)) then
  sll=0.
  Return
End if

sll(1,2)=Maxval(tlld(:,:,1),sermask)
sll(1,1)=Minval(tlld(:,:,1),sermask)
sll(2,2)=Maxval(tlld(:,:,2),sermask)
sll(2,1)=Minval(tlld(:,:,2),sermask)

posll(1,2,:)=Maxloc(tlld(:,:,1),sermask)
posll(1,1,:)=Minloc(tlld(:,:,1),sermask)
posll(2,2,:)=Maxloc(tlld(:,:,2),sermask)
posll(2,1,:)=Minloc(tlld(:,:,2),sermask)
Do i=1,2
  Do j=1,2
    ! 1.6 is assumed to span from the centre to the corner (i.e., sqrt(2) if
    ! using a square grid)
    sll(i,j)=sll(i,j)+Real(j*2-3)*grid(posll(i,j,1),posll(i,j,2))*1.6/240.
  End Do
End Do

sll(1,1)=real(int((sll(1,1)-latlon(1))*120./real(nscale)))*real(nscale)/120.+latlon(1)
sll(1,2)=real(rndup((sll(1,2)-latlon(1))*120./real(nscale)))*real(nscale)/120.+latlon(1)
sll(2,1)=-real(rndup((latlon(2)-sll(2,1))*120./real(nscale)))*real(nscale)/120.+latlon(2)
sll(2,2)=-real(int((latlon(2)-sll(2,2))*120./real(nscale)))*real(nscale)/120.+latlon(2)

! Check bounds
Do i=1,2
  If (sll(i,1).LT.templl(i,1)) sll(i,1)=templl(i,1)
  If (sll(i,2).GT.templl(i,2)) sll(i,2)=templl(i,2)
End Do

! Consistancy check
If ((sll(1,1).GT.sll(1,2)).OR.(sll(2,1).GT.sll(2,2))) then
  sll=0.
End If

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the next scale to load.  The suboutine
! attempts to minimise the number of times the sibveg data is loaded.
!

Subroutine findsmallscale(nscale,scalelimit,latlon,llstore,grid,maskn,rlld,subsec,sll,sibsize,sibdim)

Implicit None

Integer, intent(in) :: sibsize,scalelimit
Integer, dimension(1:2), intent(in) :: sibdim
Integer, intent(inout) :: nscale
Integer, intent(out) :: subsec
Integer, dimension(1:2), intent(out) :: llstore
Real, dimension(1:sibdim(1),1:sibdim(2)), intent(in) :: grid
Real, dimension(1:sibdim(1),1:sibdim(2),1:2), intent(in) :: rlld
Real, dimension(1:2), intent(in) :: latlon
Real, dimension(1:2,1:2), intent(out) :: sll
Real tscale
Integer mode,maxscale,subsecmax
Integer findfact
Logical, dimension(1:sibdim(1),1:sibdim(2)), intent(in) :: maskn

tscale=Maxval(grid,maskn)

mode=1
If (nscale.EQ.999) mode=0

maxscale=Int(0.5*real(nscale)/Real(scalelimit))*scalelimit
maxscale=findfact(21600,maxscale,-scalelimit)
If (maxscale.LT.scalelimit) maxscale=scalelimit

llstore=(/ 43200/maxscale , 21600/maxscale /)
Call searchdim(mode,sll,maxscale,tscale,latlon,llstore,grid,maskn,rlld,sibdim)
Call scaleconvert(maxscale,subsecmax,llstore,sll,sibsize)

If (subsecmax.LT.1) Then
  Write(6,*) "WARN: Cannot locate unassigned points in findsmallscale"
  mode=0
  nscale=maxscale
Else
  nscale=Int(Minval(grid,maskn)/Real(scalelimit))*scalelimit
  nscale=findfact(21600,nscale,-scalelimit)
  If (nscale.LT.scalelimit) nscale=scalelimit
  subsec=subsecmax+1
  Do While (subsec.GT.subsecmax)
    ! Get estimate of array size
    llstore=(/ 43200/nscale , 21600/nscale /)
    ! Calculate domain for search
    Call searchdim(mode,sll,nscale,tscale,latlon,llstore,grid,maskn,rlld,sibdim)
    ! Define number of points in domain and subdivide into tiles if array is too big
    Call scaleconvert(nscale,subsec,llstore,sll,sibsize)
    If (subsec.GT.subsecmax) Then
      nscale=nscale+scalelimit
      nscale=findfact(21600,nscale,scalelimit)
    End If
  End Do
End If

If (nscale.GT.maxscale) nscale=maxscale
If (nscale.LT.scalelimit) nscale=scalelimit


llstore=(/ 43200/nscale , 21600/nscale /)
! Calculate domain for search
Call searchdim(mode,sll,nscale,tscale,latlon,llstore,grid,maskn,rlld,sibdim)
! Define number of points in domain and subdivide into tiles if array is too big
Call scaleconvert(nscale,subsec,llstore,sll,sibsize)

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This function finds the nearest factor of an integer
!

Integer Function findfact(x,y,delta)

Implicit None

Integer, intent(in) :: x,y,delta
Integer z

z=y

If (z.EQ.0) Then
  findfact=1
  Return
End If

Do While (Mod(x,z).NE.0.)
  z=z+delta
  If (z.LT.1) z=1
  If (z.GT.x) z=x
End Do

findfact=z

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate rsmin
! Better to calculate rsmin for each biome and combine in parallel.
! However, since the LAI is already averaged, we simply estimate
! rsmin for the averaged biome.
!

subroutine calrsmin(landdata,laidata,sibdim,mthrng,rdata)

implicit None

integer, intent(in) :: mthrng
integer, dimension(1:2), intent(in) :: sibdim
real, dimension(1:sibdim(1),1:sibdim(2),0:42), intent(in) :: landdata
real, dimension(1:sibdim(1),1:sibdim(2),0:mthrng-1), intent(in) :: laidata
real, dimension(1:sibdim(1),1:sibdim(2),1:mthrng), intent(out) :: rdata
real, dimension(1:42) :: rsunc
integer ilon,ilat,i
real nsum,tmprs

rsunc(1:13)=(/ 350.,300.,300.,300.,300., 230.,230.,230., 150., 230., 995., 150.,9900. /) ! CASA
rsunc(14:42)=(/ 370., 330., 260., 200., 150., 130., 200., 150., 110., 160., &
                100., 120.,  90.,  90.,  80.,  90., 150.,  80., 100.,  80., &
                 80.,  80.,  80.,  60.,  60., 120.,  80., 180., 995.  /)                  

rdata=0.

do ilat=1,sibdim(2)
  do ilon=1,sibdim(1)
    nsum=sum(landdata(ilon,ilat,1:41))
    if (nsum>0.) then
      tmprs=dot_product(landdata(ilon,ilat,1:41),1./rsunc(1:41))/nsum
      do i=1,mthrng
        rdata(ilon,ilat,i)=1./(tmprs*laidata(ilon,ilat,i-1))
      end do
    end if
  end do
end do

where (rdata<=0.)
  rdata=-1.
end where

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate roughness
! Better to calculate roughness for each biome and then combine drag
! at zmin.  However, since LAI is already averaged over biomes,
! we simply estimate the average roughness length assuming zo=0.1*h,
! and then adjust for LAI.
!

Subroutine calrough(landdata,laidata,sibdim,mthrng,rdata,zmin)

Implicit None

integer, intent(in) :: mthrng
Integer, dimension(1:2), intent(in) :: sibdim
real, intent(in) :: zmin
Real, dimension(1:sibdim(1),1:sibdim(2),0:42), intent(in) :: landdata
real, dimension(1:sibdim(1),1:sibdim(2),0:mthrng-1), intent(in) :: laidata
Real, dimension(1:sibdim(1),1:sibdim(2),1:mthrng), intent(out) :: rdata
real, dimension(1:42) :: xhc
Integer ilon,ilat,i
real nsum,tmplai,tmphc
real psih,rl,usuh,xx,dh,z0,z0h
real, parameter :: cr    = 0.3          ! element drag coefficient
real, parameter :: cs    = 0.003        ! substrate drag coefficient
real, parameter :: beta  = cr/cs        ! ratio cr/cs
real, parameter :: ccd   = 15.0         ! constant in d/h equation
real, parameter :: ccw   = 2.0          ! ccw=(zw-d)/(h-d)
real, parameter :: usuhm = 0.3          ! (max of us/uh)
real, parameter :: vonk  = 0.4          ! von karman constant
real, parameter :: zomin = 0.02

psih=alog(ccw)-1.0+1.0/ccw  ! i.e. .19315
xhc(1:13)=(/ 32.,20.,20.,17.,17., 1., 1., 1., 0.5, 0.6, 0.01, 1.,0.01 /) ! CASA
!xhc(1:13)=(/ 17.,35.,15.5,20.,19.3,0.6,0.6,7.,8.,0.6,0.5,0.6,6.,0.6,0.01,0.2,0.01 /) ! IGBP
xhc(14:42)=(/ 30.00, 28.00, 25.00, 17.00, 12.00, 10.00,  9.00,  7.00,  5.50,  3.00, &
               2.50,  2.00,  1.00,  0.60,  0.50,  0.50,  0.45,  0.75,  0.60,  0.45, &
               0.40,  0.60,  0.06,  0.24,  0.25,  0.35,  0.30,  2.50,  0.01  /)

rdata=0.

do ilat=1,sibdim(2)
  do ilon=1,sibdim(1)
    nsum=sum(landdata(ilon,ilat,1:41))
    if (nsum>0.) then
      tmphc=dot_product(landdata(ilon,ilat,1:41),1./log(zmin*10./xhc(1:41)))/nsum
      tmphc=zmin*10.*exp(-1./tmphc)
      do i=1,mthrng
        rl = 0.5*laidata(ilon,ilat,i-1)
        usuh   = min(sqrt(cs+cr*rl),usuhm)
        xx     = sqrt(ccd*max(rl,0.0005))  ! sqrt(7.5*rlai)
        dh     = 1.0 - (1.0 - exp(-xx))/xx ! .5*xx -.166*xx*xx + .
        z0h    = (1.0 - dh) * exp(psih - vonk/usuh)
        rdata(ilon,ilat,i)=max(zomin,tmphc*z0h)
        !tmplai=-0.0075*max(0.5,min(5.,laidata(ilon,ilat,i-1)))
        !rdata(ilon,ilat,i)=max(zomin,tmphc*(1.-0.91*exp(tmplai))) ! sellers 1996
      end do
    end if
  end do
end do

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate greeness (moved to CCAM indata.f)
!

!subroutine calgreen(laidata,sibdim,mthrng,rdata)
!
!implicit None
!
!integer, intent(in) :: mthrng
!integer, dimension(1:2), intent(in) :: sibdim
!real, dimension(1:sibdim(1),1:sibdim(2),0:mthrng-1), intent(in) :: laidata
!real, dimension(1:sibdim(1),1:sibdim(2),1:mthrng), intent(out) :: rdata
!
!rdata(:,:,:)=min(0.98,1.-exp(-0.6*laidata(:,:,:))) ! Sellers 1996 (see also Masson 2003)
!
!! do ilon=1,sibdim(1)
!!   do ilat=1,sibdim(2)
!!     efrac=landtype(ilon,ilat,1)                        ! evergreen
!!     tfrac=sum(landtype(ilon,ilat,2:5))                 ! non-evergreen trees
!!     gfrac=landtype(ilon,ilat,7)+landtype(ilon,ilat,12) ! grasses
!!     xfrac=1.-efrac-tfrac-gfrac                         ! mixed
!!     do i=1,mthrng
!!       rdata(ilon,ilat,i)=efrac*0.98+ &
!!                          tfrac*min(0.95,1.-exp(-0.5*lai(ilon,ilat,i-1)))+ &
!!                          gfrac*min(0.95,1.-exp(-0.6*lai(ilon,ilat,i-1)))+ &
!!                          xfrac*min(0.95,1.-0.5*exp(-0.5*lai(ilon,ilat,i-1))-0.5*exp(-0.6*lai(ilon,ilat,i-1)))
!!     end do
!!   end do
!! end do
!
!return
!end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate zobler soil texture from fao
!

subroutine calsoilnear(landdata,soildata,lsdata,sibdim,tdata)

implicit none

integer, dimension(1:2), intent(in) :: sibdim
integer, dimension(1:sibdim(1),1:sibdim(2)), intent(out) :: tdata
real, dimension(1:sibdim(1),1:sibdim(2),0:42), intent(in) :: landdata
real, dimension(1:sibdim(1),1:sibdim(2),0:8), intent(in) :: soildata
real, dimension(1:sibdim(1),1:sibdim(2)), intent(in) :: lsdata
integer ilon,ilat,pos(1),i

do ilon=1,sibdim(1)
  do ilat=1,sibdim(2)
    pos=Maxloc(landdata(ilon,ilat,1:41))
    if (1-nint(lsdata(ilon,ilat))==0) then
      if (landdata(ilon,ilat,0)>=landdata(ilon,ilat,42)) then
        tdata(ilon,ilat)=0 ! water
      else
        tdata(ilon,ilat)=-1 ! in-land water  
      end if
    else if (pos(1)==13) then
      tdata(ilon,ilat)=9 ! ice
    else
      pos=Maxloc(soildata(ilon,ilat,:))
      tdata(ilon,ilat)=pos(1)
    end if
  end do
end do

return
end

!subroutine calsoilusda(landdata,soildata,lsdata,sibdim,tdata)
!
!implicit none
!
!integer, dimension(1:2), intent(in) :: sibdim
!integer, dimension(1:sibdim(1),1:sibdim(2)), intent(out) :: tdata
!real, dimension(1:sibdim(1),1:sibdim(2),0:13), intent(in) :: landdata
!real, dimension(1:sibdim(1),1:sibdim(2),0:1), intent(in) :: soildata
!real, dimension(1:sibdim(1),1:sibdim(2)), intent(in) :: lsdata
!Real, dimension(12,7,2) :: vecpos
!Integer, dimension(12) :: vecnum
!Integer convert(1:12)
!integer ilon,ilat,pos(1),i,usdasoil,intersect
!integer calintersect
!
!vecpos=-1.
!
!vecnum(1)=3 ! Sand
!vecpos(1,1:3,1)=(/ 0.,   0.1,  0. /)
!vecpos(1,1:3,2)=(/ 0.85, 0.9,  1. /)
!vecnum(2)=4 ! Loam sand
!vecpos(2,1:4,1)=(/ 0.,   0.15, 0.1,  0.   /)
!vecpos(2,1:4,2)=(/ 0.7,  0.85, 0.9,  0.85 /)
!vecnum(3)=7 ! Sand loam
!vecpos(3,1:7,1)=(/ 0.,   0.07, 0.07, 0.2,  0.2,  0.15, 0.  /)
!vecpos(3,1:7,2)=(/ 0.5,  0.43, 0.52, 0.52, 0.8,  0.85, 0.7 /)
!vecnum(4)=6 ! Silt loam
!vecpos(4,1:6,1)=(/ 0.12, 0.27, 0.27, 0.,   0.,  0.12 /)
!vecpos(4,1:6,2)=(/ 0.,   0.,   0.23, 0.5,  0.2, 0.08 /)
!vecnum(5)=5 ! Loam
!vecpos(5,1:5,1)=(/ 0.07, 0.27, 0.27, 0.2,  0.07 /)
!vecpos(5,1:5,2)=(/ 0.43, 0.23, 0.45, 0.52, 0.52 /)
!vecnum(6)=5 ! Sand clay loam
!vecpos(6,1:5,1)=(/ 0.2,  0.27, 0.35, 0.35, 0.2 /)
!vecpos(6,1:5,2)=(/ 0.52, 0.45, 0.45, 0.65, 0.8 /)
!vecnum(7)=4 ! Silt clay loam
!vecpos(7,1:4,1)=(/ 0.27, 0.4,  0.4,  0.27 /)
!vecpos(7,1:4,2)=(/ 0.,   0.,   0.2,  0.2 /)
!vecnum(8)=4 ! Clay loam
!vecpos(8,1:4,1)=(/ 0.27, 0.4,  0.4,  0.27 /)
!vecpos(8,1:4,2)=(/ 0.2,  0.2,  0.45, 0.45 /)
!vecnum(9)=3 ! Sand clay
!vecpos(9,1:3,1)=(/ 0.35, 0.55, 0.35 /)
!vecpos(9,1:3,2)=(/ 0.45, 0.45, 0.65 /)
!vecnum(10)=3 ! Silt clay
!vecpos(10,1:3,1)=(/ 0.4,  0.6,  0.4 /)
!vecpos(10,1:3,2)=(/ 0.,   0.,   0.2 /)
!vecnum(11)=5 ! Clay
!vecpos(11,1:5,1)=(/ 0.4,  0.6,  1.,   0.55, 0.40 /)
!vecpos(11,1:5,2)=(/ 0.2,  0.,   0.,   0.45, 0.45 /)
!vecnum(12)=4 ! Silt
!vecpos(12,1:4,1)=(/ 0.,   0.12, 0.12, 0. /)
!vecpos(12,1:4,2)=(/ 0.,   0.,   0.08, 0.2 /)
!
!convert=(/ 1, 1, 4, 2, 4, 7, 2, 2, 5, 6, 3, 2 /)
!
!
!do ilon=1,sibdim(1)
!  do ilat=1,sibdim(2)
!    pos=Maxloc(landdata(ilon,ilat,1:13))
!    if (1-nint(lsdata(ilon,ilat)).eq.0) then
!      tdata(ilon,ilat)=0 ! water
!    else if (pos(1).eq.13) then
!      tdata(ilon,ilat)=9 ! ice
!    else
!      intersect=0
!      usdasoil=0
!      Do While ((Mod(intersect,2).EQ.0).and.(usdasoil.lt.12))
!        usdasoil=usdasoil+1
!        intersect=0
!        Do i=1,vecnum(usdasoil)-1
!          intersect=intersect+calintersect(vecpos(usdasoil,i,:),vecpos(usdasoil,i+1,:),soildata(ilon,ilat,:))
!        End Do
!        intersect=intersect+calintersect(vecpos(usdasoil,vecnum(usdasoil),:),vecpos(usdasoil,1,:),soildata(ilon,ilat,:))
!      End Do
!      tdata(ilon,ilat)=convert(usdasoil)
!    end if
!  end do
!end do
!
!return
!end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This function determines if the vertical soil vector poss intersects
! the line from posa to posb. Returns 0 if no intersection and 1 if
! there is an intersection
!

!Integer function calintersect(posa,posb,poss)
!
!Implicit None
!
!Real, dimension(2), intent(in) :: posa,posb,poss
!Real tempa,tempb,tempc,m,c,yi
!
!calintersect=0
!
!If (All(posa.EQ.poss).OR.All(posb.EQ.poss)) then
!  calintersect=1
!  Return
!End if
!
!If (posa(1).EQ.posb(1)) then
!  ! vertical line
!  tempa=(poss(2)-posa(2))*(poss(2)-posb(2))
!  If ((tempa.LE.0.).AND.(poss(1).EQ.posa(1))) calintersect=1
!  
!Else
!  ! non-vertical line
!  m=(posa(2)-posb(2))/(posa(1)-posb(1))
!  c=(posb(2)*posa(1)-posa(2)*posb(1))/(posa(1)-posb(1))
!  yi=m*poss(1)+c
!  
!  tempa=(yi-posa(2))*(yi-posb(2))
!  tempb=(poss(1)-posa(1))*(poss(1)-posb(1))
!  tempc=(yi-poss(2))*yi
!  
!  If ((tempa.LE.0.).AND.(tempb.LE.0.).AND.(tempc.LE.0.)) calintersect=1
!  
!End if
!
!! Special case
!If ((Sum(poss).GE.1.).AND.(Sum(posa).GE.1.).AND.(Sum(posb).GE.1.)) calintersect=0
!
!Return
!End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine patches LAI data from Lee, S. data
!

subroutine ozlai(sibdim,num,dataout,rlld,month)

use ccinterp

implicit none

integer, intent(in) :: num,month
integer, dimension(2), intent(in) :: sibdim
integer, dimension(sibdim(1),sibdim(2)) :: ncount
integer, dimension(2) :: lldim,pxy
integer ix,iy,lci,lcj,nface,ii,jj,imth,mthmn,mthmx
integer ni,nj,tix,tiy,k
real, dimension(sibdim(1),sibdim(2),0:num), intent(inout) :: dataout
real, dimension(sibdim(1),sibdim(2),2), intent(in) :: rlld
real, dimension(:,:,:), allocatable :: laiin
real, dimension(2,2,0:num) :: covertemp
real bx,by,bdelta,px,py,bmnx,bmxx,bmny,bmxy
real serlon,serlat,tbx,tby,tbdelta,alci,alcj
real indexlon,indexlat,ipol
character*2 cmth
logical, dimension(sibdim(1),sibdim(2)) :: ptreg,sermask

write(6,*) 'Apply LAI patch'

mthmn=month
mthmx=month
if (num.eq.11) then
  mthmn=1
  mthmx=12
end if

write(cmth,'(I2.2)') mthmn
write(6,*) 'Open aus_lai.'//cmth
open(30,file='aus_lai.'//cmth)
read(30,*) bx,by,bdelta,ix,iy
allocate(laiin(ix,iy,0:num))
close(30)
lldim(1)=ix
lldim(2)=iy

do imth=mthmn,mthmx
  write(cmth,'(I2.2)') imth
  write(6,*) 'Open aus_lai.'//cmth
  open(30,file='aus_lai.'//cmth)
  read(30,*) tbx,tby,tbdelta,tix,tiy
  if ((tix.ne.ix).or.(tiy.ne.iy).or.(tbx.ne.bx).or.(tby.ne.by).or.(tbdelta.ne.bdelta)) then
    write(6,*) "ERROR: LAI data has different dimensions for different months"
    stop
  end if
  read(30,*) laiin(:,:,imth-mthmn)
  close(30)
end do
  
bmnx=bx
bmxx=bx+bdelta*real(ix-1)
bmny=by
bmxy=by+bdelta*real(iy-1)

if (bmnx.lt.minval(rlld(:,:,1))) then
  bx=bx+360.
  bmnx=bx
  bmxx=bx+bdelta*real(ix-1)
end if
    
if (bmxx.gt.maxval(rlld(:,:,1))) then
  bx=bx-360.
  bmnx=bx
  bmxx=bx+bdelta*real(ix-1)
end if

ptreg=rlld(:,:,1).gt.bmnx
ptreg=ptreg.and.(rlld(:,:,1).lt.bmxx)
ptreg=ptreg.and.(rlld(:,:,2).gt.bmny)
ptreg=ptreg.and.(rlld(:,:,2).lt.bmxy)

where(ptreg)
  ncount=0
elsewhere
  ncount=1
end where
 
do imth=0,num
  where(ptreg)
    dataout(:,:,imth)=0.
  end where
end do

! bin
do ii=1,ix
  px=bx+bdelta*real(ii-1)
  do jj=1,iy
    py=by+bdelta*real(jj-1)
    call lltoijmod(px,py,alci,alcj,nface)
    lci = nint(alci)
    lcj = nint(alcj)
    lcj = lcj+nface*sibdim(1)
    if (ptreg(lci,lcj)) then
      dataout(lci,lcj,:)=dataout(lci,lcj,:)+laiin(ii,jj,:)
      ncount(lci,lcj)=ncount(lci,lcj)+1
    end if
  end do
end do
  
! interpolate
do lci=1,sibdim(1)
  do lcj=1,sibdim(2)
    if ((ncount(lci,lcj).eq.0).and.ptreg(lci,lcj)) then
      px=rlld(lci,lcj,1)
      py=rlld(lci,lcj,2)
     
      serlon=indexlon(px,bx,120.*bdelta)
      serlat=indexlat(py,by,120.*bdelta)

      ii=int(serlon)
      jj=int(serlat)

      if ((ii.GE.1).AND.(ii.LT.ix).AND.(jj.GE.1).AND.(jj.LT.iy)) Then
	      
        ni=ii+1 ! 4 point interpolation
        nj=jj+1
	      
        serlon=serlon-real(ii)
        serlat=serlat-real(jj)
	
        call realfill(covertemp,laiin,lldim,ii,ni,jj,nj,num)
   
        do k=0,num
          dataout(lci,lcj,k)=ipol(covertemp(:,:,k),serlon,serlat)
        end do
        ncount(lci,lcj)=1
      end if
        
    end if
  end do
end do

! near nbr
sermask(:,:)=ncount(:,:).gt.0
if (Any(sermask)) then
  do lci=1,sibdim(1)
    do lcj=1,sibdim(2)
      if (ncount(lci,lcj).eq.0.and.ptreg(lci,lcj)) then
        call findnear(pxy,lci,lcj,sermask,rlld,sibdim)
        dataout(lci,lcj,:)=dataout(pxy(1),pxy(2),:)
        ncount(lci,lcj)=1
      end if
    end do
  end do
else
  write(6,*) 'WARN: Cannot find any non-trivial near nbrs'
  write(6,*) '      Assume data is trivial'
  dataout=0.
  ncount=1
end if	

deallocate(laiin)

do imth=0,num
  where(ptreg.and.ncount.gt.0)
    dataout(:,:,imth)=dataout(:,:,imth)/real(ncount)
  end where
  where(dataout(:,:,imth).le.0..and.ptreg)
    dataout(:,:,imth)=0.1
  end where
end do

return
end
