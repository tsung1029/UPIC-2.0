      module mextdriv3
      
      use in3
      use minit3
      implicit none
      private
      public :: donothing,laguerregaussiandriv,gaussiandriv
      public :: planewavedriv

      contains

      subroutine donothing()
      implicit none
      return
      end subroutine

      subroutine laguerregaussiandriv(exyze, t, nx, nxe, ny, nypmx, nz,&
              &nzpmx, nvpy, nvpz, idproc)
      implicit none
      integer :: nx,nxe,ny,nypmx,nz,nzpmx,nvpy,nvpz,idproc
      real :: t
      real, dimension(:,:,:,:) :: exyze
      integer :: i,j,k,yproc,zproc
      integer,dimension(:,:,:),allocatable,save :: proc_pos
      real :: xpos, ypos, zpos, wavek
      real :: tempamp, tfac
      real :: r2, spot_size2,phi_l,l_sum,l_diff,alpha

      wavek = 6.283185307/real(nx)
      tfac = 1.

      if (.not. allocated(proc_pos)) then
              allocate(proc_pos(nvpy,nvpz,2))
              do i=0, nvpy-1
                   do j=0, nvpz-1
                        proc_pos(i+1,j+1,1) = ny / nvpy
                        proc_pos(i+1,j+1,1) = proc_pos(i+1,j+1,1)*i
                        proc_pos(i+1,j+1,2) = nz / nvpz
                        proc_pos(i+1,j+1,2) = proc_pos(i+1,j+1,2)*j
                   enddo
              enddo
      endif

      if (t < timerise) then
              tfac = t / timerise
      else if (t < timerise + timeflat + timefall) then
              tfac = 1. - (t - (timerise+timeflat))/timefall
      endif

      if (t < timerise + timeflat + timefall) then
              yproc = mod(idproc,nvpy) + 1
              zproc = idproc / nvpy + 1

              do k=1,nzpmx
                   do j=1,nypmx
                         ypos = real(proc_pos(yproc,zproc,1)+j-ny/2)
                         zpos = real(proc_pos(yproc,zproc,2)+k-nz/2)
                         r2 = ypos ** 2 + zpos ** 2
                         l_diff = real(abs(l1 - l2))
                         phi_l = l_diff * atan2(zpos,ypos)
                         l_sum = real(l1+l2)
                         spot_size2 = spot_size ** 2
                         
                         if (r2 /= 0) then
                                 tempamp = amp * tfac * (2. *&
                                & r2/spot_size2) ** (l_sum / 2)&
                                & * exp(-2.*r2/spot_size2)
                         else
                                 tempamp = amp * tfac
                         endif

                         !if (r2 > 4 * spot_size2) then
                          !       tempamp = tempamp * &
                           !     &exp(-r2**2/spot_size2**2)
                         !endif
                                                  
                         !if (l_sum /= 0) then
                          !      tempamp = tempamp / wavek * (2/l_sum)&
                           !    & **(l_sum / 2) * exp(l_sum / 2)
                         !else
                         !       tempamp = tempamp / wavek
                         !endif
                        
                         do i=2,nx
                                xpos = real(i) - 2.
                                alpha=phi_l-wavek*xpos+wavew*t
                                if (i == 16) then
                                        print*, exyze(1,i,j,k)
                                endif
                                exyze(1,i,j,k) = exyze(1,i,j,k)&
                               &-wavek*tempamp*sin(alpha)
                                if (i == 16) then
                                        print*, exyze(1,i,j,k)
                                endif
                                if (r2 /= 0) then
                                        exyze(2,i,j,k) =exyze(2,i,j,k) &
                                        &+tempamp * (-zpos*l_diff /&
                                        &r2*sin(alpha)+ypos*(-l_sum/&
                                        &r2+4./spot_size2)*cos(alpha))
                                        
                                        exyze(3,i,j,k)=exyze(3,i,j,k) &
                                        &+tempamp * (ypos*l_diff/r2 * &
                                        &sin(alpha) + zpos*(-l_sum/r2+&
                                        &4./spot_size2)*cos(alpha))
                                 endif
                         enddo
                   enddo
              enddo
      endif

      end subroutine

      subroutine gaussiandriv(exyze,t, nx,nxe,ny,nypmx,nz,nzpmx,nvpy,&
                &nvpz,idproc)
      implicit none
      integer :: nx,nxe,ny,nypmx,nz,nzpmx,nvpy,nvpz,idproc
      real :: t
      real, dimension(:,:,:,:) :: exyze
      integer :: i,j,k,yproc,zproc
      integer,dimension(:,:,:),allocatable,save :: proc_pos
      real :: xpos, ypos, zpos, wavek
      real :: tempamp, tfac
      real :: r2, spot_size2,alpha

      wavek = 6.283185307/real(nx)
      tfac = 1.

      if (.not. allocated(proc_pos)) then
              allocate(proc_pos(nvpy,nvpz,2))
              do i=0, nvpy-1
                   do j=0, nvpz-1
                        proc_pos(i+1,j+1,1) = ny / nvpy
                        proc_pos(i+1,j+1,1) = proc_pos(i+1,j+1,1)*i
                        proc_pos(i+1,j+1,2) = nz / nvpz
                        proc_pos(i+1,j+1,2) = proc_pos(i+1,j+1,2)*j
                   enddo
              enddo
      endif

      if (t < timerise) then
              tfac = t / timerise
      else if (t < timerise + timeflat + timefall) then
              tfac = 1. - (t - (timerise+timeflat))/timefall
      endif

      if (t < timerise + timeflat + timefall) then
              yproc = mod(idproc,nvpy) + 1
              zproc = idproc / nvpy + 1

              do k=1,nzpmx
                   do j=1,nypmx
                         ypos = real(proc_pos(yproc,zproc,1)+j-ny/2)
                         zpos = real(proc_pos(yproc,zproc,2)+k-nz/2)
                         r2 = ypos ** 2 + zpos ** 2
                         spot_size2 = spot_size ** 2
                         
                         tempamp = amp*tfac*exp(-2.*r2/spot_size2)

                         do i=2,nx
                                xpos = real(i) - 2.
                                alpha=-wavek*xpos+wavew*t
                                exyze(1,i,j,k) = exyze(1,i,j,k)&
                               &-wavek*tempamp*sin(alpha)
                         enddo
                   enddo
              enddo
      endif

      end subroutine

      subroutine planewavedriv(exyze,t,nx,nxe,ny,nypmx,nz,nzpmx,nvpy&
                &,nvpz,idproc)
      implicit none
      integer :: nx,nxe,ny,nypmx,nz,nzpmx,nvpy,nvpz,idproc
      real :: t
      real, dimension(:,:,:,:) :: exyze
      integer :: i
      real :: xpos, wavek
      real :: tempamp, tfac
      real :: alpha

      wavek = 6.283185307/real(nx)
      tfac = 1.

      if (t < timerise) then
              tfac = t / timerise
      else if (t < timerise + timeflat + timefall) then
              tfac = 1. - (t - (timerise+timeflat))/timefall
      endif

      if (t < timerise + timeflat + timefall) then
              tempamp = tfac * amp     
              do i=2,nxe
                      xpos = real(i) - 2.
                      alpha = wavew*t-wavek*xpos
                      exyze(1,i,:,:)=-wavek*tempamp*sin(alpha)
              enddo
      endif

      end subroutine

      end module
