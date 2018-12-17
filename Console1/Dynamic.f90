subroutine dynamic_b
!---------------------------------------------------------------------------

    use Define_Varibles
    implicit double precision(a-h,o-z)

	character	::line*240
	allocatable	::nump(:),x0(:,:),z0(:,:),xdef(:,:),zdef(:,:),dx(:,:),dz(:,:)
	allocatable	::qmp(:,:),qabx(:),qabz(:)

	nm=3

!	读入翼型离散点位移

	allocate(nump(nm))
	open(10,file="initial_flap.dat")
	read(10,*)	nump(1)!读入点的个数
	close(10)
	open(10,file="initial_main.dat")
	read(10,*)	nump(2)
	close(10)
	if(nm.eq.3)	then
	open(10,file="initial_slat.dat")
	read(10,*)	nump(3)
	close(10)
	endif

	allocate(x0(nm,maxval(nump)),z0(nm,maxval(nump)))
	allocate(xdef(nm,maxval(nump)),zdef(nm,maxval(nump)))
	allocate(dx(nm,maxval(nump)+3),dz(nm,maxval(nump)+3))

	open(10,file="initial_flap.dat")
	read(10,*)	nump(1)
	do i=1,nump(1)
	read(10,*)	x0(1,i),z0(1,i)!读出每个点的位置坐标
	enddo
	close(10)

	open(10,file="initial_main.dat")
	read(10,*)	nump(2)
	do i=1,nump(2)
	read(10,*)	x0(2,i),z0(2,i)
	enddo
	close(10)

	open(10,file="def_flap.dat")
	read(10,*)	nump(1)
	do i=1,nump(1)
	read(10,*)	xdef(1,i),zdef(1,i)
	enddo
	close(10)

	open(10,file="def_main.dat")
	read(10,*)	nump(2)
	do i=1,nump(2)
	read(10,*)	xdef(2,i),zdef(2,i)
	enddo
	close(10)

!	三段翼情形
	if(nm.eq.3)	then
	open(10,file="initial_slat.dat")
	read(10,*)	nump(3)
	do i=1,nump(3)
	read(10,*)	x0(3,i),z0(3,i)
	enddo
	close(10)

	open(10,file="def_slat.dat")
	read(10,*)	nump(3)
	do i=1,nump(3)
	read(10,*)	xdef(3,i),zdef(3,i)
	enddo
	close(10)
	endif

	do n=1,nm
	do i=1,nump(n)
	dx(n,i)=xdef(n,i)-x0(n,i)!计算每个方向上的变形量
	dz(n,i)=zdef(n,i)-z0(n,i)
	enddo
	enddo

    do ii=1,nmp
    do k=1,number(ii,2)
    do i=1,number(ii,1)
        deformation(ii,i,k,1)=0.0
        deformation(ii,i,k,2)=0.0
    enddo
    enddo
    enddo

	allocate(qmp(maxval(nump)+3,maxval(nump)+3),qabx(maxval(nump)+3),qabz(maxval(nump)+3))

	do n=1,nm

!	获取插值系数矩阵

!	检查相邻近点
	dd=1.0e-6
	nnn0=nump(n)     !每段翼型表面的点数
    do i=1,nnn0-1
	  j0=i 
	  lgd=0 
	  do while(j0.le.nnn0)
	    if(lgd.eq.1) then       !当两点离得较近时，lgd会变成1，从而使得j0的值保持不变
	      j0=j0
	    else
		  j0=j0+1
	    endif
	    if(j0.gt.nnn0) exit
		x1=x0(n,i)
        y1=z0(n,i)
        x2=x0(n,i+1)
        y2=z0(n,i+1)
	    d=sqrt((x2-x1)**2+(y2-y1)**2)
	    if(d.lt.dd) then		  	
		  do k=j0,nnn0-1
              x0(n,k)=x0(n,k+1)
              z0(n,k)=z0(n,k+1)
              dx(n,k)=dx(n,k+1)
              dz(n,k)=dz(n,k+1)
		  end do
		  nnn0=nnn0-1
	      lgd=1
		else
	      lgd=0
		end if
	  end do
	end do
	nin=nnn0

    do i=nin+1,nin+3        !因为b的系数矩阵是四阶，丁力论文13页
    dx(n,i)=0.
    dz(n,i)=0.
    enddo

	qmp=0.d0

!	RBF method

	do j=1,nin
	do i=1,nin
	rij=sqrt((x0(n,i)-x0(n,j))**2+(z0(n,i)-z0(n,j))**2)!丁力论文12页
	qmp(i,j)=rij**3
	enddo
	enddo

	do i=1,nin
	qmp(i,nin+1)=1.0d0
	qmp(i,nin+2)=x0(n,i)
	qmp(i,nin+3)=z0(n,i)
	qmp(nin+1,i)=1.0d0
	qmp(nin+2,i)=x0(n,i)
	qmp(nin+3,i)=z0(n,i)
	enddo 

	epsinv=1.d-20
	call dinvr(maxval(nump)+3,nin+3,epsinv,qmp)
	
	qabx=0.d0;qabz=0.d0
	do i=1,nin+3
		do j=1,nin+3
		qabx(i)=qabx(i)+qmp(i,j)*dx(n,j)!求出的是系数矩阵中a矩阵
		qabz(i)=qabz(i)+qmp(i,j)*dz(n,j)!求出的是系数矩阵中b矩阵
		enddo
	enddo

	do ii=1,nmp

	im=number(ii,1)-1
	km=number(ii,2)-1
	
	do k=1,km+1
	do i=1,im+1
	if(kboundary(ii,i,k).eq.-1.and.kwalltype(ii,i,k).eq.n)	then
	xx=coord_initial(ii,i,k,1);zz=coord_initial(ii,i,k,2)
	deformation(ii,i,k,1)=0.d0;deformation(ii,i,k,2)=0.d0
		do j=1,nin
		rij=sqrt((xx-x0(n,j))**2+(zz-z0(n,j))**2)
		deformation(ii,i,k,1)=deformation(ii,i,k,1)+qabx(j)*rij**3
		deformation(ii,i,k,2)=deformation(ii,i,k,2)+qabz(j)*rij**3
        
		enddo
		deformation(ii,i,k,1)=deformation(ii,i,k,1)+qabx(nin+1)+qabx(nin+2)*xx+qabx(nin+3)*zz
		deformation(ii,i,k,2)=deformation(ii,i,k,2)+qabz(nin+1)+qabz(nin+2)*xx+qabz(nin+3)*zz
	endif
	enddo
	enddo

	enddo

	enddo

	deallocate(qmp,qabx,qabz)
	deallocate(nump,x0,z0,xdef,zdef,dx,dz)
    return
    end