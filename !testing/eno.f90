	program ENO
	dimension ii(50,50),jj(50,50),du(4,50),u(50,50),u1d(10),ii1d(10),kk1d(10),du1d(4,10)
	ni = 10
	n  = 10
	nx = 50
	ny = 50
	print *, "Hello from ENO code!!!"
	
	do 40 i=1,50
		do 30 j=1,50
			if(i.le.25) u(i,j)=1.  
			if(i.gt.25) u(i,j)=.125	
30		continue
40  continue

	do 45 i=1,10
		if(i.le.5) u1d(i)=1.  
45		if(i.gt.5) u1d(i)=.125	
	do 47 i=1,10
		du1d(1,i)=u1d(i)
		kk1d(i)=i
47		ii1d(i)=i
	print *,"ii1d="
	print *, ii1d 
	print *,"u1d="
	print *, u1d 
	print *, "Starting 'ninter' procedure"

	call ninter(kk1d,du1d,4,10)
	do 48 i=1,10
48		ii1d(i)=i-kk1d(i)
		

	print *,"du1d="
	print *, du1d 
	print *,"ii1d="
	print *, ii1d 

	call rctble	


!	print *,"u="
!	do 50 i=1,50
!		print *,u(i,10)
!50	continue
!	print *,"ii(1,10)=",ii(1,10)
!	print *,"ii(10,1)=",ii(10,1)
!	print *,"jj(1,10)=",jj(1,10)
!	print *,"jj(10,1)=",jj(10,1)
!	print *,"ii(10,1)=",ii(10,1)
!	print *,"jj(10,1)=",jj(10,1)

!	call stencil(u,ii,jj,nx,ny)

!print *,"u="
!	do 60 i=1,50
!		print *,u(i,10)
!60	continue
!	print *,"ii(1,10)=",ii(1,10)
!	print *,"ii(10,1)=",ii(10,1)
!	print *,"jj(1,10)=",jj(1,10)
!	print *,"jj(10,1)=",jj(10,1)
!	print *,"jj(1,10)=",jj(1,10)
!	print *,"ii(10,1)=",ii(10,1)
!	print *,"jj(10,1)=",jj(10,1)


	end program ENO
	

!	Procedure 1 is a FORTRAN subprogram to determine the best 1-d stencil. It is
!	invoked repeatedly, in Procedure 2, to compute the best stencil in each direction. We will
!   see in the next section how such Stencils may be used in 2-d reconstructions.

!	Procedure 3 presents a FORTRAN subprogram to develop a table of coefficients to be
!	used to efficiently compute two-dimensional reconstruction at the center of a cell. 
!	Procedure 4 presents a subprogram to perform 2-d reconstructions using the stencils generated
!	in Procedures 1 and 2. 



!	Procedure 1
	subroutine ninter(ii,du,ni,n)
!   На входе:
!   ii -- одномерный массив из целых точек
!   du -- двумерный массим разностей, первая размерность -- номер разности, вторая -- количество точек
!   ni -- некоторый целый параметр (в вызове из stencil ni=4 -- порядок интерполяции?)
!   n  -- некоторый целый параметр (в вызове из stencil это переменные nx и ny, видимо, это количество ячеек)
!   На выходе:
!   в ii пишет какие-то номера (номера шаблонов?)
!	в du пишет разности
	 

	dimension u(50), ii(50), du(4,50)
!	****************************
!	nonoscillatory interpolation
!	****************************

!!!
!	print *,"ii="
!	print *,ii
!	print *,"du="
!	print *,du
!	print *,"ni=",ni,"        n=",n
!!!



	do 30 m=2, ni
		do 35 i=1, n-1
35		du(m,i)=du(m-1,i+1)-du(m-1,i)
	
		!du(m,n)=du(m-1,  1)-du(m-1,n) ! Периодические гран. условия, вычитает из начала
		du(m,n)=du(m-1,  n)-du(m-1,n) ! Прозрачные гран. условия, вычитаем из того же значения

		do 40 i=1,n-1
			i0=ii(i)
			ip=i0
			if(ip.le.0)ip=ip+n ! Периодические гран. условия, прибавляет к числу точек в сетке
			im=i0-1
			if(im.le.0)im=im+n ! Периодические гран. условия, прибавляет к числу точек в сетке
			ii(i)=i0+imn(du(m,im),du(m,ip))	! Если первая разность строго больше второй, то ноль, если меньше-равна, то -1
40		continue
30	continue	
	!print *, ii(1)
	!print *, uu(1)
	return 
	end
	
	function imn(x,y)
	imn=0
	if(abs(x).le.abs(y))imn=-1
	return
	end



!	Procedure 2
	subroutine stencil(u,ii,jj,nx,ny)
	dimension u(50,50),ii(50,50),jj(50,50)
	dimension du(4,50),kk(50)

!!!!
	print *, "Starting 'stencil' procedure"

!!!!


	do 10 i=1,nx
		do 5 j=1,ny
			du(1,j)=u(i,j)
5			kk(j)=j
		call ninter(kk,du,4,ny)
		
		do 7 j=1,ny
7		jj(i,j)=j-kk(j)
10  continue

	do 20 j=1,ny
		do 25 i=1,nx
			du(1,i)=u(i,j)
25			kk(i)=i
		call ninter(kk,du,4,nx)
		
		do 27 i=i,nx
27		ii(i,j)=i-kk(i)
20	continue

	return
	end



!	Procedure 3
	subroutine rctble
	common/rc/dlprc(4,4)
	al=0.5
	do 20 k=1,4
	zk=k+al-1.
	do 20 j=1,4
	pr=1.
	sm=0.
	do 15 m=0,4
	if(m.eq.j)goto 15
	anw=zk-m
	pr=pr*anw/float(j-m)
	sm=sm+1./anw
15	continue
20	dlprc(j,k)=pr*sm
	write (*,*) 'dlprc'
	do 30 j=1,4		
	write (*,*) (dlprc(j,k),k=1,4)
30	continue
	return
	end							  ;



!	Procedure 4
	subroutine rcnst(u,uu,nx,ny)
	common/rc/dlprc(4,4)
	common/d/a,b,almdx,almdy
	dimension u(50,50),uu(50,50)
	dimension wgs(50,50)
	dimension ii(50,50),jj(50,50)
	call stencil(u,ii,jj,nx,ny)
	do 10 i=1,nx
	do 10 j=1,ny
	kj=jj(i,j)+1
	jbg=j-kj
	do 10 l=1,2
	sm=0.
	smu=0.
	do 5 jlg=1,4
	jv=jbg+jlg
	if(jv.gt.ny)jv=jv-ny
	if(jv.lt.1)jv=jv+ny
	smu=smu+u(i,jv)
5	sm=sm+smu*dlprc(jlg,kj)
10  wgs(i,j)=sm
	do 20 i=1,nx
	do 20 j=1,ny
	ki=ii(i,j)+1
	ibg=i-ki
	sm=0.
	smu=0.
	do 15 ilg=1,4
	iv=ibg+ilg
	if(iv.gt.nx)iv=iv-nx
	if(iv.lt.1)iv=iv+nx
	smu=smu+wgs(iv,j)
15  sm=sm+smu*dlprc(ilg,ki)
	uu(i,j)=sm
20  continue
	return
	end