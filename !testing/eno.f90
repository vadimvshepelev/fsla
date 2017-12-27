	subroutine ninter(ii,du,ni,n)
	dimension u(50), ii(50), du(4,50)
c	****************************
c	nonoscillatory interpolation
c	****************************
	do 30 m=2, ni
	do 35 i=1, n-1
35	du(m,i)=du(m-1,i+1)-du(m-1,i)
	du(m,n)=du(m-1,  1)-du(m-1,n)
	do 40 i=1,n-1
	i0=ii(i)
	ip=i0
	if(ip.le.0)ip=ip+n
	im=i0-1
	if(im.le.0)im=im+n
	ii(i)=i0+imn(du(m,im),du(m,ip))
40  continue
30	continue
	return 
	end
	function imn(x,y)
	imn=0
	if(abs(x).le.abs(y))imn=-1
	return
	end
