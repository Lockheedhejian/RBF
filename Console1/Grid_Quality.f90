subroutine Grid_Quality
    

    use Define_Variables
    implicit double precision(a-h,o-z)
    call cpu_time(time_begin)
    
!---------------------------------------------------------------
    !该子程序用来计算网格的质量
!---------------------------------------------------------------   
    do ii=1,Num_Block
    do k=1,num(ii,3)
    do j=1,num(ii,2)-1
    do i=1,num(ii,1)-1
        X_Min=Points_Position_initial(ii,i,j,k,1)
        x2=Points_Position_initial(ii,i+1,j,k,1)
        x3=Points_Position_initial(ii,i+1,j+1,k,1)
        x4=Points_Position_initial(ii,i,j+1,k,1)
        Y_Min=Points_Position_initial(ii,i,j,k,2)
        y2=Points_Position_initial(ii,i+1,j,k,2)
        y3=Points_Position_initial(ii,i+1,j+1,k,2)
        y4=Points_Position_initial(ii,i,j+1,k,2)
        xl1x=x2-X_Min
        xl1y=y2-Y_Min
        xl2x=x3-x2
        xl2y=y3-y2
        xl3x=x4-x3
        xl3y=y4-y3
        xl4x=X_Min-x4
        xl4y=Y_Min-y4
        cd1=sqrt(xl1x**2+xl1y**2)
        cd2=sqrt(xl2x**2+xl2y**2)
        cd3=sqrt(xl3x**2+xl3y**2)
        cd4=sqrt(xl4x**2+xl4y**2)
        dj1=-(xl1x*xl4x+xl1y*xl4y)/(cd1*cd4)
        dj2=-(xl2x*xl1x+xl2y*xl1y)/(cd2*cd1)
        dj3=-(xl3x*xl2x+xl3y*xl2y)/(cd3*cd2)
        dj4=-(xl4x*xl3x+xl4y*xl3y)/(cd4*cd3)
        Q1=acos(dj1)*180.0/pi
        Q2=acos(dj2)*180.0/pi
        Q3=acos(dj3)*180.0/pi
        Q4=acos(dj4)*180.0/pi
        Qmin=min(Q1,Q2,Q3,Q4)
        Qmax=max(Q1,Q2,Q3,Q4)
        quality(ii,i,j,k)=max(((Qmax-Qe)/(180.0-Qe)),((Qe-Qmin)/Qe))     
    
    enddo
    enddo
    enddo
    enddo
    
    do ii=1,Num_Block
    do k=1,num(ii,3)
    do j=1,num(ii,2)-1
    do i=1,num(ii,1)-1
        vv=area_1(ii,i,j,k)
        ave_quality_new=ave_quality_new+quality(ii,i,j,k)*vv
        sum=sum+vv
    enddo
    enddo
    enddo
    enddo
    !write(*,*)ave_quality_new/sum
    call cpu_time(time_end)
    Grid_Quality_Time=time_end-time_begin
    
    end
