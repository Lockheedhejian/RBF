!---------------------------------------------------------------------------
    subroutine geometry
!---------------------------------------------------------------------------
!计算网格的面积
!---------------------------------------------------------------------------
    use Define_Variables
    implicit double precision(A_Matrix-h,o-z)

    areatotal=0.d0
    
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
        xl11=x3-X_Min
        xl12=y3-Y_Min
        xl21=x4-x2
        xl22=y4-y2
        area(ii,i,j,k)=0.5*(xl11*xl22-xl12*xl21)
        if(area(ii,i,j,k).lt.0.0)then
            area(ii,i,j,k)=-area(ii,i,j,k)
        endif
        if(area_min.gt.abs(area(ii,i,j,k))) area_min=abs(area(ii,i,j,k))
        if(area_max.lt.abs(area(ii,i,j,k))) area_max=abs(area(ii,i,j,k))
        
        vv=(abs(area(ii,i,j,k))-area_min)/(area_max-area_min)
        area_1(ii,i,j,k)=1.0-vv**2
		areatotal=areatotal+area(ii,i,j,k)
        
  !---------------------------------------------------------------------------------------------
 !       计算网格边的长度
  !----------------------------------------------------------------------------------------------
        Edge_distance(ii,i,j,k,1)=sqrt((x4-X_Min)**2+(y4-Y_Min)**2)
        Edge_distance(ii,i,j,k,2)=sqrt((x2-X_Min)**2+(y2-Y_Min)**2)
        
        if(i.eq.(num(ii,1)-1))then
            do jx=1,num(ii,2)-1
            Edge_distance(ii,num(ii,1),jx,k,1)=sqrt((Points_Position_initial(ii,num(ii,1),jx+1,k,1)&
                                                 -Points_Position_initial(ii,i,jx,k,1))**2+(Points_Position_initial(ii,num(ii,1),jx+1,k,2)-Points_Position_initial(ii,i,jx,k,2))**2)
            end do
        end if
               
        if(j.eq.(num(ii,2)-1))then
            do ix=1,num(ii,1)-1
            Edge_distance(ii,ix,num(ii,2),k,1)=sqrt((Points_Position_initial(ii,ix+1,num(ii,2),k,1)&
                                                -Points_Position_initial(ii,ix,j,k,1))**2+(Points_Position_initial(ii,ix+1,num(ii,2),k,2)-Points_Position_initial(ii,ix,j,k,2))**2)
            end do
        end if
        
        
    !---------------------------------------------------------------------    

        
    enddo
    enddo
    enddo
    enddo

	

    do ii=1,Num_Block
    do k=1,num(ii,3)
    do j=1,num(ii,2)-1
    do i=1,num(ii,1)-1
        vv=(abs(area(ii,i,j,k))-area_min)/(area_max-area_min)
        area_1(ii,i,j,k)=1.0-vv**2
		areatotal=areatotal+area(ii,i,j,k)
    enddo
    enddo
    enddo
    enddo

    end