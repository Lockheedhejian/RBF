subroutine RBF
    use Define_Variables
    
    
     nbc=0

    do ii=1,Edge_Count
    do k=1,Wall_Points_Number(ii,3)
    do j=1,Wall_Points_Number(ii,2)
    do i=1,Wall_Points_Number(ii,1)

            nbc=nbc+1
            Points_Position_RBF(nbc,1)=Wall_Points_Position(ii,i,j,k,1)
            Points_Position_RBF(nbc,2)=Wall_Points_Position(ii,i,j,k,2)
            kpoint_RBF(nbc,1)=ii
            kpoint_RBF(nbc,2)=i
            kpoint_RBF(nbc,3)=j
    enddo
    enddo
    enddo
    enddo

   NBC_3=nbc+3
   Radius=1.0
    
    allocate (C_Matrix(nbc3,nbc3),A_Matrix(nbc3,nbc3))
	allocate (B_Matrix(nbc3),x(nbc3))
   
   
    do i=1,nbc
        x1=Points_Position_RBF(i,1)
        y1=Points_Position_RBF(i,2)
        do j=i+1,nbc
            x2=Points_Position_RBF(j,1)
            y2=Points_Position_RBF(j,2)
            ss=sqrt((x2-x1)**2+(y2-y1)**2)
            ss=ss/radius
            
            if(ss.le.1.0e-6) then
                tmp=0.0                
                    do k=j+1,nbc
                        Points_Position_RBF(k-1,1)=Points_Position_RBF(k,1)
                        Points_Position_RBF(k-1,2)=Points_Position_RBF(k,2)
                        kpoint_RBF(k-1,1)=kpoint_RBF(k,1)
                        kpoint_RBF(k-1,2)=kpoint_RBF(k,2)
                        kpoint_RBF(k-1,3)=kpoint_RBF(k,3)
                    enddo   
                nbc=nbc-1
            else               
                tmp=ss**3
				!tmp=ss**2*log10(ss)
                !径向基函数的两种不同的形式
            endif
            c(i,j)=tmp
        enddo
        C_Matrix(i,nbc3-2)=1.0
        C_Matrix(i,nbc3-1)=x1
        C_Matrix(i,nbc3)=y1
    enddo
    
    do i=nbc3-2,nbc3
    do j=i,nbc3
        C_Matrix(i,j)=0.0
    enddo
    enddo

    do i=1,nbc3
    do j=i,nbc3
        C_Matrix(j,i)=C_Matrix(i,j)
    enddo
    enddo
    
    end
    
    


   