subroutine TFI_Parameter
!---------------------------------------------------------------------------
!该子程序用来计算每个网格边的长度，以及TFI方法中的alpha和beta
    use Define_Variables
    implicit double precision(A_Matrix-h,o-z)

    do ii=1,Num_Block
    do k=1,num(ii,3)
    do j=1,num(ii,2)
    do i=1,num(ii,1)
        alpha=0.0
        alpha0=0.0
        beta=0.0
        beta0=0.0
        do i1=1,num(ii,1)-1
            alpha0=alpha0+Edge_distance(ii,i1,j,k,2)
            if(i1.lt.i) alpha=alpha+Edge_distance(ii,i1,j,k,2)
        enddo
        
        do j1=1,num(ii,2)-1
            beta0=beta0+Edge_distance(ii,i,j1,k,1)
            if(j1.lt.j) beta=beta+Edge_distance(ii,i,j1,k,1)
        enddo
               
        alpha=alpha/alpha0
        beta=beta/beta0
        coefficient(ii,i,j,k,2)=alpha
        coefficient(ii,i,j,k,1)=beta
        
    enddo
    enddo
    enddo
    enddo
    

    return
    end