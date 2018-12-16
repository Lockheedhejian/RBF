subroutine Read_Boundary_Condition
    
    use Define_Variables
    implicit double precision(A_Matrix-h,o-z)
    
    character::temp*30

    
    call cpu_time(time_begin)
    !call Cal_WallNumber
!---------------------------------------------------------------
    !读取网格边界信息
!---------------------------------------------------------------

    open(1,file="BC_Type.fvbnd")
    
    read(1,*)

    do ii=1,11
        read(1,*)temp
        !write(*,*)temp
        if(temp=="BOUNDARIES")then 
            go to 11  
        end if      
    end do
11      continue    
       
        ii=0
13       continue 
        read(1,*,end=14)BC_Temp(1),BC_Temp(2),BC_Temp(3),BC_Temp(4),BC_Temp(5),BC_Temp(6),BC_Temp(7),BC_Temp(8)
            
        

        !-------------------------------------------------------------------------------
            !主要用来统计远场类型的边界
        !-------------------------------------------------------------------------------

        !-------------------------------------------------------------------------------
            if(BC_Temp(1).eq.2)then
               
       
                ii=ii+1
               
                do k=1,Z_Max
                do j=BC_Temp(5),BC_Temp(6)
                do i=BC_Temp(3),BC_Temp(4)
                    Wall_Points_Number(ii,1)=BC_Temp(4)
                    Wall_Points_Number(ii,2)=BC_Temp(6)    
                    Wall_Points_Number(ii,3)=BC_Temp(8)                      
                Wall_Points_Position(ii,i,j,k,1)=Points_Position_initial(BC_Temp(2),i,j,k,1)
                Wall_Points_Position(ii,i,j,k,2)=Points_Position_initial(BC_Temp(2),i,j,k,2)
                Wall_Points_Position(ii,i,j,k,3)=Points_Position_initial(BC_Temp(2),i,j,k,3)
                
               
                if (Wall_Points_Position(ii,i,j,k,1).le.X_Min)then
                    X_Min=Wall_Points_Position(ii,i,j,k,1)
                end if
                if (Wall_Points_Position(ii,i,j,k,2).le.Y_Min)then
                    Y_Min=Wall_Points_Position(ii,i,j,k,2)
                end if
              
                end do
                end do
                end do
                !write(*,*) BC_Temp(1),BC_Temp(2),BC_Temp(3),BC_Temp(4),BC_Temp(5),BC_Temp(6),BC_Temp(7),BC_Temp(8)
                 
       end if
       ! !-------------------------------------------------------------------------------
            if(BC_Temp(2).le.Num_Block)then
                go to 13
            end if
        !-------------------------------------------------------------------------------
 
14      continue     

     close(1)   
    call cpu_time(time_end)
    Read_Boundary_Condition_Time=time_end-time_begin
    end 

!--------------------------------------------------------------------------------------    
    
    
                