    program read_point_position
    
    use Define_Variables
    allocate(BC_Temp(8))
     
     call CPU_TIME(time_begin)
     call Count_Max
     
     open(10,file="BC_Type.dat")
     read(10,*)Num_Block
     write(*,"(a25,i5)")"The num of Block is :",Num_Block
     call Cal_WallNumber
     !��������ά��
     
     allocate(Points_Position_initial(Num_Block,X_max,Y_max,Z_max,3))
     allocate(Wall_Points_Position(Num_Block,X_max,Y_max,Z_max,3))
     allocate(num(Num_Block,3)) 
     allocate(Area(Num_Block,X_max,Y_max,Z_max),quality(Num_Block,X_max,Y_max,Z_max),area_1(Num_Block,X_max,Y_max,Z_max))
     allocate(Edge_distance(Num_Block,X_max,Y_max,Z_max,2))
     allocate(coefficient(Num_Block,X_max,Y_max,Z_max,3))
     allocate(Wall_Points_Number(4,3))
     allocate(Wall_Points_Rotate(4,X_max,Y_max,Z_max,3))
     allocate(Points_Position_RBF(Num_Block*X_max*Y_max*Z_max,3))
     allocate(kpoint_RBF(Num_Block*X_max*Y_max*Z_max,3))
     allocate(BC_Type(Edge_Count,3))
     !ά���������
      
     
     
     do ii=1,Num_Block
            read(10,*)num(ii,1),num(ii,2),num(ii,3)
            KMAX=num(ii,3)
            JMAX=num(ii,2)
            IMAX=num(ii,1)
            write(*,"(a5,i4,a3,a5,i4,a3,a5,i4)")"IMAX=",IMAX,"   ","JMAX=",JMAX,"   ","KAMX=",KMAX
            read(10,*)(((Points_Position_initial(ii,i,j,k,1),i=1,IMAX),j=1,JMAX),k=1,KMAX)
            read(10,*)(((Points_Position_initial(ii,i,j,k,2),i=1,IMAX),j=1,JMAX),k=1,KMAX)
            read(10,*)(((Points_Position_initial(ii,i,j,k,3),i=1,IMAX),j=1,JMAX),k=1,KMAX) 
            !write(*,*)(((Points_Position_initial(ii,i,j,k,1),i=1,IMAX),j=1,JMAX),k=1,KMAX)
        end do
        call geometry
        call Grid_Quality
        call TFI_Parameter
        call Read_Boundary_Condition
        call rotate(0.0,0.0,10.0)
        call RBF    
        
    
    
    call cpu_time(time_end)
    
    !call Geometry
    !�������ʱcpu���е�ʱ��
       write(*,"(a27,f10.7,a1)")"The Program used time is :",time_end-time_begin,"s"
       write(*,"(a27,f10.7,a1)")"The Read_Boundary_Condition used time is :",Read_Boundary_Condition_Time,"s"
       write(*,"(a27,f10.7,a1)")"The Count_Max used time is :",Count_Max_Time,"s"
       write(*,"(a27,f10.7,a1)")"The Cal_WallNumber_Time used time is :",Cal_WallNumber_Time,"s"
       write(*,"(a27,f10.7,a1)")"The Grid_Quality_Time used time is :",Cal_WallNumber_Time,"s"
       
       
       !����������е�ʱ��
write(*,*)"I want to change world"
    end
