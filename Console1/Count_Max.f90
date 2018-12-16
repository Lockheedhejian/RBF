subroutine Count_Max
        use Define_Variables
        implicit double precision(a-h,o-z)
        
        call cpu_time(time_begin)
        
        open(11,file="BC_Type.x")
        read(11,*)number_region
        !write(*,*)number_region
        allocate(number(number_region,3))
        temp_num=0
        do ii=1,number_region
            read(11,*)number(ii,1),number(ii,2),number(ii,3)!读取.x中每一行的数据
                                                    !在导出的数据中，会输出分区边界的网格点数，以及三个方向上的节点数
            if(number(ii,2).eq.1.and.number(ii,3).eq.1)then
                temp_num=temp_num+1
                go to 11
            end if
11  continue
            if(number(ii,2).ne.1)then
                X_max=max(X_max,number(ii,1))
                Y_max=max(Y_max,number(ii,2))
                Z_max=max(Z_max,number(ii,3))
                C_MAX=max(X_max,Y_max,Z_max)
            end if    
           
        end do 
        call cpu_time(time_end)
        !-------------------------------------------------------------------------------
        Count_Max_Time=time_end-time_begin
        
        
         
    end