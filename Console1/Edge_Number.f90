subroutine Cal_WallNumber
        use Define_Variables
        implicit double precision(A_Matrix-h,o-z)
        character::temp*30
        
        
        call cpu_time(time_begin)
        open(12,file="BC_Type.fvbnd")
    
        read(12,*)
        

        do ii=1,11
            read(12,*)temp
            if(temp=="BOUNDARIES")then 
            go to 111 
            end if      
        end do
111      continue    
!       
        Edge_Count=0
113       continue 
        read(12,*,end=114)BC_Temp(1),BC_Temp(2),BC_Temp(3),BC_Temp(4),BC_Temp(5),BC_Temp(6),BC_Temp(7),BC_Temp(8)

            if(BC_Temp(1).eq.2)then
               Edge_Count=Edge_Count+1!用来记录网格中边界类型为Wall的数量
            end if
            if(BC_Temp(2).lt.Num_Block)then
                go to 113
            end if            
            
114         continue   
            call cpu_time(time_end)
            Cal_WallNumber_Time=time_end-time_begin
           
    end