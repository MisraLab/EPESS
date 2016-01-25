
    function [check] = belongs(phi, angle_slice) 
             in_slice = zeros(length(angle_slice));
    
             for i=1:2:length(angle_slice)
                 if phi>= angle_slice(i) && phi<= angle_slice(i+1)
                    in_slice(i) = 1;
            
                 else 
                
                    in_slice(i) = 0;
                 end
            
             end
        
             check = any(inslice);
     end
