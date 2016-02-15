
function [check] = belongs(phi, angle_slice) 

% This function is just a check to see if a proposed point belongs to some
% angle range

         in_slice = zeros(length(angle_slice),1);


         for i=1:2:length(angle_slice)

             if phi>= angle_slice(i) && phi<= angle_slice(i+1)
                in_slice(i,1) = 1;

             else 

                in_slice(i,1) = 0;
             end

         end

         check = any(in_slice);
end
