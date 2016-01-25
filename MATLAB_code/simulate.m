function [phi] = simulate(slice_range) 

% Given the angle slice range, this code samples a point unformly from it.

 diff = zeros(length(slice_range),1);
%  in_slice_min = zeros(length(slice_range),1);
%  in_slice_max = zeros(length(slice_range),1);
% 
%  
%      for i=1:2:length(slice_range)
%                   
%          
%          if phi_min>= slice_range(i) && phi_min<= slice_range(i+1)
%                     in_slice_min(i,1) = 1;
%          end
%          
%      end
%      
%      for i=1:2:length(slice_range)
%                   
%          if phi_max>= slice_range(i) && phi_max<= slice_range(i+1)
%                     in_slice_max(i,1) = 1;
%          end
%               
%      end
%      
%      
%      
%      if any(in_slice_min) ==1
%          k1 = find(in_slice_min);
%          slice_range(k1) = phi_min;
%      end
%      
%      
%      if  any(in_slice_max) ==1
%          k2 = find(in_slice_max);
%          slice_range(k2+1) = phi_max;
%      end
     
   for i=1:2:length(slice_range)
         
         diff(i,1) = slice_range(i+1) - slice_range(i);
         
   end

  total = sum(diff);
  diff = diff/total;
  index = discretesample(diff,1);
  phi = rand*(slice_range(index+1) - slice_range(index)) + slice_range(index);
  
  
end