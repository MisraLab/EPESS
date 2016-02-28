function [phi] = simulate(slice_range) 

% Given the angle slice range, this code samples a point unformly from it.

 diff = zeros(length(slice_range),1);
     
   for i=1:2:length(slice_range)
         
         diff(i,1) = slice_range(i+1) - slice_range(i);
         
   end

  total = sum(diff);
  diff = diff/total;
  index = discretesample(diff,1);
  phi = rand*(slice_range(index+1) - slice_range(index)) + slice_range(index);
  
  
end