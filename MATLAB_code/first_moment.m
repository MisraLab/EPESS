function [value] = first_moment(xx, nu, slice_range)

value = 0;
for i=1:2:length(slice_range)
    k1 = sin(i+1) - sin(i);
    k2 = -(cos(i+1) - cos(i));
    value = value + xx*k1 + nu*k2;
end

end