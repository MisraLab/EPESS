function mvndensity  = mnvdensity(x, mu, scalar)

mvndensity = exp(- (x-mu) * (x-mu)' / (2*scalar) ) / (sqrt(2 * pi * abs(scalar) ));

end