 # Read about the Van der Pol oscillator at
 #   http://en.wikipedia.org/wiki/Van_der_Pol_oscillator.

 function [vyd] = fvanderpol(vt, vy)
   #mu = varargin;
   vyd = [vy(2); 3 .* (1 - vy(1).^2) .* vy(2) - vy(1)];
 endfunction 
