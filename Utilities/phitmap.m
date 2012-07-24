function y=phitmap(image,te,coeff,b0)
%  PHITMAP(IMAGE,TEMPERATURE SENSITIVITY COEFFICIENT (ppm/deg. Celsius), TE(ms))
%  
%  IMAGE = A phase difference map of arbitrary size.
%  TEMP. COEFF. = Sensitivity of PRF to temperature (ppm/deg. C)
%  TE = Echo time in ms
%
%  This function assumes that the input 'image' is a phase difference MR image
%  and creates a thermal map based on the phase difference, time to echo (TE)
%  and the temperature sensitivity coefficient.
%
%  The temperature sensitivity coefficient of water is ~ .01 ppm/degC.  The 
%  established range is from about .007-.011 ppm/degC for a variety of tissue.
%  Often used is .009
%
%

%
% Author: R. Jason Stafford
% Date: 3/98
% Revision:
%
if nargin==2
   y = image./(2*3.14159*.01*te*1e-3*63.868756);
end


if nargin==4,
   for kk=1:size(image,3)
      y(:,:,kk) = image(:,:,kk)./(2*3.14159*coeff*te*1e-3*b0);
   end
else
   y = image./(2*3.14159*coeff*te*1e-3*63.868756);
end
