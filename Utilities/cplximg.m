function cimg=cplximg(data)
%CIMG combines real and imaginary parts into a complex number. 
%
%CIMAGE = CPLXIMG(DATA,DATA2)
%
%CPLXIMG takes: 
%A single 3D variable as its argument (DATA).
%DATA is assumed to be interleaved REAL and IMAGINARY
%images (respectively, i.e., real1, imag1, real1, imag2 ...)
%The data is combined to form complex images.
%
%
%Author: R. Jason Stafford
%Date: 06/98
%

list = [0:size(data,3)-1];
masqr = ~mod(list,2);
masqi = ~mod(list-1,2);
cimg = data(:,:,masqr) + i*data(:,:,masqi);
