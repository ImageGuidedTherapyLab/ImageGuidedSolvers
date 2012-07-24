%  Matlab routine that analyzes CSI DATA PixelWise
function MapPixel=MapCSI(RealPixel,ImagPixel,EchoSpacing,Larmor);
data = RealPixel + i * ImagPixel;
MapPixel=zeros(6,1);
if abs(data)>50;   % Only processes data with signal
  % Timing-correction for first echo
  data=data.*exp(-2*pi*i*31.25e3*4e-6);  
  [p q]=stmcb(data,1,2);
  rtz=roots(q);
  ppm=imag(log(rtz))/2/pi/(1e-3*EchoSpacing*Larmor);
  MapPixel(1)=ppm(1);
  MapPixel(2)=ppm(2);
  t2star=-1./real(log(rtz)/EchoSpacing)/1000;
  MapPixel(3)=((p(1).*rtz(1)+p(2))/((2.*q(1).*rtz(1))+q(2)));
  MapPixel(4)=((p(1).*rtz(2)+p(2))/((2.*q(1).*rtz(2))+q(2)));
  MapPixel(5)=(t2star(1));
  MapPixel(6)=(t2star(2));
end
