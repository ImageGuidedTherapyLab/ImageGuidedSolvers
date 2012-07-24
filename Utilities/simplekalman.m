% kalman filter for relative temperature difference
clear all
close all

% imaging parameters
dirid = 1145831;
exampath = './s%d/i%d.MRDC.%d'
exampath = '/home/dfuentes/DDDAS/data/mdacc/canine_oct08/s%d/i%d.MRDC.%d'
ntime = 30;
nslice = 5;
alpha = -0.0097;
te=15;
x0 = 1;
xf = 256;
y0 = 1;
yf = 256;
%x0 = 80;
%xf = 130;
%y0 = 120;
%yf = 170;
xpixel = xf-x0+1;
ypixel = yf-y0+1;
tmap =zeros(xpixel, ypixel, nslice, ntime);
dat  =zeros(xpixel, ypixel, nslice, ntime);

% Define the temperature to measure the temperature itself:
s.H = speye(xpixel*ypixel);
% Define a measurement error (stdev) of .5 degC:
s.R = .5^2 * speye(xpixel*ypixel); % variance, hence stdev^2

% initial conditions
s.x = zeros(xpixel*ypixel,1); 
%s.P = inv(s.H)*s.R*inv(s.H'); 
s.P = s.R; 
% build FD state transition matrix and load
s.B = speye(xpixel*ypixel);
s.u = zeros(xpixel*ypixel,1);
s.A = speye(xpixel*ypixel);
% Define a process noise 
s.Q = .5^2 * speye(xpixel*ypixel); %(stdev) of .5 degC, hence stdev^2
s.Q = 0; % no noise in FD prediction

for itime = 1:ntime
  for kk = 1:nslice
    realid = 2* nslice *(itime-1) + (2 *(kk-1) + 1)       ; 
    imagid = 2* nslice *(itime-1) + (2 *(kk-1) + 1) + 1   ; 
    disp(sprintf('itime = %d nslice = %d  realid = %d imagid = %d', itime,kk,realid,imagid));
    timg(:,:,1)=dicomread( sprintf(exampath,dirid, dirid + realid ,realid ) );
    timg(:,:,2)=dicomread( sprintf(exampath,dirid, dirid + imagid ,imagid ) );
    dat(:,:,kk,itime) = cplximg(double(timg(x0:xf,y0:yf,:)));
    % create tmap
    if itime>1,
      phasediff(:,:,kk,itime) = angle(dat(:,:,kk,itime-1).*conj(dat(:,:,kk,itime)));
      tmap(:,:,kk,itime) = (tmap(:,:,kk,itime-1)+phitmap(phasediff(:,:,kk,itime),te,alpha));
      %s(end).z = reshape(tmap(:,:,kk,itime),[xpixel*ypixel, 1]);
      %s(end+1)=kalmanf(s(end)); % perform a Kalman filter iteration
    end
  end
end

imagesc(tmap(:,:,2,27))
colormap 'hot'
caxis([0 20])

