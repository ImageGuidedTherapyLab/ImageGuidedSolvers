
clear all 
close all

dirid = 1145831;
exampath = './s%d/i%d.MRDC.%d'
ntime = 30;
dirid = 12345;
exampath = './knownphase/s%d/i%d.MRDC.%d'
ntime = 10;
nslice = 1;
snrthresh=5;
alpha = 1/2/pi;
te=1/(1e-3*63.868756);
tmap      =zeros(4, 4, ntime, nslice);
phasediff =zeros(4, 4, ntime, nslice);
dat       =zeros(4, 4, ntime, nslice);

for itime = 1:ntime
  for kk = 1:nslice
    realid = 2* nslice *(itime-1) + (2 *(kk-1) + 1)       ; 
    imagid = 2* nslice *(itime-1) + (2 *(kk-1) + 1) + 1   ; 
    disp(sprintf('itime = %d nslice = %d  realid = %d imagid = %d', itime,kk,realid,imagid));
    timg(:,:,1)=dicomread( sprintf(exampath,dirid, dirid + realid ,realid ) );
    timg(:,:,2)=dicomread( sprintf(exampath,dirid, dirid + imagid ,imagid ) );
    dat(:,:,itime,kk) = cplximg(double(timg));

    if itime>1,
    phasediff(:,:,itime,kk)=angle(dat(:,:,itime-1,kk).*conj(dat(:,:,itime,kk)));
    tmap(:,:,itime,kk) = (tmap(:,:,itime-1,kk)+phitmap(phasediff(:,:,itime,kk),te,alpha));
    end
  end
end
imagesc(tmap(:,:,7,1))
colormap 'hot'
caxis([0 20])
