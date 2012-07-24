clear all
close all


syms x_lb x_ub x gamma


phi = exp(gamma*(x-x_ub)/(x_ub-x_lb))+exp(-gamma*(x-x_lb)/(x_ub-x_lb));

dphidx = simple(diff(phi,x));

pretty(phi);
pretty(dphidx);
%xmin= .05; xmax= .15;
%xmid= (xmin+xmax)/2;
%lenx=xmax-xmin;
%
%x=[xmin-lenx/10:lenx/100:xmax+lenx/10];
%
%y1=exp(1000/lenx*(x-xmax));
%y2=exp(-1000/lenx*(x-xmin)); 
%y=y1+y2;
%
%
%plot(x,y)
%hold
%plot(x,y1,'r')
%plot(x,y2,'g')

