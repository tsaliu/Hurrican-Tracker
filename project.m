clc; clear; close all;
warning off;
%load raw;
%[num,txt,raw] = xlsread('hurricane harvy 2017.xlsx');
% raw = textread('data.txt','%f');
% data=[];
% data=reshape(raw,[],121)';
% save data;
load data;
%create the map axes. See http://www.mathworks.com/help/map/the-map-frame.html
axesm('MapProjection','robinson',...
     'MapLatLimit',[10 40],'MapLonLimit',[-100 -50], ...
    'Frame','on','Grid','on', 'MeridianLabel', 'on', 'ParallelLabel', 'on')

%load the coast data and plot. See http://www.mathworks.com/help/map/create-a-world-map.html
figure(1)
load coast
[latcells, loncells] = polysplit(lat, long);
numel(latcells);
axis off;
tightmap;


plotm(lat, long);   %plot land
date=data(:,1);
time=data(:,2);
latd=data(:,3);
lond=data(:,4);
wind=data(:,5);
press=data(:,6);


hold on
% pcolorm(latd,lond,val);
h=plotm(latd,lond,wind,'-xr');
xmap=h.XData;   %lon
ymap=h.YData;   %lat
% 
% xmap(1:10)=[];
% ymap(1:10)=[];
% time(1:10)=[];
% date(1:10)=[];

%plot(xmap,ymap,'o');


%calc mapping ------ nonlinear (this is linear, not working)
% syms a b;
% eq1 = a*lond(1)+b == xmap(1);
% eq2 = a*lond(112)+b ==  xmap(112);
% [A,B]=equationsToMatrix([eq1, eq2],[a,b]);
% sol1=linsolve(A,B);
% 
% syms a2 b2;
% eq11 = a2*latd(1)+b2 == ymap(1);
% eq22 = a2*latd(112)+b2 == ymap(112);
% [A2,B2]=equationsToMatrix([eq11,eq22],[a2,b2]);
% sol2=linsolve(A2,B2);
% mapping(1,1)=1*sol1(1);
% mapping(1,2)=1*sol1(2);
% mapping(2,1)=1*sol2(1);
% mapping(2,2)=1*sol2(2);
% save mapping;
% load mapping;
% 
% % x=lond*0.0143859+1.08506;
% % y=latd*0.01676+0.000144;
% 
% x=lond*sol1(1)+sol1(2);
% y=latd*sol2(1)+sol2(2);
%   
% plot(x,y,'o');



[x,dx,y,dy,T]=xymodel(date,time,xmap,ymap);
attime=cumsum(T);
orgattime=attime;
attime=[0 attime];
figure(2)
subplot(2,1,1)
plot(attime,x);
hold on
plot(attime,y),grid;
ylabel('position (m)');xlabel('time (h)')
legend('x','y');
title('x,y');


for kk=1:length(attime)-1
    for k=attime(kk):attime(kk+1)
        k=fix(k);
        attimec(k+1)=k;
        dxc(k+1)=dx(kk);
        dyc(k+1)=dy(kk);
    end
end
subplot(2,1,2)
% plot(attime,dx);
% hold on
% plot(attime,dy);
% plot(attime,dz),grid;
plot(attimec,dxc);
hold on
plot(attimec,dyc),grid;
ylabel('velocity (m/h)');xlabel('time (h)');
legend('x','y');
title('velocity');

T=round(T);
P00=zeros([4 4]);
%P00=(1e6/3)^2*eye(6);
% P00=1e6*eye(6);
covp=(1e-5/3)^2;
covv=5e-4;
P00(1,1)=covp;
P00(3,3)=covp;
P00(2,2)=covv;
P00(4,4)=covv;

winsize=10;
% [hxk1k1,Pk1k1,nis,nees,crlb]=kalman(P00,T,x,y,dx,dy);
% [hxk1k1,Pk1k1,nis,nees,crlb]=slidekalman(P00,T,x,y,dx,dy,winsize);
[hxk1k1,Pk1k1,nis,nees,crlb,eak,prob]=fadekalman(P00,T,x,y,dx,dy);
figure(3)
subplot(2,1,1)
plot(attime,eak),grid;
title('normalized innovation square')
subplot(2,1,2)
plot(attime,prob),grid;
title('tail probability')

figure(4)
subplot(2,2,1)
plot(orgattime,hxk1k1(1,(1:4:length(hxk1k1)))), grid;
hold on
plot(orgattime,hxk1k1(3,(3:4:length(hxk1k1))));
title('pos');

subplot(2,2,2)
plot(orgattime,hxk1k1(2,(2:4:length(hxk1k1)))), grid;
hold on
plot(orgattime,hxk1k1(4,(4:4:length(hxk1k1))));
title('vel');

subplot(2,2,3)
plot(orgattime,Pk1k1(1,(1:4:length(Pk1k1)))),grid;
hold on
plot(orgattime,Pk1k1(3,(3:4:length(Pk1k1))));
title('cov pos');

subplot(2,2,4)
plot(orgattime,Pk1k1(2,(2:4:length(Pk1k1)))),grid;
hold on
plot(orgattime,Pk1k1(4,(4:4:length(Pk1k1))));
title('cov vel');

% [dlat,dlon]=covlatlon(hxk1k1(1,(1:6:length(hxk1k1))),hxk1k1(3,(3:6:length(hxk1k1))),hxk1k1(5,(5:6:length(hxk1k1))));
% xp=hxk1k1(1,(1:6:length(hxk1k1)))+Pk1k1(1,(1:6:length(Pk1k1)));
% yp=hxk1k1(3,(2:6:length(hxk1k1)))+Pk1k1(3,(1:6:length(Pk1k1)));
% zp=hxk1k1(5,(1:6:length(hxk1k1)))+Pk1k1(5,(1:6:length(Pk1k1)));
% yp=hxk1k1(3,(3:6:length(hxk1k1)));
% zp=hxk1k1(5,(5:6:length(hxk1k1)));
% [dlatp,dlonp]=covlatlon(xp,yp,zp);
% covrad=real(sqrt((dlatp-dlat).^2+(dlonp-dlon).^2))./1000;

tx=hxk1k1(1,(1:4:length(hxk1k1)));
ty=hxk1k1(3,(3:4:length(hxk1k1)));
figure(1)
h2=plot(tx,ty,'-xk');

% h2=plotm(dlat,dlon,'-xk');
% legend([h,h2],'Actual','Estimate');
% hold on
% drawcov(h2.XData,h2.YData,covrad);


figure(5)
subplot(2,1,1)
plot(orgattime,nis(1,(1:4:length(nis)))),grid;
hold on
plot(orgattime,nis(3,(3:4:length(nis))));
title('NIS pos');

subplot(2,1,2)
plot(orgattime,nis(2,(2:4:length(nis)))),grid;
hold on
plot(orgattime,nis(4,(4:4:length(nis))));
title('NIS vel');

figure(6)
subplot(2,1,1)
plot(orgattime,nees(1,(1:4:length(nees)))),grid;
hold on
plot(orgattime,nees(3,(3:4:length(nees))));
title('NEES pos');

subplot(2,1,2)
plot(orgattime,nees(2,(2:4:length(nees)))),grid;
hold on
plot(orgattime,nees(4,(4:4:length(nees))));
title('NEES vel');

% xb=hxk1k1(1,(1:6:length(hxk1k1)));
% yb=hxk1k1(3,(3:6:length(hxk1k1)));
% zb=hxk1k1(5,(5:6:length(hxk1k1)));
% dxb=hxk1k1(2,(2:6:length(hxk1k1)));
% dyb=hxk1k1(4,(4:6:length(hxk1k1)));
% dzb=hxk1k1(6,(6:6:length(hxk1k1)));
% xbb=mean(xb);
% ybb=mean(yb);
% zbb=mean(zb);
% dxbb=mean(dxb);
% dybb=mean(dyb);
% dzbb=mean(dzb);
% rmsex=sqrt((1/length(xb))*sum((xb-xbb).^2));

figure(1)

Nruns=10000;
[rnis,rnees,rmse,hxk1k1b,Pk1k1b,svar,varsmean,varsvar,reak,rprob]=slideruns(Nruns,P00,T,x,y,dx,dy);
figure(7)
subplot(2,1,1)
plot(orgattime,rnis(1,(1:4:length(rnis)))),grid;
hold on
plot(orgattime,rnis(3,(3:4:length(rnis))));
title([num2str(Nruns) ' Runs NIS pos']);

subplot(2,1,2)
plot(orgattime,rnis(2,(2:4:length(rnis)))),grid;
hold on
plot(orgattime,rnis(4,(4:4:length(rnis))));
title([num2str(Nruns) ' Runs NIS vel']);
    
figure(8)
subplot(2,1,1)
plot(orgattime,rnees(1,(1:4:length(rnees)))),grid;
hold on
plot(orgattime,rnees(3,(3:4:length(rnees))));
title([num2str(Nruns) ' Runs NEES pos']);

subplot(2,1,2)
plot(orgattime,rnees(2,(2:4:length(rnees)))),grid;
hold on
plot(orgattime,rnees(4,(4:4:length(rnees))));
title([num2str(Nruns) ' Runs NEES vel']);
figure(9)
subplot(2,1,1)
plot(orgattime,rmse(1,:)),grid;
hold on
plot(orgattime,rmse(2,:));
title(['RMSE of ' num2str(Nruns) ' Runs pos']);

subplot(2,1,2)
plot(orgattime,rmse(3,:)),grid;
hold on
plot(orgattime,rmse(4,:));
title(['RMSE of ' num2str(Nruns) ' Runs vel']);

figure(10)
subplot(2,2,1)
plot(orgattime,hxk1k1b(1,:));
hold on
plot(orgattime,hxk1k1b(2,:)),grid;
title([num2str(Nruns) ' Runs pos']);
xlabel('time (h)')
legend('x','y');

subplot(2,2,2)
plot(orgattime,hxk1k1b(3,:));
hold on
plot(orgattime,hxk1k1b(4,:)),grid;
title([num2str(Nruns) ' Runs vel']);
xlabel('time (h)');
legend('x','y');

subplot(2,2,3)
plot(orgattime,Pk1k1b(1,:));
hold on
plot(orgattime,Pk1k1b(3,:)),grid;
title([num2str(Nruns) ' Runs pos cov']);
xlabel('time (h)')
legend('x','y');

subplot(2,2,4)
plot(orgattime,Pk1k1b(2,:));
hold on
plot(orgattime,Pk1k1b(4,:)),grid;
title([num2str(Nruns) ' Runs vel cov']);
xlabel('time (h)')
legend('x','y');

figure(14)
subplot(2,1,1)
plot(attime,reak),grid;
title([num2str(Nruns) ' Runs normalized innovation square'])
subplot(2,1,2)
plot(attime,rprob),grid;
title([num2str(Nruns) ' Runs tail probability'])


% [latb,lonb]=covlatlon(hxk1k1b(1,:),hxk1k1b(2,:),hxk1k1b(3,:));
rxb=hxk1k1b(1,:);
ryb=hxk1k1b(2,:);
figure(1)
h3=plot(rxb,ryb,'-xb');
% legend([h,h2,h3],'Actual','Estimate','N Runs Estimate');
legend([h,h2,h3],'Actual','Single Estimate',[num2str(Nruns) ' Estimate']);
% drawcov(rxb,ryb,Pk1k1b(1,:),Pk1k1b(3,:));

figure(11)
subplot(2,1,1)
plot(orgattime,svar(1,:)),grid;
hold on
plot(orgattime,svar(2,:));
title(['sample var of ' num2str(Nruns) ' Estimation pos']);

subplot(2,1,2)
plot(orgattime,svar(3,:)),grid;
hold on
plot(orgattime,svar(4,:));
title(['sample var of ' num2str(Nruns) ' Estimation vel']);


figure(12)
subplot(2,2,1)
plot(orgattime,varsmean(1,:)),grid;
hold on
plot(orgattime,varsmean(2,:));
title(['var of sample mean for ' num2str(Nruns) ' Estimation pos']);

subplot(2,2,2)
plot(orgattime,varsmean(3,:)),grid;
hold on
plot(orgattime,varsmean(4,:));
title(['var of sample mean for ' num2str(Nruns) ' Estimation vel']);


subplot(2,2,3)
plot(orgattime,varsvar(1,:)),grid;
hold on
plot(orgattime,varsvar(2,:));
title(['var of sample var for ' num2str(Nruns) ' Estimation pos']);

subplot(2,2,4)
plot(orgattime,varsvar(3,:)),grid;
hold on
plot(orgattime,varsvar(4,:));
title(['var of sample var for ' num2str(Nruns) ' Estimation vel']);

numest=500;
numrun=100;
[meannest,varnest]=sliderunsruns(numest,P00,T,x,y,dx,dy,numrun);
figure(13)
subplot(4,2,1)
plot(attime,x),grid;
hold on
plot(orgattime,meannest(1,:));
legend('meas',['mean of ' num2str(numest) ' estimate for ' num2str(numrun) ' runs']);
title(['mean of ' num2str(numest) ' estimate for ' num2str(numrun) ' runs for pos x']);

subplot(4,2,3)
plot(attime,y),grid;
hold on 
plot(orgattime, meannest(2,:));
legend('meas',['mean of ' num2str(numest) ' estimate for ' num2str(numrun) ' runs']);
title(['mean of ' num2str(numest) ' estimate for ' num2str(numrun) ' runs for pos y']);

subplot(4,2,5)
plot(attimec,dxc),grid;
hold on 
plot(orgattime, meannest(3,:));
legend('meas',['mean of ' num2str(numest) ' estimate for ' num2str(numrun) ' runs']);
title(['mean of ' num2str(numest) ' estimate for ' num2str(numrun) ' runs for vel x']);


subplot(4,2,7)
plot(attimec,dyc),grid;
hold on 
plot(orgattime, meannest(4,:));
legend('meas',['mean of ' num2str(numest) ' estimate for ' num2str(numrun) ' runs']);
title(['mean of ' num2str(numest) ' estimate for ' num2str(numrun) ' runs for vel y']);

subplot(4,2,2)
plot(orgattime, varnest(1,:)),grid;
title(['var of ' num2str(numest) ' estimate for ' num2str(numrun) ' runs for pos x']);

subplot(4,2,4)
plot(orgattime, varnest(2,:)),grid;
title(['var of ' num2str(numest) ' estimate for ' num2str(numrun) ' runs for pos y']);

subplot(4,2,6)
plot(orgattime, varnest(3,:)),grid;
title(['var of ' num2str(numest) ' estimate for ' num2str(numrun) ' runs for vel x']);

subplot(4,2,8)
plot(orgattime, varnest(1,:)),grid;
title(['var of ' num2str(numest) ' estimate for ' num2str(numrun) ' runs for vel y']);
