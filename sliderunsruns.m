function [meannest,varnest]=sliderunsruns(nruns,P00,T,x,y,dx,dy,M)

xxt=[];
yyt=[];
dxxt=[];
dyyt=[];

meannest=[];
varnest=[];

for i=1:M
    [rnis,rnees,rmse,hxk1k1b,Pk1k1b,svar,varsmean,varsvar,eak,prob]=slideruns(nruns,P00,T,x,y,dx,dy);
    
    xx=hxk1k1b(1,:);
    yy=hxk1k1b(2,:);
    
    dxx=hxk1k1b(3,:);
    dyy=hxk1k1b(4,:);
    
    xxt=[xxt; xx];
    yyt=[yyt; yy];
    dxxt=[dxxt; dxx];
    dyyt=[dyyt; dyy];
    
end

for ii=1:length(xx)
   xxb=mean(xxt(:,ii));
   yyb=mean(yyt(:,ii));
   
   dxxb=mean(dxxt(:,ii));
   dyyb=mean(dyyt(:,ii));
   
   meannest(1,ii)=xxb;
   meannest(2,ii)=yyb;
   meannest(3,ii)=dxxb;
   meannest(4,ii)=dyyb;
   
   varnest(1,ii)=mean((xxt(:,ii)-x(ii)).^2);
   varnest(2,ii)=mean((yyt(:,ii)-y(ii)).^2);
   varnest(3,ii)=mean((dxxt(:,ii)-dx(ii)).^2);
   varnest(4,ii)=mean((dyyt(:,ii)-dy(ii)).^2);
   
end

end