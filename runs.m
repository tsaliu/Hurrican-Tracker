function [rnis,rnees,rmse,hxk1k1b,Pk1k1b,svar,varsmean,varsvar]=runs(nruns,P00,T,x,y,dx,dy,N)
nis=[];
nees=[];
rnis=[];
rnees=[];
xbt=[];
ybt=[];

dxbt=[];
dybt=[];

pxt=[];
pyt=[];
pdxt=[];
pdyt=[];

xbb=[];
ybb=[];

dxbb=[];
dybb=[];

rmse=[];
hxk1k1b=[];
Pk1k1b=[];
% crlbb=[];

svar=[];
varsmean=[];
varsvar=[];
for i=1:nruns
%     [hxk1k1,Pk1k1,nis,nees,crlb]=kalman(P00,T,x,y,dx,dy);
    [hxk1k1,Pk1k1,nis,nees,crlb]=slidekalman(P00,T,x,y,dx,dy,N);
    nis=nis+nis;
    nees=nees+nees;
    
    xb=hxk1k1(1,(1:4:length(hxk1k1)));
    yb=hxk1k1(3,(3:4:length(hxk1k1)));

    dxb=hxk1k1(2,(2:4:length(hxk1k1)));
    dyb=hxk1k1(4,(4:4:length(hxk1k1)));

    xbt=[xbt; xb];
    ybt=[ybt; yb];

    dxbt=[dxbt; dxb];
    dybt=[dybt; dyb];

    px=Pk1k1(1,(1:4:length(hxk1k1)));
    pdx=Pk1k1(2,(2:4:length(hxk1k1)));
    py=Pk1k1(3,(3:4:length(hxk1k1)));
    pdy=Pk1k1(4,(4:4:length(hxk1k1)));
    
    pxt=[pxt; px];
    pyt=[pyt; py];
    pdxt=[pdxt; pdx];
    pdyt=[pdyt; pdy];
end
    rnis=nis./i;
    rnees=nees./i;
    

    for ii=1:length(xb)
        xbb=mean(xbt(:,ii));
        ybb=mean(ybt(:,ii));

        dxbb=mean(dxbt(:,ii));
        dybb=mean(dybt(:,ii));
        
        pxb=mean(pxt(:,ii));
        pdxb=mean(pdxt(:,ii));
        pyb=mean(pyt(:,ii));
        pdyb=mean(pdyt(:,ii));
 
        
        rmse(1,ii)=sqrt((1/i).*sum((xbt(:,ii)-xbb).^2));
        rmse(2,ii)=sqrt((1/i).*sum((ybt(:,ii)-ybb).^2));
 
        rmse(3,ii)=sqrt((1/i).*sum((dxbt(:,ii)-dxbb).^2));
        rmse(4,ii)=sqrt((1/i).*sum((dybt(:,ii)-dybb).^2));

        
        hxk1k1b(1,ii)=xbb;
        hxk1k1b(2,ii)=ybb;

        hxk1k1b(3,ii)=dxbb;
        hxk1k1b(4,ii)=dybb;
        
        Pk1k1b(1,ii)=pxb;
        Pk1k1b(2,ii)=pdxb;
        Pk1k1b(3,ii)=pyb;
        Pk1k1b(4,ii)=pdyb;
    
        svar(1,ii)=(1/(i-1))*sum((xbt(:,ii)-xbb).^2);
        svar(2,ii)=(1/(i-1))*sum((ybt(:,ii)-ybb).^2);
        svar(3,ii)=(1/(i-1))*sum((dxbt(:,ii)-dxbb).^2);
        svar(4,ii)=(1/(i-1))*sum((dybt(:,ii)-dybb).^2);
        
        varsmean(1,ii)=mean((xbt(:,ii)-x(ii)).^2);
        varsmean(2,ii)=mean((ybt(:,ii)-y(ii)).^2);
        varsmean(3,ii)=mean((dxbt(:,ii)-dx(ii)).^2);
        varsmean(4,ii)=mean((dybt(:,ii)-dy(ii)).^2);
        
        varsvar(1,ii)=mean((pxt(:,ii)-pxb).^2);
        varsvar(2,ii)=mean((pyt(:,ii)-pyb).^2);
        varsvar(3,ii)=mean((pdxt(:,ii)-pdxb).^2);
        varsvar(4,ii)=mean((pdyt(:,ii)-pdyb).^2);
    end
 

    

end