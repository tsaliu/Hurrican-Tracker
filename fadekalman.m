function [xout,pout,nis,nees,crlb,eak,prob]=fadekalman(P00,T,x,y,dx,dy)
xout=[];
pout=[];
nis=[];
nees=[];
xtilt=zeros([4 4]);
crlb=[];
u=0;
sigv=sqrt(P00(1,1));
H=[1 0 1 0];%[1 0 1 0 1 0]
% x00=[sqrt(P00(1,1))*randn(1)+x(1) 0 0 0 0 0;
%     0 sqrt(P00(2,2))*randn(1)+dx(1) 0 0 0 0;
%     0 0 sqrt(P00(3,3))*randn(1)+y(1) 0 0 0;
%     0 0 0 sqrt(P00(4,4))*randn(1)+dy(1) 0 0
%     0 0 0 0 sqrt(P00(5,5))*randn(1)+z(1) 0
%     0 0 0 0 0 sqrt(P00(6,6))*randn(1)+dz(1)];
x00=[x(1) 0 0 0;
    0 dx(1) 0 0;
    0 0 y(1) 0;
    0 0 0 dy(1)];

%R=1*eye(6);
sigr=0.1;
R=sigr^2;

a=0.8;
ea=0;
eak=[];
eak(1)=ea;
nz=2;
na=nz*(1+a)/(1-a);
prob=[];
prob(1)=0;
upinterval=0.95;
lowinterval=0.05;

for k=1:length(T)
    if prob(k)<upinterval && prob(k)>lowinterval
        [F,G,Q]=FGfromTWNA(T(k),sigv,2);
    %     obs=[H;H*F;H*F^2;H*F^3]
    %     rrrr=rank(obs)
        if k == 1
    %        xx=x00+normrnd(0,sigv,length(x00),length(x00));
           xx=x00+mvnrnd(zeros([length(Q) length(Q)]),sqrt(Q));
           hxkk=xx;
           zk1=H*xx+normrnd(0,sigr,1,length(H));
           Pkk=P00;
        else
            xx=[x(k) 0 0 0;
            0 dx(k) 0 0;
            0 0 y(k) 0;
            0 0 0 dy(k)];
    %         xx=xx+normrnd(0,sigv,length(xx),length(xx));
            xx=xx+mvnrnd(zeros([length(Q) length(Q)]),sqrt(Q));
            zk1=H*xx+normrnd(0,sigr,1,length(H));
            hxkk=hxk1k1;
            Pkk=Pk1k1;
        end
    else 
            [F,G,Q]=FGfromTWNA(T(k),sigv,2);
            xx=[x(k) 0 0 0;
            0 dx(k) 0 0;
            0 0 y(k) 0;
            0 0 0 dy(k)];
    %         xx=xx+normrnd(0,sigv,length(xx),length(xx));
            xx=xx+mvnrnd(zeros([length(Q) length(Q)]),sqrt(Q));
            zk1=H*xx+normrnd(0,sigr,1,length(H));
            hxkk=xx;
            Pkk=P00;
    end
    

        hxk1k=F*hxkk;   %+G*u
        hzk1k=H*hxk1k;
        vk1=zk1-hzk1k;

        Pk1k=F*Pkk*F'+Q;
        Sk1=R+H*Pk1k*H';
        Wk1=Pk1k*H'*(Sk1^-1);
        Pk1k1=Pk1k-Wk1*Sk1*Wk1';
        hxk1k1=hxk1k+Wk1*vk1;

        nisk=vk1'*(Sk1^-1)*vk1;

        xtil=xx-hxkk;
        neesk=xtil'*(Pkk^-1)*xtil;

        xtilt=xtilt+xtil;

        crlbb=xtilt.^2./k;

        e=vk1'*(Sk1^-1)*vk1;
        ee=[e(1,1) e(2,2) e(3,3) e(4,4)];
        ek=max(ee);
        eak(k+1)=a.*eak(k)+ek;
        prob(k+1)=chi2cdf(eak(k+1),na);
        

        xout=[xout hxk1k1];
        pout=[pout Pk1k1];
        nis=[nis nisk];
        nees=[nees neesk];
        crlb=[crlb crlbb];
    
end
    
end
