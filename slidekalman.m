function [hxk1k1,Pk1k1,nis,nees,crlb]=slidekalman(P00,T,x,y,dx,dy,N)

hxk1k1=[];
Pk1k1=[];
nis=[];
nees=[];
crlb=[];
x(1)=[];
y(1)=[];
T=reshape(T,N,[])';
x=reshape(x,N,[])';
y=reshape(y,N,[])';
dx=reshape(dx,N,[])';
dy=reshape(dy,N,[])';
[a,b]=size(T);

for k=1:a
    [hxk1k1t,Pk1k1t,nist,neest,crlbt]=kalman(P00,T(k,:),x(k,:),y(k,:),dx(k,:),dy(k,:));
    hxk1k1=[hxk1k1 hxk1k1t];
    Pk1k1=[Pk1k1 Pk1k1t];
    nis=[nis nist];
    nees=[nees neest];
    crlb=[crlb crlbt];
end



end