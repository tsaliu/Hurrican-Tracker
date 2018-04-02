function [x,y,z,dx,dy,dz,T]=latloncov(date,time,lat,lon)

    latr=lat*pi/180;
    lonr=lon*pi/180;
    r=6378100;   %raduis of earth in m

    rho=r*sin(latr);
    x=rho.*cos(lonr);
    y=rho.*sin(lonr);
    z=r*cos(latr);

    
    for i=1:size(time)-1
        if ((date(i+1)-date(i))*100)>1
            date(i)=date(i+1)-0.01;
        end
        if (time(i+1)-time(i)) > 0
            T(i)=(date(i+1)-date(i))*100*24+time(i+1)-time(i);
        end
        if (time(i+1)-time(i)) < 0
            T(i)=(date(i+1)-date(i))*100*24+time(i+1)-time(i);
        end
        dx(i)=(x(i+1)-x(i))./T(i);
        dy(i)=(y(i+1)-y(i))./T(i);
        dz(i)=(z(i+1)-z(i))./T(i);
    end
%     for j=1:size(x)-1
%         dx=(x(j+1)-x(j))./T;
%         dy=(y(j+1)-y(j))./T;
%         dz=(z(j+1)-z(j))./T;
%     end
end