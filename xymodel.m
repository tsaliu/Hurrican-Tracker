function [x,dx,y,dy,T]=xymodel(date,time,xmap,ymap)

x=xmap;
y=ymap;

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
    end
end