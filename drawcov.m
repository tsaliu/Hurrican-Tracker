function h = drawcov(x,y,xr,yr)
hold on
% load mapping;
th = 0:pi/50:2*pi;
    for i=1:1:length(x)
        xunit = xr(i) * cos(th) + x(i);
        yunit = yr(i) * sin(th) + y(i);
        h = plot(xunit, yunit,'b');
    end

end