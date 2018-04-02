function [lat,lon]=covlatlon(x,y,z)

r=6378100;  %radius of earth in m
latr=asin(z./r);
lonr=atan(y./x);

lat=90-(latr.*180./pi);
lon=lonr.*180./pi;

end
