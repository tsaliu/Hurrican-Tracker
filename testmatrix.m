clc; clear; close all;

a=[1 2;6 7];
b=[5 6;7 8];
a' ^-1*b^-1*a^-1
(a*b*a')^-1

z=[1;2;3;4];
% z^-1
if 1 == 1
   disp('=='); 
else
    disp('2');
end

x=1:120;
x=reshape(x,30,[])';
