clear
clc
close all
%%
a11=1;
a12=2;
a13=3;
a21=4;
a22=5;
a23=6;
a31=7;
a32=8;
a33=9;
%%
mu1=10;
mu2=11;
mu3=12;
%%
El1=13;
El2=14;
El3=15;
%%
syms t1 t2 t3 
R1=[cos(t1) sin(t1) 0; -sin(t1) cos(t1) 0;0 0 1];
R2=[cos(t2) 0 sin(t2);0 1 0;-sin(t2) 0 cos(t2)];
R3=[1 0 0;0 cos(t3) sin(t3);0 -sin(t3) cos(t3)];
R=R1*R2*R3;
%%
pol=[a11 a12 a13;a21 a22 a23;a31 a32 a33];
MU=[mu1 mu2 mu3];
E0=[El1;El2;El3];
EL=R*E0;
EL_tra=EL.';
%%
syms j1 j2 j3 i1 i2 i3
H1=((2*j1^2)/(2*i1))+((2*j2^2)/(2*i2))+((2*j3^2)/(2*i3));
H2 =-0.5*(EL_tra*pol*EL);
H3=(-MU*EL);
H=H2+H3;
F=exp(H);
fo=@(t1,t2,t3) exp(((a13.*(El2.*(cos(t3).*sin(t1) - cos(t1).*sin(t2).*sin(t3)) + El3.*(sin(t1).*sin(t3) + cos(t1).*cos(t3).*sin(t2)) + El1.*cos(t1).*cos(t2)) - a33.*(El1.*sin(t2) - El3.*cos(t2).*cos(t3) + El2.*cos(t2).*sin(t3)) + a23.*(El2.*(cos(t1).*cos(t3) + sin(t1).*sin(t2).*sin(t3)) + El3.*(cos(t1).*sin(t3) - cos(t3).*sin(t1).*sin(t2)) - El1.*cos(t2).*sin(t1))).*(El1.*sin(t2) - El3.*cos(t2).*cos(t3) + El2.*cos(t2).*sin(t3)))/2 + mu3.*(El1.*sin(t2) - El3.*cos(t2).*cos(t3) + El2.*cos(t2).*sin(t3)) - ((El2.*(cos(t3).*sin(t1) - cos(t1).*sin(t2).*sin(t3)) + El3.*(sin(t1).*sin(t3) + cos(t1).*cos(t3).*sin(t2)) + El1.*cos(t1).*cos(t2)).*(a11.*(El2.*(cos(t3).*sin(t1) - cos(t1).*sin(t2).*sin(t3)) + El3.*(sin(t1).*sin(t3) + cos(t1).*cos(t3).*sin(t2)) + El1.*cos(t1).*cos(t2)) - a31.*(El1.*sin(t2) - El3.*cos(t2).*cos(t3) + El2.*cos(t2).*sin(t3)) + a21.*(El2.*(cos(t1).*cos(t3) + sin(t1).*sin(t2).*sin(t3)) + El3.*(cos(t1).*sin(t3) - cos(t3).*sin(t1).*sin(t2)) - El1.*cos(t2).*sin(t1))))/2 - ((El2.*(cos(t1).*cos(t3) + sin(t1).*sin(t2).*sin(t3)) + El3.*(cos(t1).*sin(t3) - cos(t3).*sin(t1).*sin(t2)) - El1.*cos(t2).*sin(t1)).*(a12.*(El2.*(cos(t3).*sin(t1) - cos(t1).*sin(t2).*sin(t3)) + El3.*(sin(t1).*sin(t3) + cos(t1).*cos(t3).*sin(t2)) + El1.*cos(t1).*cos(t2)) - a32.*(El1.*sin(t2) - El3.*cos(t2).*cos(t3) + El2.*cos(t2).*sin(t3)) + a22.*(El2.*(cos(t1).*cos(t3) + sin(t1).*sin(t2).*sin(t3)) + El3.*(cos(t1).*sin(t3) - cos(t3).*sin(t1).*sin(t2)) - El1.*cos(t2).*sin(t1))))/2 - mu1.*(El2.*(cos(t3).*sin(t1) - cos(t1).*sin(t2).*sin(t3)) + El3.*(sin(t1).*sin(t3) + cos(t1).*cos(t3).*sin(t2)) + El1.*cos(t1).*cos(t2)) - mu2.*(El2.*(cos(t1).*cos(t3) + sin(t1).*sin(t2).*sin(t3)) + El3.*(cos(t1).*sin(t3) - cos(t3).*sin(t1).*sin(t2)) - El1.*cos(t2).*sin(t1)))

q=integral3(fo,0,2*pi,0,2*pi,0,2*pi);