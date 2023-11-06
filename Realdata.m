clc
clear
tic
m=17;
n=1544;
M=ones(1,m);
N=ones(n,1);
load pisa
% pisa is the data
a0=ones(1,m);
b0=zeros(1,m);
c0=b0+0.2;
d0=b0+0.85; 
% a0,b0,c0,d0 are the initial values of EEM
pri_ab=inv([1,0;0,1]);
% EM4PNO is the EEM for 4PNO
% EM3PNO is the EEM for 3PNO
% EM2PNO is the EEM for 2PNO
[ra1,rb1,rc1,rd1,rth1,l1,Infor4]=EM4PNO1(u,a0,b0,c0,d0,pri_ab,[5,17],[17,5],n,30,1000);
[ra2,rb2,rc2,rth2,l2,Infor3]=EM3PNO1(u,a0,b0,c0,pri_ab,[5,17],30,1000);
[ra3,rb3,rth3,l3,Infor2]=EM2PNO1(u,a0,b0,pri_ab,30,1000);




% AIC4 is the AIC for 4PNO
% AIC3 is the AIC for 3PNO
% AIC2 is the AIC for 2PNO
% SXX4 is the SX2 for 4PNO
% SXX3 is the SX2 for 3PNO
% SXX2 is the SX2 for 2PNO
SXX4=SX2(u,ra1,rb1,rc1,rd1,17,30);
SXX3=SX2(u,ra2,rb2,rc2,rd1*0+1,17,30);
SXX2=SX2(u,ra3,rb3,rb2*0,rb2*0+1,17,30);
AIC4=-2*l1+4*17;
AIC3=-2*l2+3*17;
AIC2=-2*l3+2*17;

se4=sqrt(abs(diag(inv(Infor4))));
se3=sqrt(abs(diag(inv(Infor3))));
se2=sqrt(abs(diag(inv(Infor2))));






