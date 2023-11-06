function [w,p]=wrnd1(u,a,b,c,d,th,n)
N=ones(n,1);
c1=N*c;
d1=N*d;
th1=th*a;
pro=normcdf(th1+b,0,1);
 %pro(pro>0.9999)=0.9999;
pro1=d1.*pro+c1.*(1-pro);
 %pro1(pro1>0.9999)=0.9999;
% pro1(pro1<0.0001)=0.0001;
pro2=1-pro1;
p1=(d1.*pro)./pro1;
p2=((1-d1).*pro)./pro2;
p=p1.*u+p2.*(1-u);
w=binornd(1,p);

