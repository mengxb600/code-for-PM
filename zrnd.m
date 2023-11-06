function z=zrnd(u,a,b,th0)
[N,n]=size(u);
mu=th0*a+b;
sigma=1;
y=unifrnd(0,1,N,n);
k=normcdf(0,mu,sigma);
k(k>0.9999)=0.9999;
k(k<0.0001)=0.0001;

x1=k+y.*(1-k);
x2=k.*y;
x=x1.*u+x2.*(1-u);
z=norminv(x,mu,sigma);
end