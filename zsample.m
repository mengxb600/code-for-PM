function z=zsample(u,a,b,c,d,th0)
[m,n]=size(u);
mu=th0*a+b;
p0=normcdf(mu,0,1);
p=c+(d-c).*p0;
y=unifrnd(0,1,m,n);
t=(1-u-c).*(1-p0)./(1-u-p);
k=y>t;
z=norminv(((1-u-p).*y+k.*(c-d).*(1-p0))./(1-u-(c+k.*(d-c))))+mu;

