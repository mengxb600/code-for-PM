function [Ra,Rb,Rc,TH1,L,Infor]=EM3PNO1(u,a,b,c,abprior,cprior,grid,ntime)
[n,m]=size(u);
N=ones(n,1);
I=1;
indice=1;
SIGMA=abprior;
MU=[0;0];
LL(1)=1;
d=ones(1,length(a));
while indice==1&&I<ntime

    [pr,Z1,Z2,TH1,TH2,L]=E2(u,a,b,c,grid);
    
    S12=sum(TH1);
    S11=sum(TH2);
    S1=[S11,S12;S12,n];
    S2=[sum(Z2);sum(Z1)];
    s6=sum((1-pr).*u);
    s7=sum(1-pr);

    SS=(S1+SIGMA)^(-1)*(S2+SIGMA*MU);
    at=SS(1,:);
    at=at.*(at>0);
    bt=SS(2,:);
    ct=(s6+cprior(1)-1)./(s7+sum(cprior)-2);     
    ct(ct<=0)=0.0001;
    ct(ct>=1)=0.999;
    LL(I+1)=L;


    if  max(max(abs([at;bt;ct]-[a;b;c])))<10^(-4)   
         
         indice=0;
         a=at;
         b=bt;
         c=ct;
     else
         a=at;
         b=bt;
         c=ct;
         I=I+1;
     end   
end
Ra=a;
Rb=b;
Rc=c;

[pr,Z1,Z2,TH1,TH2,L]=E2(u,a,b,c,grid);

S12=sum(TH1);
S11=sum(TH2);
S1=[S11,S12;S12,n];
S2=[sum(Z2);sum(Z1)];
s6=sum((1-pr).*u);
s7=sum(1-pr);



lab=S2-S1*[a;b];
cl1=1./(c.*(1-c));
cl2=1./(c-1);
lc=sum([cl1;cl2].*[s6;s7]);
B1=[lab;lc];
B10=B1(1:end);
In3=B10'*B10;


m=length(a);
lcc=sum([-1./(c.^2)+1./((1-c).^2);-1./((1-c).^2)].*[s6;s7]);
In1=zeros(3*m);
for i=1:m
    In1(3*(i-1)+1:3*(i-1)+2,3*(i-1)+1:3*(i-1)+2)=-S1;
    In1(3*(i-1)+3,3*(i-1)+3)=lcc(i);
end


B0=IB3(u,a,b,c,20000);
In2=B0;
Infor=-In1-In2+In3;
raa=ones(1,m);
prab=[raa;raa;4./(c.^2)+16./((1-c).^2)];
prab(1:end)
Infor=Infor+diag(prab(1:end));
