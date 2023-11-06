function [Ra,Rb,Rc,L]=EM3PNO(u,a,b,c,abprior,cprior,grid,ntime)
[n,m]=size(u);
N=ones(n,1);
I=1;
indice=1;
SIGMA=abprior;
MU=[2;0];
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
  %  at=at.*(at>0);
    bt=SS(2,:);
    ct=(s6+cprior(1)-1)./(s7+sum(cprior)-2);     
    ct(ct<=0)=0.0001;
    ct(ct>=1)=0.9999;
    LL(I+1)=L;
    if  abs(L-LL(end-1))<10^(-4)   
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


