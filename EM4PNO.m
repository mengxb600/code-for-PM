function [Ra,Rb,Rc,Rd,L]=EM4PNO(u,a,b,c,d,pri_ab,priorc,priord,n,grid,nn)
N=ones(n,1);
I=1;
indice=1;
SIGMA=pri_ab;
MU=[2;0];
LL(1)=1;
while indice==1&&I<nn
    [pr,Z1,Z2,TH1,TH2,L]=E1(u,a,b,c,d,grid);

    S12=sum(TH1);
    S11=sum(TH2);
    S1=[S11,S12;S12,n];
    S2=[sum(Z2);sum(Z1)];
    s6=sum((1-pr).*u);
    s7=sum(1-pr);
    s8=sum(pr);
    s9=sum(u.*pr);  
    
    SS=(S1+SIGMA)^(-1)*(S2+SIGMA*MU);
    at=SS(1,:);
    at=at.*(at>0);
    bt=SS(2,:);

    ct=(s6+priorc(1)-1)./(s7+sum(priorc)-2);
    dt=(s9+priord(1)-1)./(s8+sum(priord)-2);
    Lo=dt<ct;
    ct(Lo)=(s6(Lo)+s9(Lo)+priorc(1)+priord(1)-2)./(s7(Lo)+sum(priorc)+s8(Lo)+sum(priord)-4);
    dt(Lo)=ct(Lo);
    dt(dt>=0.999)=0.999;
    ct(ct<=0)=0.001;
    LL(I+1)=L;

  
 
%     ct(1)=0;
%     dt(1)=1;
%     at(1)=1;
%     bt(1)=0;
    if abs(L-LL(end-1))<10^(-4)   
         indice=0;
         a=at;
         b=bt;
         c=ct;
         d=dt;
     else
         a=at;
         b=bt;
         c=ct;
         d=dt;
         I=I+1;
     end   
end
Ra=a;
Rb=b;
Rc=c;
Rd=d;

