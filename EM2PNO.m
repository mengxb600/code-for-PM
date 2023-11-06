function [Ra,Rb,I]=EM2PNO(u,a,b,abprior,grid,nn)
[n,m]=size(u);
I=1;
indice=1;
SIGMA=abprior;
MU=[1;0];
LL(1)=1;
while indice==1&&I<nn

    [Z1,Z2,TH1,TH2,L]=E3(u,a,b,grid);
    
    S12=sum(TH1);
    S11=sum(TH2);

    S1=[S11,S12;S12,n];
    S2=[sum(Z2);sum(Z1)];


    SS=(S1+SIGMA)^(-1)*(S2+SIGMA*MU);
    at=SS(1,:);
    at=at.*(at>0);
    bt=SS(2,:);
    
    LL(I+1)=L;


    if   abs(LL(end)-LL(end-1))<10^(-4)
         indice=0;
         a=at;
         b=bt;

     else
         a=at;
         b=bt;

         I=I+1;
     end   
end
Ra=a;
Rb=b;



