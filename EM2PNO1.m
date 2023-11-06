function [Ra,Rb,TH1,L,Infor]=EM2PNO1(u,a,b,abprior,grid,ntime)
[n,m]=size(u);
I=1;
indice=1;
SIGMA=abprior;
MU=[0;0];
LL(1)=1;
while indice==1&&I<ntime

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


    if   max(max(abs([at;bt]-[a;b])))<10^(-4)
        
        %abs(LL(end)-LL(end-1))<10^(-4)
       
        
         
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

[Z1,Z2,TH1,TH2,L]=E3(u,a,b,grid);

S12=sum(TH1);
S11=sum(TH2);
S1=[S11,S12;S12,n];

S2=[sum(Z2);sum(Z1)];


lab=S2-S1*[a;b];
B1=[lab];
B10=B1(1:end);
In3=B10'*B10;


m=length(a);
In1=zeros(2*m);
for i=1:m
    In1(2*(i-1)+1:2*(i-1)+2,2*(i-1)+1:2*(i-1)+2)=-S1;
end


In2=IB2(u,a,b,2000);
Infor=-In1-In2+In3;



