function [Ra,Rb,Rc,Rd,TH1,L,Infor]=EM4PNO1(u,a,b,c,d,pri_ab,priorc,priord,n,grid,nn)
N=ones(n,1);
I=1;
indice=1;
SIGMA=pri_ab;
MU=[0;0];
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

    ct=(s6+priorc(1)-1+0.05)./(s7+sum(priorc)-2+0.1);
    dt=(s9+priord(1)-1+0.05)./(s8+sum(priord)-2+0.1);
%     Lo=dt<ct;
%     ct(Lo)=(s6(Lo)+s9(Lo)+priorc(1)+priord(1)-2)./(s7(Lo)+sum(priorc)+s8(Lo)+sum(priord)-4);
%     dt(Lo)=ct(Lo);
    dt(dt>=0.999)=0.999;
    ct(ct<=0)=0.001;
    LL(I+1)=L;

  

    if max(max(abs([at;bt;ct;dt]-[a;b;c;d])))<10^(-4)   
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

[pr,Z1,Z2,TH1,TH2,L]=E1(u,a,b,c,d,grid);

 S12=sum(TH1);
 S11=sum(TH2);
 S1=[S11,S12;S12,n];
 S2=[sum(Z2);sum(Z1)];
 s6=sum((1-pr).*u);
 s7=sum(1-pr);
 s8=sum(pr);
 s9=sum(u.*pr);





lab=S2-S1*[a;b];
cl1=1./(c.*(1-c));
cl2=1./(c-1);
lc=sum([cl1;cl2].*[s6;s7]);
dl1=1./(d.*(1-d));
dl2=1./(d-1);
ld=sum([dl1;dl2].*[s9;s8]);
B1=[lab;lc;ld];
B10=B1(1:end);
In3=B10'*B10;
m=length(a);
lcc=sum([-1./(c.^2)+1./((1-c).^2);-1./((1-c).^2)].*[s6;s7]);
ldd=sum([-1./(d.^2)+1./((1-d).^2);-1./((1-d).^2)].*[s9;s8]);
In1=zeros(4*m);
for i=1:m
    In1(4*(i-1)+1:4*(i-1)+2,4*(i-1)+1:4*(i-1)+2)=-S1;
    In1(4*(i-1)+3,4*(i-1)+3)=lcc(i);
    In1(4*(i-1)+4,4*(i-1)+4)=ldd(i);
end
In2=In(u,a,b,c,d,20000);
Infor=-In1-In2+In3;
raa=ones(1,m);
prab=[raa;raa;4./(c.^2)+16./((1-c).^2);16./(d.^2)+4./((1-d).^2)];
Infor=Infor+diag(prab(1:end));
