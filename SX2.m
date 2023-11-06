function [SX,E,MM]=SX2(u,a,b,c,d,n,nn)
S=SG(a,b,c,d,n,nn);
S1=SE(a,b,c,d,n,nn);
S2=S(2:end-1);
E=S1(:,1:end-1)./S2;
u1=sum(u,2);
SX=0;


for k=2:n-1
    loca=find(u1==k);
    N=length(loca);
    if N==0
        M=zeros(1,17);
    else
        M=mean(u(loca,:),1);
    end
    MM(k,:)=M;
    SX=SX+N*((M'-E(:,k)).^2)./(E(:,k).*(1-E(:,k)));
end
