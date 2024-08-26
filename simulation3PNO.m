clc
clear
% m: the test length
% n:the sample size of test takers
m=40;
n=500;
M=ones(1,m);
N=ones(n,1);
rng(202302)
% a, b,c: are the true values of item parameters
% th: the true values of persons' abilities. 
a=unifrnd(0.5,3,1,m);
b=normrnd(0,1,1,m);
c=unifrnd(0,0.4,1,m);
% prioab:the prior for a and b 
prioab=[1,0;0,1];
prioab2=0*prioab;
pric=[5,17];
pric2=[1,1];
a1=N*a;
b1=N*b;
c1=N*c;
th=normrnd(0,1,n,1);
th1=th*M;
% Generating Probability
p=c1+(1-c1).*normcdf(a1.*th1+b1,0,1);

% Starting points
a0=ones(1,m);
b0=zeros(1,m);
c0=b0+0.2;
for i=1:100
    i
    u=binornd(1,p); 
  % informative prior
    tic
    [ra1,rb1,rc1]=EM3PNO(u,a0,b0,c0,prioab,pric,n,30,500);
    ti(i,1)=toc;
  % non-informative prior
    tic
    [ra2,rb2,rc2]=EM3PNO(u,a0,b0,c0,prioab2,pric2,n,30,500);
    ti(i+100,1)=toc;
    
    
    
    tic
    [ra3,rb3,rc3]=NMSM3(u,a0,b0,c0,prioab,pric,n,1000,30,500);
    ti(i,2)=toc;
    tic
    [ra4,rb4,rc4]=NMSM3(u,a0,b0,c0,prioab2,pric2,n,1000,30,500);
    ti(i+100,2)=toc;
    
    
%     
%     tic
     [ra5,rb5,rc5]=MCEM3(u,a0,b0,c0,prioab,pric,n,30,500,20);
%     ti(i,3)=toc;
%     tic
     [ra6,rb6,rc6]=MCEM3(u,a0,b0,c0,prioab2,pric2,n,30,500,20);
%     ti(i+100,3)=toc;
    
    
     CC=corrcoef(a,ra1);
     coa1(i)=CC(2);
     CC=corrcoef(b,rb1);
     cob1(i)=CC(2);
     CC=corrcoef(c,rc1);
     coc1(i)=CC(2);
     

     CC=corrcoef(a,ra2);
     coa2(i)=CC(2);
     CC=corrcoef(b,rb2);
     cob2(i)=CC(2);
     CC=corrcoef(c,rc2);
     coc2(i)=CC(2);
     
     CC=corrcoef(a,ra3);
     coa3(i)=CC(2);
     CC=corrcoef(b,rb3);
     cob3(i)=CC(2);
     CC=corrcoef(c,rc3);
     coc3(i)=CC(2);
     
     
     CC=corrcoef(a,ra4);
     coa4(i)=CC(2);
     CC=corrcoef(b,rb4);
     cob4(i)=CC(2);
     CC=corrcoef(c,rc4);
     coc4(i)=CC(2);
     
     
     CC=corrcoef(a,ra5);
     coa5(i)=CC(2);
     CC=corrcoef(b,rb5);
     cob5(i)=CC(2);
     CC=corrcoef(c,rc5);
     coc5(i)=CC(2);
     
     
     CC=corrcoef(a,ra6);
     coa6(i)=CC(2);
     CC=corrcoef(b,rb6);
     cob6(i)=CC(2);
     CC=corrcoef(c,rc6);
     coc6(i)=CC(2);
    
    RA(i,:)=a;
    RB(i,:)=b;
    RC(i,:)=c;

    
    RA1(i,:)=ra1;
    RB1(i,:)=rb1;
    RC1(i,:)=rc1;
     
        
    RA2(i,:)=ra2;
    RB2(i,:)=rb2;
    RC2(i,:)=rc2;
     
     
    RA3(i,:)=ra3;
    RB3(i,:)=rb3;
    RC3(i,:)=rc3;
     
     
    RA4(i,:)=ra4;
    RB4(i,:)=rb4;
    RC4(i,:)=rc4;
    
    RA5(i,:)=ra5;
    RB5(i,:)=rb5;
    RC5(i,:)=rc5;
    
    RA6(i,:)=ra6;
    RB6(i,:)=rb6;
    RC6(i,:)=rc6;
end
 armse1=sqrt(mean((RA-RA1).^2));
 brmse1=sqrt(mean((RB-RB1).^2));
 crmse1=sqrt(mean((RC-RC1).^2));

 armse2=sqrt(mean((RA-RA2).^2));
 brmse2=sqrt(mean((RB-RB2).^2));
 crmse2=sqrt(mean((RC-RC2).^2));
 
 
 armse3=sqrt(mean((RA-RA3).^2));
 brmse3=sqrt(mean((RB-RB3).^2));
 crmse3=sqrt(mean((RC-RC3).^2));
 
 
 armse4=sqrt(mean((RA-RA4).^2));
 brmse4=sqrt(mean((RB-RB4).^2));
 crmse4=sqrt(mean((RC-RC4).^2));
 
 armse5=sqrt(mean((RA-RA5).^2));
 brmse5=sqrt(mean((RB-RB5).^2));
 crmse5=sqrt(mean((RC-RC5).^2));
 
 
 armse6=sqrt(mean((RA-RA6).^2));
 brmse6=sqrt(mean((RB-RB6).^2));
 crmse6=sqrt(mean((RC-RC6).^2));
 
coa=[coa1;coa2;coa3;coa4;coa5;coa6];
cob=[cob1;cob2;cob3;cob4;cob5;cob6];
coc=[coc1;coc2;coc3;coc4;coc5;coc6];
% The mean of RMSE across items
armse=[armse1;armse2;armse3;armse4;armse5;armse6];
brmse=[brmse1;brmse2;brmse3;brmse4;brmse5;brmse6];
crmse=[crmse1;crmse2;crmse3;crmse4;crmse5;crmse6];
are=mean(armse');
bre=mean(brmse');
cre=mean(crmse');
% The mean of correlations across simultions. 
ca=mean(coa');
cb=mean(cob');
cc=mean(coc');
