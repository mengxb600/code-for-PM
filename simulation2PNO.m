clc
clear
tic
% m: the test length
% n:the sample size of test takers
m=40;
n=500;
M=ones(1,m);
N=ones(n,1);
rng(1)
% a, b: are the true values of item parameters
% th: the true values of persons' abilities. 
a=unifrnd(0.5,3,1,m);
b=normrnd(0,1,1,m);
a1=N*a;
b1=N*b;
% prioab:the prior for a and b 
prioab=[1,0;0,1];
th=normrnd(0,1,n,1);
th1=th*M;
% Generating Probability
p=normcdf(a1.*th1+b1,0,1);
% Starting points
a0=ones(1,m);
b0=zeros(1,m);
rng(202302)
SIGMA=[1,0;0,1];
for i=1:100
   i
    u=binornd(1,p);  
% informative prior
    tic
    [ra1,rb1]=EM2PNO(u,a0,b0,prioab,n,30,500);
    ti(i,1)=toc;
% non-informative prior
    tic
    [ra2,rb2]=EM2PNO(u,a0,b0,0*prioab,n,30,500);
    ti(i+100,1)=toc;
    
    tic
    [ra3,rb3]=NMSM2(u,a0,b0,SIGMA,n,30,1000,500);
    ti(i,2)=toc;
    tic
    [ra4,rb4]=NMSM2(u,a0,b0,0*SIGMA,n,30,1000,500);
    ti(i+100,2)=toc;
    
     tic
     [ra5,rb5]=MCEM2(u,a0,b0,SIGMA,n,30,500,20);
    ti(i,3)=toc;
     tic
     [ra6,rb6]=MCEM2(u,a0,b0,0*SIGMA,n,30,500,20);
     ti(i+100,3)=toc;

     CC=corrcoef(a,ra1);
     coa1(i)=CC(2);
     CC=corrcoef(b,rb1);
     cob1(i)=CC(2);

     

     CC=corrcoef(a,ra2);
     coa2(i)=CC(2);
     CC=corrcoef(b,rb2);
     cob2(i)=CC(2);

     
     CC=corrcoef(a,ra3);
     coa3(i)=CC(2);
     CC=corrcoef(b,rb3);
     cob3(i)=CC(2);

     
     CC=corrcoef(a,ra4);
     coa4(i)=CC(2);
     CC=corrcoef(b,rb4);
     cob4(i)=CC(2);

     
     
      CC=corrcoef(a,ra5);
      coa5(i)=CC(2);
      CC=corrcoef(b,rb5);
      cob5(i)=CC(2);
% 
%      
%      
      CC=corrcoef(a,ra6);
      coa6(i)=CC(2);
      CC=corrcoef(b,rb6);
      cob6(i)=CC(2);

    
    RA(i,:)=a;
    RB(i,:)=b;


    
    RA1(i,:)=ra1;
    RB1(i,:)=rb1;

     
        
    RA2(i,:)=ra2;
    RB2(i,:)=rb2;
 
     
     
    RA3(i,:)=ra3;
    RB3(i,:)=rb3;

     
     
    RA4(i,:)=ra4;
    RB4(i,:)=rb4;

    
     RA5(i,:)=ra5;
     RB5(i,:)=rb5;
 
     
     RA6(i,:)=ra6;
     RB6(i,:)=rb6;


end
coa=[coa1;coa2;coa3;coa4;coa5;coa6];
cob=[cob1;cob2;cob3;cob4;cob5;cob6];


% The RMSE
armse1=sqrt(mean((RA-RA1).^2));
brmse1=sqrt(mean((RB-RB1).^2));


armse2=sqrt(mean((RA-RA2).^2));
brmse2=sqrt(mean((RB-RB2).^2));

armse3=sqrt(mean((RA-RA3).^2));
brmse3=sqrt(mean((RB-RB3).^2));


armse4=sqrt(mean((RA-RA4).^2));
brmse4=sqrt(mean((RB-RB4).^2));

armse5=sqrt(mean((RA-RA5).^2));
brmse5=sqrt(mean((RB-RB5).^2));

armse6=sqrt(mean((RA-RA6).^2));
brmse6=sqrt(mean((RB-RB6).^2));

armse=[armse1;armse2;armse3;armse4;armse5;armse6];
brmse=[brmse1;brmse2;brmse3;brmse4;brmse4;brmse6];

% The mean of RMSE across items
are=mean(armse');
bre=mean(brmse');

% The mean of correlations across simultions. 
ca=mean(coa');
cb=mean(cob');
toc

