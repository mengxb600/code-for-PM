clc
clear
tic
% m: the test length
% n:the sample size of test takers
m=20;
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
prioab=[1,0;0,1]^(-1);


th=normrnd(0,1,n,1);
th1=th*M;
% Generating Probability
p=normcdf(a1.*th1+b1,0,1);
% Starting points
a0=ones(1,m);
b0=zeros(1,m);
rng(202302)
for i=1:100
   i

    u=binornd(1,p);  
% informative prior
    [ra1,rb1]=EM2PNO(u,a0,b0,prioab,30,1000);
% non-informative prior
    [ra2,rb2]=EM2PNO(u,a0,b0,0*prioab,30,1000);



     CC=corrcoef(a,ra1);
     coa1(i)=CC(2);
     CC=corrcoef(b,rb1);
     cob1(i)=CC(2);
   
     CC=corrcoef(a,ra2);
     coa2(i)=CC(2);
     CC=corrcoef(b,rb2);
     cob2(i)=CC(2);

     

    RA(i,:)=a;
    RB(i,:)=b;


    
    RA1(i,:)=ra1;
    RB1(i,:)=rb1;
    
    RA2(i,:)=ra2;
    RB2(i,:)=rb2;
end
coa=[coa1;coa2];
cob=[cob1;cob2];


% The RMSE
armse1=sqrt(mean((RA-RA1).^2));
brmse1=sqrt(mean((RB-RB1).^2));


armse2=sqrt(mean((RA-RA2).^2));
brmse2=sqrt(mean((RB-RB2).^2));

armse=[armse1;armse2];
brmse=[brmse1;brmse2];

% The mean of RMSE across items
are=mean(armse');
bre=mean(brmse');

% The mean of correlations across simultions. 
ca=mean(coa');
cb=mean(cob');
toc

