clc
clear
% m: the test length
% n:the sample size of test takers
m=20;
n=500;
M=ones(1,m);
N=ones(n,1);

rng(202304)
% a, b, c, d: are the true values of item parameters
% th: the true values of persons' abilities. 
a=unifrnd(0.5,3,1,m);
b=normrnd(0,1,1,m);
c=unifrnd(0,0.4,1,m);
d=1-unifrnd(0,0.3,1,m);
th=normrnd(0,1,n,1);
% pri_ab: the prior for (a, b)
pri_ab=inv([1,0;0,1]);
a1=N*a;
b1=N*b; 
c1=N*c;
d1=N*d;
th1=th*M;
% Generating Probability
p=c1+(d1-c1).*normcdf(a1.*th1+b1,0,1);
% Starting points
a0=ones(1,m);
b0=zeros(1,m);
c0=b0+0.2;
d0=b0+0.8;

for i=1:100
    i
    u=binornd(1,p);
    % Under informative priors
    [ra1,rb1,rc1,rd1,I1]=EM4PNO(u,a0,b0,c0,d0,pri_ab,[5,17],[17,5],n,30,1000);
    % Under non-informative priors
    [ra2,rb2,rc2,rd2,I2(i)]=EM4PNO(u,a0,b0,c0,d0,0*pri_ab,[1,1],[1,1],n,30,1000);

 
    RA(i,:)=a;
    RB(i,:)=b;
    RC(i,:)=c;
    RD(i,:)=d;

    
    
     RA1(i,:)=ra1;
     RB1(i,:)=rb1;
     RC1(i,:)=rc1;
     RD1(i,:)=rd1;

    
    
    
    RA2(i,:)=ra2;
    RB2(i,:)=rb2;
    RC2(i,:)=rc2;
    RD2(i,:)=rd2;
    
    % Correlation coefficients between the estimates and the true values
    CC=corrcoef(a,ra1);
    coa1(i)=CC(2);
    CC=corrcoef(b,rb1);
    cob1(i)=CC(2);
    CC=corrcoef(c,rc1);
    coc1(i)=CC(2);
    CC=corrcoef(d,rd1);
    cod1(i)=CC(2);
 
    
    CC=corrcoef(a,ra2);
    coa2(i)=CC(2);
    CC=corrcoef(b,rb2);
    cob2(i)=CC(2);
    CC=corrcoef(c,rc2);
    coc2(i)=CC(2);
    CC=corrcoef(d,rd2);
    cod2(i)=CC(2);
       
       


end
coa=[coa1;coa2];
cob=[cob1;cob2];
coc=[coc1;coc2];
cod=[cod1;cod2];

% The RMSE
armse1=sqrt(mean((RA-RA1).^2));
brmse1=sqrt(mean((RB-RB1).^2));
crmse1=sqrt(mean((RC-RC1).^2));
drmse1=sqrt(mean((RD-RD1).^2));

armse2=sqrt(mean((RA-RA2).^2));
brmse2=sqrt(mean((RB-RB2).^2));
crmse2=sqrt(mean((RC-RC2).^2));
drmse2=sqrt(mean((RD-RD2).^2));

  


armse=[armse1;armse2];
brmse=[brmse1;brmse2];
crmse=[crmse1;crmse2];
drmse=[drmse1;drmse2];
% The mean of RMSE across items
are=mean(armse');
bre=mean(brmse');
cre=mean(crmse');
dre=mean(drmse');
% The mean of correlations across simultions. 
ca=mean(coa');
cb=mean(cob');
cc=mean(coc');
cd=mean(cod');


