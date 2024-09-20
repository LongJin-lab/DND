clear all
close all
global ita1;
global ita2;
global L;
global Lambda;
global n;
global y00
global c1;
global w;
global l;
global c0;
global c2
global c3
global basea
global y0jl
global d

d=0.2;
c0=50;
ts=10;
c1=100; 
c2=20000;
c3=200;
n=6;
Aa=eye(42);
P=eye(42);
jla=zeros(n,n);
    for j=1:n
          jla(j,j) = 0;
    end
jla(1,2) = 1;jla(2,3) = 1;jla(3,6) = 1;jla(5,6) = 1;jla(4,5) = 1;jla(1,4) = 1; 
jla(2,1) = 1;jla(3,2) = 1;jla(6,3) = 1;jla(6,5) = 1;jla(5,4) = 1;jla(4,1) = 1;     
L=diag(jla*ones(n,1))-jla;   
Lambda=zeros(n,n);
Lambda(1,1)=1;

Lambda(4,4)=1;    
basea(1,:)=[0;0;0;0];
basea(2,:)=[0;0;0;0];
basea(3,:)=[0;0;0;0];
basea(4,:)=[0;0;0;0];
basea(5,:)=[0;0;0;0];
basea(6,:)=[0;0;0;0];
y0jl=[0.201459192549988
0.325000174653958
0.0118213136360045
0.867043422712064
0.115312160175345
-0.541666790724753
0.0146884702391163
0.206847237638700
0.324066510260604
0.000958329240331745
0.870614819763741
0.0460444987881901
-0.534954857657347
0.0585912355105916
0.225211111587345
0.324042747925226
-0.000222545870472021
0.870345803826861
0.0401216889347861
-0.535527818203203
0.0301602931126907
0.204817426485698
0.324636838280901
0.00555192302382549
0.865529736295158
0.0759030209588581
-0.545192246159201
0.0599797772530693
0.234864171351792
0.323824317260063
-0.00180001625518178
0.870631526199001
0.0324674746167191
-0.534972629147689
0.0596173629011659
0.216523163078061
0.323667520457398
-0.00149328011255287
0.871013104821733
0.0329062933410458
-0.534194978625258
0.0371586702128454];

Y011(:,:,1)=[y0jl+0*rand(7*n,1);zeros(28*n,1)];
dt=0.001;
t_store=0:dt:ts;
myTheta(:,:,1)=Y011(1:7*n,:,1);
mw(:,:,1)=Y011(7*n+1:10*n,:,1);
myDTheta(:,:,1)=Y011(10*n+1:17*n,:,1);
myDw(:,:,1)=Y011(17*n+1:20*n,:,1);
lambda1(:,:,1)=Y011(20*n+1:23*n,:,1);
lambda2(:,:,1)=Y011(23*n+1:26*n,:,1);
lambda3(:,:,1)=Y011(26*n+1:28*n,:,1);
% lambda4(:,:,1)=Y011(40*n+1:42*n,:,1);
myDDTheta(:,:,1)=Y011(28*n+1:35*n,:,1);
ddTheta(:,:,1)=zeros(18*n,1);
sumee(:,:,1)=zeros(18*n,1);
Y01=zeros(42*n,0);
tau_store=zeros(7*n,0);
error_Estore=zeros(2*n,0);
x1(:,:,1)=zeros(42,1);
error_tau=zeros(42,1);
error_tau_n=zeros(42,1);
 for m= 1:length(t_store)-1
     m
TT=t_store(m);
T=2*pi; 

a=4;
b=1;
rx=0.005*(6*(cos(pi/5*a*TT/T*10-pi/2))+5*cos(pi/5*b*TT/T*10-pi/2));
ry=0.005*(-6*(sin(pi/5*a*TT/T*10-pi/2))+5*sin(pi/5*b*TT/T*10-pi/2)-1);
rz=0.785104;
rd=[rx;ry;rz];
myPostionstore=zeros(3,0);

for i=1:n
    theta1=myTheta((7*(i-1)+1),:,m);
    theta2=myTheta((7*(i-1)+2),:,m);
    theta3=myTheta((7*(i-1)+3),:,m);
    theta4=myTheta((7*(i-1)+4),:,m);
    theta5=myTheta((7*(i-1)+5),:,m);
    theta6=myTheta((7*(i-1)+6),:,m);
    theta7=myTheta((7*(i-1)+7),:,m);
    
    jlq=myTheta((7*(i-1)+1):(7*(i-1)+7),:,m) ;
    jldq=myDTheta((7*(i-1)+1):(7*(i-1)+7),:,m) ;
    jlqdq=[jlq;jldq];
    dotA=dotj(jlqdq);
    dotB=dotjw(jlqdq);
    theta1_0=y0jl((7*(i-1)+1));
    theta2_0=y0jl((7*(i-1)+2));
    theta3_0=y0jl((7*(i-1)+3));
    theta4_0=y0jl((7*(i-1)+4));
    theta5_0=y0jl((7*(i-1)+5));
    theta6_0=y0jl((7*(i-1)+6));
    theta7_0=y0jl((7*(i-1)+7));
                                                                                                                                                                                                                                                                                                                                                            
A =...
[   (11*sin(theta6)*(sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) - cos(theta4)*sin(theta1)*sin(theta2)))/125 - (79*sin(theta1)*sin(theta2))/250 - (107*cos(theta6)*(sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) - cos(theta4)*sin(theta1)*sin(theta2)))/1000 - (33*cos(theta1)*sin(theta3))/400 + (33*cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)))/400 - (11*cos(theta6)*(cos(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta1)*sin(theta2)*sin(theta4)) + sin(theta5)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))))/125 + (48*sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)))/125 - (107*sin(theta6)*(cos(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta1)*sin(theta2)*sin(theta4)) + sin(theta5)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))))/1000 - (48*cos(theta4)*sin(theta1)*sin(theta2))/125 + (33*sin(theta1)*sin(theta2)*sin(theta4))/400 - (33*cos(theta2)*cos(theta3)*sin(theta1))/400, (79*cos(theta1)*cos(theta2))/250 - (107*cos(theta6)*(cos(theta1)*cos(theta2)*cos(theta4) + cos(theta1)*cos(theta3)*sin(theta2)*sin(theta4)))/1000 + (11*sin(theta6)*(cos(theta1)*cos(theta2)*cos(theta4) + cos(theta1)*cos(theta3)*sin(theta2)*sin(theta4)))/125 + (11*cos(theta6)*(cos(theta5)*(cos(theta1)*cos(theta2)*sin(theta4) - cos(theta1)*cos(theta3)*cos(theta4)*sin(theta2)) + cos(theta1)*sin(theta2)*sin(theta3)*sin(theta5)))/125 + (107*sin(theta6)*(cos(theta5)*(cos(theta1)*cos(theta2)*sin(theta4) - cos(theta1)*cos(theta3)*cos(theta4)*sin(theta2)) + cos(theta1)*sin(theta2)*sin(theta3)*sin(theta5)))/1000 + (48*cos(theta1)*cos(theta2)*cos(theta4))/125 - (33*cos(theta1)*cos(theta3)*sin(theta2))/400 - (33*cos(theta1)*cos(theta2)*sin(theta4))/400 + (33*cos(theta1)*cos(theta3)*cos(theta4)*sin(theta2))/400 + (48*cos(theta1)*cos(theta3)*sin(theta2)*sin(theta4))/125, (33*cos(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)))/400 - (33*cos(theta3)*sin(theta1))/400 + (48*sin(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)))/125 + (11*cos(theta6)*(sin(theta5)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta4)*cos(theta5)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))))/125 + (107*sin(theta6)*(sin(theta5)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta4)*cos(theta5)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))))/1000 - (107*cos(theta6)*sin(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)))/1000 + (11*sin(theta4)*sin(theta6)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)))/125 - (33*cos(theta1)*cos(theta2)*sin(theta3))/400, (11*sin(theta6)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta1)*sin(theta2)*sin(theta4)))/125 - (107*cos(theta6)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta1)*sin(theta2)*sin(theta4)))/1000 + (48*cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))/125 - (33*sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))/400 - (48*cos(theta1)*sin(theta2)*sin(theta4))/125 + (11*cos(theta5)*cos(theta6)*(sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + cos(theta1)*cos(theta4)*sin(theta2)))/125 + (107*cos(theta5)*sin(theta6)*(sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + cos(theta1)*cos(theta4)*sin(theta2)))/1000 - (33*cos(theta1)*cos(theta4)*sin(theta2))/400,   (11*cos(theta6)*(sin(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta1)*sin(theta2)*sin(theta4)) - cos(theta5)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))))/125 + (107*sin(theta6)*(sin(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta1)*sin(theta2)*sin(theta4)) - cos(theta5)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))))/1000, (11*cos(theta6)*(sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + cos(theta1)*cos(theta4)*sin(theta2)))/125 + (107*sin(theta6)*(sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + cos(theta1)*cos(theta4)*sin(theta2)))/1000 - (107*cos(theta6)*(cos(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta1)*sin(theta2)*sin(theta4)) + sin(theta5)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))))/1000 + (11*sin(theta6)*(cos(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta1)*sin(theta2)*sin(theta4)) + sin(theta5)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))))/125, 0;...
 (79*cos(theta1)*sin(theta2))/250 - (33*sin(theta1)*sin(theta3))/400 - (107*cos(theta6)*(sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + cos(theta1)*cos(theta4)*sin(theta2)))/1000 + (11*sin(theta6)*(sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + cos(theta1)*cos(theta4)*sin(theta2)))/125 - (11*cos(theta6)*(cos(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta1)*sin(theta2)*sin(theta4)) + sin(theta5)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))))/125 + (33*cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))/400 - (107*sin(theta6)*(cos(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta1)*sin(theta2)*sin(theta4)) + sin(theta5)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))))/1000 + (48*sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))/125 - (33*cos(theta1)*sin(theta2)*sin(theta4))/400 + (33*cos(theta1)*cos(theta2)*cos(theta3))/400 + (48*cos(theta1)*cos(theta4)*sin(theta2))/125, (11*cos(theta6)*(cos(theta5)*(cos(theta2)*sin(theta1)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta1)*sin(theta2)) + sin(theta1)*sin(theta2)*sin(theta3)*sin(theta5)))/125 + (79*cos(theta2)*sin(theta1))/250 + (107*sin(theta6)*(cos(theta5)*(cos(theta2)*sin(theta1)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta1)*sin(theta2)) + sin(theta1)*sin(theta2)*sin(theta3)*sin(theta5)))/1000 - (107*cos(theta6)*(cos(theta2)*cos(theta4)*sin(theta1) + cos(theta3)*sin(theta1)*sin(theta2)*sin(theta4)))/1000 + (11*sin(theta6)*(cos(theta2)*cos(theta4)*sin(theta1) + cos(theta3)*sin(theta1)*sin(theta2)*sin(theta4)))/125 - (33*cos(theta3)*sin(theta1)*sin(theta2))/400 - (33*cos(theta2)*sin(theta1)*sin(theta4))/400 + (48*cos(theta2)*cos(theta4)*sin(theta1))/125 + (33*cos(theta3)*cos(theta4)*sin(theta1)*sin(theta2))/400 + (48*cos(theta3)*sin(theta1)*sin(theta2)*sin(theta4))/125, (33*cos(theta1)*cos(theta3))/400 - (33*cos(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)))/400 - (48*sin(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)))/125 - (11*cos(theta6)*(sin(theta5)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) - cos(theta4)*cos(theta5)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))))/125 - (107*sin(theta6)*(sin(theta5)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) - cos(theta4)*cos(theta5)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))))/1000 - (33*cos(theta2)*sin(theta1)*sin(theta3))/400 + (107*cos(theta6)*sin(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)))/1000 - (11*sin(theta4)*sin(theta6)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)))/125, (107*cos(theta6)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta1)*sin(theta2)*sin(theta4)))/1000 - (11*sin(theta6)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta1)*sin(theta2)*sin(theta4)))/125 - (48*cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)))/125 + (33*sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)))/400 - (33*cos(theta4)*sin(theta1)*sin(theta2))/400 - (48*sin(theta1)*sin(theta2)*sin(theta4))/125 - (11*cos(theta5)*cos(theta6)*(sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) - cos(theta4)*sin(theta1)*sin(theta2)))/125 - (107*cos(theta5)*sin(theta6)*(sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) - cos(theta4)*sin(theta1)*sin(theta2)))/1000, - (11*cos(theta6)*(sin(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta1)*sin(theta2)*sin(theta4)) - cos(theta5)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))))/125 - (107*sin(theta6)*(sin(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta1)*sin(theta2)*sin(theta4)) - cos(theta5)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))))/1000, (107*cos(theta6)*(cos(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta1)*sin(theta2)*sin(theta4)) + sin(theta5)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))))/1000 - (107*sin(theta6)*(sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) - cos(theta4)*sin(theta1)*sin(theta2)))/1000 - (11*cos(theta6)*(sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) - cos(theta4)*sin(theta1)*sin(theta2)))/125 - (11*sin(theta6)*(cos(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta1)*sin(theta2)*sin(theta4)) + sin(theta5)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))))/125, 0;...
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                 (33*sin(theta2)*sin(theta4))/400 - (33*cos(theta2)*cos(theta3))/400 - (48*cos(theta4)*sin(theta2))/125 - (79*sin(theta2))/250 - (11*cos(theta6)*(cos(theta5)*(sin(theta2)*sin(theta4) + cos(theta2)*cos(theta3)*cos(theta4)) - cos(theta2)*sin(theta3)*sin(theta5)))/125 - (107*sin(theta6)*(cos(theta5)*(sin(theta2)*sin(theta4) + cos(theta2)*cos(theta3)*cos(theta4)) - cos(theta2)*sin(theta3)*sin(theta5)))/1000 + (107*cos(theta6)*(cos(theta4)*sin(theta2) - cos(theta2)*cos(theta3)*sin(theta4)))/1000 - (11*sin(theta6)*(cos(theta4)*sin(theta2) - cos(theta2)*cos(theta3)*sin(theta4)))/125 + (33*cos(theta2)*cos(theta3)*cos(theta4))/400 + (48*cos(theta2)*cos(theta3)*sin(theta4))/125,                                                                                                                                                                                                                                                                                                                                                                                (33*sin(theta2)*sin(theta3))/400 + (11*cos(theta6)*(cos(theta3)*sin(theta2)*sin(theta5) + cos(theta4)*cos(theta5)*sin(theta2)*sin(theta3)))/125 + (107*sin(theta6)*(cos(theta3)*sin(theta2)*sin(theta5) + cos(theta4)*cos(theta5)*sin(theta2)*sin(theta3)))/1000 - (33*cos(theta4)*sin(theta2)*sin(theta3))/400 - (48*sin(theta2)*sin(theta3)*sin(theta4))/125 + (107*cos(theta6)*sin(theta2)*sin(theta3)*sin(theta4))/1000 - (11*sin(theta2)*sin(theta3)*sin(theta4)*sin(theta6))/125,                                                                                                                                                                                                                                                                                                                         (107*cos(theta6)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)))/1000 - (48*cos(theta2)*sin(theta4))/125 - (33*cos(theta2)*cos(theta4))/400 - (11*sin(theta6)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)))/125 - (33*cos(theta3)*sin(theta2)*sin(theta4))/400 + (11*cos(theta5)*cos(theta6)*(cos(theta2)*cos(theta4) + cos(theta3)*sin(theta2)*sin(theta4)))/125 + (107*cos(theta5)*sin(theta6)*(cos(theta2)*cos(theta4) + cos(theta3)*sin(theta2)*sin(theta4)))/1000 + (48*cos(theta3)*cos(theta4)*sin(theta2))/125,                                                                                                                                                                                         - (11*cos(theta6)*(sin(theta5)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) - cos(theta5)*sin(theta2)*sin(theta3)))/125 - (107*sin(theta6)*(sin(theta5)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) - cos(theta5)*sin(theta2)*sin(theta3)))/1000,                                                                                                                                                                                                                                                                                                 (107*cos(theta6)*(cos(theta5)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) + sin(theta2)*sin(theta3)*sin(theta5)))/1000 - (11*sin(theta6)*(cos(theta5)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) + sin(theta2)*sin(theta3)*sin(theta5)))/125 + (11*cos(theta6)*(cos(theta2)*cos(theta4) + cos(theta3)*sin(theta2)*sin(theta4)))/125 + (107*sin(theta6)*(cos(theta2)*cos(theta4) + cos(theta3)*sin(theta2)*sin(theta4)))/1000, 0];
B= [ 0, - sin(theta7)*(sin(theta5)*(sin(theta2)*sin(theta4) + cos(theta2)*cos(theta3)*cos(theta4)) + cos(theta2)*cos(theta5)*sin(theta3)) - cos(theta7)*(cos(theta6)*(cos(theta5)*(sin(theta2)*sin(theta4) + cos(theta2)*cos(theta3)*cos(theta4)) - cos(theta2)*sin(theta3)*sin(theta5)) + sin(theta6)*(cos(theta4)*sin(theta2) - cos(theta2)*cos(theta3)*sin(theta4))),   cos(theta7)*(cos(theta6)*(cos(theta3)*sin(theta2)*sin(theta5) + cos(theta4)*cos(theta5)*sin(theta2)*sin(theta3)) - sin(theta2)*sin(theta3)*sin(theta4)*sin(theta6)) - sin(theta7)*(cos(theta3)*cos(theta5)*sin(theta2) - cos(theta4)*sin(theta2)*sin(theta3)*sin(theta5)), sin(theta5)*sin(theta7)*(cos(theta2)*cos(theta4) + cos(theta3)*sin(theta2)*sin(theta4)) - cos(theta7)*(sin(theta6)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) - cos(theta5)*cos(theta6)*(cos(theta2)*cos(theta4) + cos(theta3)*sin(theta2)*sin(theta4))), sin(theta7)*(cos(theta5)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) + sin(theta2)*sin(theta3)*sin(theta5)) - cos(theta6)*cos(theta7)*(sin(theta5)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) - cos(theta5)*sin(theta2)*sin(theta3)), -cos(theta7)*(sin(theta6)*(cos(theta5)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) + sin(theta2)*sin(theta3)*sin(theta5)) - cos(theta6)*(cos(theta2)*cos(theta4) + cos(theta3)*sin(theta2)*sin(theta4))),   cos(theta7)*(sin(theta5)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) - cos(theta5)*sin(theta2)*sin(theta3)) - sin(theta7)*(cos(theta6)*(cos(theta5)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) + sin(theta2)*sin(theta3)*sin(theta5)) + sin(theta6)*(cos(theta2)*cos(theta4) + cos(theta3)*sin(theta2)*sin(theta4)));...
0,   sin(theta7)*(cos(theta6)*(cos(theta5)*(sin(theta2)*sin(theta4) + cos(theta2)*cos(theta3)*cos(theta4)) - cos(theta2)*sin(theta3)*sin(theta5)) + sin(theta6)*(cos(theta4)*sin(theta2) - cos(theta2)*cos(theta3)*sin(theta4))) - cos(theta7)*(sin(theta5)*(sin(theta2)*sin(theta4) + cos(theta2)*cos(theta3)*cos(theta4)) + cos(theta2)*cos(theta5)*sin(theta3)), - sin(theta7)*(cos(theta6)*(cos(theta3)*sin(theta2)*sin(theta5) + cos(theta4)*cos(theta5)*sin(theta2)*sin(theta3)) - sin(theta2)*sin(theta3)*sin(theta4)*sin(theta6)) - cos(theta7)*(cos(theta3)*cos(theta5)*sin(theta2) - cos(theta4)*sin(theta2)*sin(theta3)*sin(theta5)), sin(theta7)*(sin(theta6)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) - cos(theta5)*cos(theta6)*(cos(theta2)*cos(theta4) + cos(theta3)*sin(theta2)*sin(theta4))) + cos(theta7)*sin(theta5)*(cos(theta2)*cos(theta4) + cos(theta3)*sin(theta2)*sin(theta4)), cos(theta7)*(cos(theta5)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) + sin(theta2)*sin(theta3)*sin(theta5)) + cos(theta6)*sin(theta7)*(sin(theta5)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) - cos(theta5)*sin(theta2)*sin(theta3)),  sin(theta7)*(sin(theta6)*(cos(theta5)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) + sin(theta2)*sin(theta3)*sin(theta5)) - cos(theta6)*(cos(theta2)*cos(theta4) + cos(theta3)*sin(theta2)*sin(theta4))), - sin(theta7)*(sin(theta5)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) - cos(theta5)*sin(theta2)*sin(theta3)) - cos(theta7)*(cos(theta6)*(cos(theta5)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) + sin(theta2)*sin(theta3)*sin(theta5)) + sin(theta6)*(cos(theta2)*cos(theta4) + cos(theta3)*sin(theta2)*sin(theta4)))];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T1=...
[ cos(theta1), -sin(theta1), 0,        0;...
 sin(theta1),  cos(theta1), 0,        0;...
      0,       0, 1, 333/1000;...
      0,       0, 0,        1];
  
T2=...
[cos(theta2),  -sin(theta2), 0, 0;...
      0,       0, 1, 0;...
-sin(theta2), -cos(theta2), 0, 0;...
      0,       0, 0, 1];
   T3=...
 [  cos(theta3),  -sin(theta3),0, 0;...
 0, 0, -1, -79/250;...
 sin(theta3), cos(theta3), 0, 0;...
  0, 0, 0, 1];
        T4=...
[cos(theta4), -sin(theta4),0, 33/400;...
 0, 0,  -1,      0;...
sin(theta4),  cos(theta4), 0,   0;...
0,0,0,      1];
T5=...
[  cos(theta5),  -sin(theta5), 0,-33/400;...
 0, 0,  1, 48/125;...
  -sin(theta5), -cos(theta5), 0, 0;...
  0, 0,0, 1];
T6=...
[  cos(theta6),-sin(theta6), 0, 0;...
 0, 0,-1, 0;...
 sin(theta6), cos(theta6), 0, 0;...
0, 0, 0, 1];
      
      T7=...
[ cos(theta7), -sin(theta7), 0, 11/125;...
 0, 0, -1,0;...
sin(theta7),  cos(theta7), 0, 0;...
0, 0,0,1];
      T8=...   
         [1,         0,         0,         0;...
         0,    1,         0,         0;...
         0,         0,    1,    0.1070;...
         0,         0,         0,    1];
 T1_0=...
[ cos(theta1_0), -sin(theta1_0), 0,        0;...
 sin(theta1_0),  cos(theta1_0), 0,        0;...
      0,       0, 1, 333/1000;...
      0,       0, 0,        1];
  
T2_0=...
[cos(theta2_0),  -sin(theta2_0), 0, 0;...
      0,       0, 1, 0;...
-sin(theta2_0), -cos(theta2_0), 0, 0;...
      0,       0, 0, 1];
   T3_0=...
 [  cos(theta3_0),  -sin(theta3_0),0, 0;...
 0, 0, -1, -79/250;...
 sin(theta3_0), cos(theta3_0), 0, 0;...
  0, 0, 0, 1];
        T4_0=...
[cos(theta4_0), -sin(theta4_0),0, 33/400;...
 0, 0,  -1,      0;...
sin(theta4_0),  cos(theta4_0), 0,   0;...
0,0,0,      1];
T5_0=...
[  cos(theta5_0),  -sin(theta5_0), 0,-33/400;...
 0, 0,  1, 48/125;...
  -sin(theta5_0), -cos(theta5_0), 0, 0;...
  0, 0,0, 1];
T6_0=...
[  cos(theta6_0),-sin(theta6_0), 0, 0;...
 0, 0,-1, 0;...
 sin(theta6_0), cos(theta6_0), 0, 0;...
0, 0, 0, 1];
      
      T7_0=...
[ cos(theta7_0), -sin(theta7_0), 0, 11/125;...
 0, 0, -1,0;...
sin(theta7_0),  cos(theta7_0), 0, 0;...
0, 0,0,1];
      T8_0=...   
         [1,         0,         0,         0;...
         0,    1,         0,         0;...
         0,         0,    1,    0.1070;...
         0,         0,         0,    1];
       
Basepos=[0;0;0;1]+basea(i,:)';
Pos8=T1*T2*T3*T4*T5*T6*T7*T8*Basepos;
myPosition=Pos8(1:3);
myPostionstore=[myPostionstore;myPosition];  
myJ(:,:,i)=A;
mydotJ(:,:,i)=dotA;
myJw(:,:,i)=B;
mydotJw(:,:,i)=dotB;
myww(:,:,i)=T1*T2*T3*T4*T5*T6*T7*T8;
myw((2*(i-1)+1):(2*i),1)=[myww(3,1,i);myww(3,2,i)];
myww_0(:,:,i)=T1_0*T2_0*T3_0*T4_0*T5_0*T6_0*T7_0*T8_0;
myw_0((2*(i-1)+1):(2*i),1)=[myww_0(3,1,i);myww_0(3,2,i)];
J((3*(i-1)+1):(3*i),(7*(i-1)+1):(7*i))=myJ(:,:,i);
dotJ((3*(i-1)+1):(3*i),(7*(i-1)+1):(7*i))=mydotJ(:,:,i);
Jw((2*(i-1)+1):(2*i),(7*(i-1)+1):(7*i))=myJw(:,:,i);
dotJw((2*(i-1)+1):(2*i),(7*(i-1)+1):(7*i))=mydotJw(:,:,i);
 
       

end

H(:,:,m)=[    eye(7*n,7*n),     zeros(7*n,3*n),    J',             zeros(7*n,3*n),    Jw';
              zeros(3*n,7*n),   zeros(3*n,3*n),    zeros(3*n,3*n), eye(3*n,3*n),      zeros(3*n,2*n);
              J,                zeros(3*n,3*n),    zeros(3*n,3*n), zeros(3*n,3*n),    zeros(3*n,2*n);
              zeros(3*n,7*n),   eye(3*n,3*n),      zeros(3*n,3*n), zeros(3*n,3*n),    zeros(3*n,2*n);             
              Jw,               zeros(2*n,3*n),    zeros(2*n,3*n), zeros(2*n,3*n),    zeros(2*n,2*n)];

bjl(:,:,m)=-c3*kron(L,eye(3))*myPostionstore+c3*kron(L',eye(3))*mw(:,:,m)+c3*kron(Lambda,eye(3))*(kron(ones(n,1),rd)-myPostionstore);
cj1(:,:,m)=-100*(myw-myw_0);
ujl(:,:,m)=[-c1*(myTheta(:,:,m)-myTheta(:,:,1));zeros(3*n,1);bjl(:,:,m);-kron(L,eye(3))*myPostionstore;cj1(:,:,m)];
err(:,:,m)=H(:,:,m)*ddTheta(:,:,m)-ujl(:,:,m);
sumee(:,:,m+1)=err(:,:,m)+sumee(:,:,m);    
h_H(:,:,m)=H(:,:,m)'*pinv(H(:,:,m)*H(:,:,m)');
noise=0.5*m*dt+0.1;
error_end(:,:,m)=myw-myw_0;





if m==1 
    ddTheta(:,:,m+1)=ddTheta(:,:,m)-h_H(:,:,m)*((H(:,:,m)*ddTheta(:,:,m)-ujl(:,:,m))+noise);
else
ddTheta(:,:,m+1)=ddTheta(:,:,m)+h_H(:,:,m)*(-(H(:,:,m)-H(:,:,m-1))*ddTheta(:,:,m)-dt*c1*err(:,:,m)+(ujl(:,:,m)-ujl(:,:,m-1))-dt*dt*c2*sumee(:,:,m)+noise);

end

myDTheta(:,:,m+1)=ddTheta(1:7*n,:,m+1);
myDw(:,:,m+1)=ddTheta(7*n+1:10*n,:,m+1);
lambda1(:,:,m+1)=ddTheta(10*n+1:13*n,:,m+1);
lambda2(:,:,m+1)=ddTheta(13*n+1:16*n,:,m+1);
lambda3(:,:,m+1)=ddTheta(16*n+1:18*n,:,m+1);
mw(:,:,m+1)=mw(:,:,m)-dt*(kron(L,eye(3))*myPostionstore);
myTheta(:,:,m+1)=myTheta(:,:,m)+dt*myDTheta(:,:,m+1);
myDDTheta(:,:,m+1)=myDTheta(:,:,m+1)-myDTheta(:,:,m);

Y011(1:7*n,:,m+1)=myTheta(:,:,m+1);
Y011(7*n+1:10*n,:,m+1)=mw(:,:,m+1);
Y011(10*n+1:17*n,:,m+1)=myDTheta(:,:,m+1);
Y011(17*n+1:20*n,:,m+1)=myDw(:,:,m+1);
Y011(20*n+1:23*n,:,m+1)=lambda1(:,:,m+1);
Y011(23*n+1:26*n,:,m+1)=lambda2(:,:,m+1);
Y011(26*n+1:28*n,:,m+1)=lambda3(:,:,m+1);
YY=Y011(:,:,m);
Y01=[Y01,YY];



theta0=y0jl;
d_theta(:,:,m)=Y011(10*n+1:17*n,:,m+1);

g(:,:,m+1) = get_GravityVector(theta0);
c(:,:,m+1) = get_CoriolisVector(theta0,d_theta(:,:,m));
M(:,:,m+1)=get_MassMatrix(theta0);

MMM(:,:,m+1)=[M(:,:,m+1),zeros(7,7),zeros(7,7),zeros(7,7),zeros(7,7),zeros(7,7);
            zeros(7,7),M(:,:,m+1),zeros(7,7),zeros(7,7),zeros(7,7),zeros(7,7);
            zeros(7,7),zeros(7,7),M(:,:,m+1),zeros(7,7),zeros(7,7),zeros(7,7);
            zeros(7,7),zeros(7,7),zeros(7,7),M(:,:,m+1),zeros(7,7),zeros(7,7);
            zeros(7,7),zeros(7,7),zeros(7,7),zeros(7,7),M(:,:,m+1),zeros(7,7);
            zeros(7,7),zeros(7,7),zeros(7,7),zeros(7,7),zeros(7,7),M(:,:,m+1)];
            
            
cc(:,:,m+1)=[c(:,:,m+1);c(:,:,m+1);c(:,:,m+1);c(:,:,m+1);c(:,:,m+1);c(:,:,m+1)];
gg(:,:,m+1)=[g(:,:,m+1);g(:,:,m+1);g(:,:,m+1);g(:,:,m+1);g(:,:,m+1);g(:,:,m+1)];
tau(:,:,m+1)=MMM(:,:,m+1)*myDDTheta(:,:,m+1)+ cc(:,:,m+1) + gg(:,:,m+1);
TAU=tau(:,:,m+1);
tau_store=[tau_store,TAU];
error_E=error_end(:,:,m);
error_Estore=[error_Estore,error_E];









 end

save Y01;
save tau_store;
save error_Estore；


myvisualization12;  
