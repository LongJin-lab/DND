global w;
global l;
global d;
global basea;
global n;
x00= [1 2 3];
y00=[1 3];
[x001,y001]=meshgrid(x00,y00);

all_error_epplot = zeros(0,0);
all_speed = zeros(0,0);

for i=1:n  
RobotPos=zeros(4,1);
RobotPos(1)=x001(i);
RobotPos(2)=y001(i);
myPostionstore=zeros(3,0);
Pos1store=zeros(4,0);
Pos2store=zeros(4,0);
Pos3store=zeros(4,0);
Pos4store=zeros(4,0);
Pos5store=zeros(4,0);
Pos6store=zeros(4,0);
Pos7store=zeros(4,0);
Pos8store=zeros(4,0);
Pos8store1=zeros(4,0);
Pos8store2=zeros(4,0);
Pos8store3=zeros(4,0);
Pos8store4=zeros(4,0);
Pos8store5=zeros(4,0);
Pos8store6=zeros(4,0);
vstore=zeros(0,3); 
for m=1:length(t_store)-1
    theta1=Y01((7*(i-1)+1),m);
    theta2=Y01((7*(i-1)+2),m);
    theta3=Y01((7*(i-1)+3),m);
    theta4=Y01((7*(i-1)+4),m);
    theta5=Y01((7*(i-1)+5),m);
    theta6=Y01((7*(i-1)+6),m);
    theta7=Y01((7*(i-1)+7),m);
    
        dtheta1=Y01(10*n+(7*(i-1)+1),m);
        dtheta2=Y01(10*n+(7*(i-1)+2),m);
        dtheta3=Y01(10*n+(7*(i-1)+3),m);
        dtheta4=Y01(10*n+(7*(i-1)+4),m);
        dtheta5=Y01(10*n+(7*(i-1)+5),m);
        dtheta6=Y01(10*n+(7*(i-1)+6),m);
        dtheta7=Y01(10*n+(7*(i-1)+7),m);
    

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

 
Basepos=[0;0;0;1]+basea(i,:)';
Pos1=T1*Basepos+RobotPos;
Pos2=T1*T2*Basepos+RobotPos;
Pos3=T1*T2*T3*Basepos+RobotPos;
Pos4=T1*T2*T3*T4*Basepos+RobotPos;
Pos5=T1*T2*T3*T4*T5*Basepos+RobotPos;
Pos6=T1*T2*T3*T4*T5*T6*Basepos+RobotPos;
Pos7=T1*T2*T3*T4*T5*T6*T7*Basepos+RobotPos;
Pos8=T1*T2*T3*T4*T5*T6*T7*T8*Basepos+RobotPos;
Pos81=T1*T2*T3*T4*T5*T6*T7*T8*([0;0;0;1]+basea(1,:)')+RobotPos;
Pos82=T1*T2*T3*T4*T5*T6*T7*T8*([0;0;0;1]+basea(2,:)')+RobotPos;
Pos83=T1*T2*T3*T4*T5*T6*T7*T8*([0;0;0;1]+basea(3,:)')+RobotPos;
Pos84=T1*T2*T3*T4*T5*T6*T7*T8*([0;0;0;1]+basea(4,:)')+RobotPos;
Pos85=T1*T2*T3*T4*T5*T6*T7*T8*([0;0;0;1]+basea(5,:)')+RobotPos;
Pos86=T1*T2*T3*T4*T5*T6*T7*T8*([0;0;0;1]+basea(6,:)')+RobotPos;
myPosition=Pos8(1:3);


Pos1store= [Pos1store,Pos1];
Pos2store= [Pos2store,Pos2]; 
Pos3store= [Pos3store,Pos3]; 
Pos4store= [Pos4store,Pos4]; 
Pos5store= [Pos5store,Pos5]; 
Pos6store= [Pos6store,Pos6]; 
Pos7store= [Pos7store,Pos7]; 
Pos8store= [Pos8store,Pos8]; 
Pos8store1= [Pos8store1,Pos81]; 
Pos8store2= [Pos8store2,Pos82]; 
Pos8store3= [Pos8store3,Pos83]; 
Pos8store4= [Pos8store4,Pos84]; 
Pos8store5= [Pos8store5,Pos85]; 
Pos8store6= [Pos8store6,Pos86]; 
myPostionstore= [myPostionstore,myPosition];   


A =...
[   (11*sin(theta6)*(sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) - cos(theta4)*sin(theta1)*sin(theta2)))/125 - (79*sin(theta1)*sin(theta2))/250 - (107*cos(theta6)*(sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) - cos(theta4)*sin(theta1)*sin(theta2)))/1000 - (33*cos(theta1)*sin(theta3))/400 + (33*cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)))/400 - (11*cos(theta6)*(cos(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta1)*sin(theta2)*sin(theta4)) + sin(theta5)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))))/125 + (48*sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)))/125 - (107*sin(theta6)*(cos(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta1)*sin(theta2)*sin(theta4)) + sin(theta5)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))))/1000 - (48*cos(theta4)*sin(theta1)*sin(theta2))/125 + (33*sin(theta1)*sin(theta2)*sin(theta4))/400 - (33*cos(theta2)*cos(theta3)*sin(theta1))/400, (79*cos(theta1)*cos(theta2))/250 - (107*cos(theta6)*(cos(theta1)*cos(theta2)*cos(theta4) + cos(theta1)*cos(theta3)*sin(theta2)*sin(theta4)))/1000 + (11*sin(theta6)*(cos(theta1)*cos(theta2)*cos(theta4) + cos(theta1)*cos(theta3)*sin(theta2)*sin(theta4)))/125 + (11*cos(theta6)*(cos(theta5)*(cos(theta1)*cos(theta2)*sin(theta4) - cos(theta1)*cos(theta3)*cos(theta4)*sin(theta2)) + cos(theta1)*sin(theta2)*sin(theta3)*sin(theta5)))/125 + (107*sin(theta6)*(cos(theta5)*(cos(theta1)*cos(theta2)*sin(theta4) - cos(theta1)*cos(theta3)*cos(theta4)*sin(theta2)) + cos(theta1)*sin(theta2)*sin(theta3)*sin(theta5)))/1000 + (48*cos(theta1)*cos(theta2)*cos(theta4))/125 - (33*cos(theta1)*cos(theta3)*sin(theta2))/400 - (33*cos(theta1)*cos(theta2)*sin(theta4))/400 + (33*cos(theta1)*cos(theta3)*cos(theta4)*sin(theta2))/400 + (48*cos(theta1)*cos(theta3)*sin(theta2)*sin(theta4))/125, (33*cos(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)))/400 - (33*cos(theta3)*sin(theta1))/400 + (48*sin(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)))/125 + (11*cos(theta6)*(sin(theta5)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta4)*cos(theta5)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))))/125 + (107*sin(theta6)*(sin(theta5)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta4)*cos(theta5)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))))/1000 - (107*cos(theta6)*sin(theta4)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)))/1000 + (11*sin(theta4)*sin(theta6)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3)))/125 - (33*cos(theta1)*cos(theta2)*sin(theta3))/400, (11*sin(theta6)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta1)*sin(theta2)*sin(theta4)))/125 - (107*cos(theta6)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta1)*sin(theta2)*sin(theta4)))/1000 + (48*cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))/125 - (33*sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))/400 - (48*cos(theta1)*sin(theta2)*sin(theta4))/125 + (11*cos(theta5)*cos(theta6)*(sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + cos(theta1)*cos(theta4)*sin(theta2)))/125 + (107*cos(theta5)*sin(theta6)*(sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + cos(theta1)*cos(theta4)*sin(theta2)))/1000 - (33*cos(theta1)*cos(theta4)*sin(theta2))/400,   (11*cos(theta6)*(sin(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta1)*sin(theta2)*sin(theta4)) - cos(theta5)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))))/125 + (107*sin(theta6)*(sin(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta1)*sin(theta2)*sin(theta4)) - cos(theta5)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))))/1000, (11*cos(theta6)*(sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + cos(theta1)*cos(theta4)*sin(theta2)))/125 + (107*sin(theta6)*(sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + cos(theta1)*cos(theta4)*sin(theta2)))/1000 - (107*cos(theta6)*(cos(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta1)*sin(theta2)*sin(theta4)) + sin(theta5)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))))/1000 + (11*sin(theta6)*(cos(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta1)*sin(theta2)*sin(theta4)) + sin(theta5)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))))/125, 0;...
 (79*cos(theta1)*sin(theta2))/250 - (33*sin(theta1)*sin(theta3))/400 - (107*cos(theta6)*(sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + cos(theta1)*cos(theta4)*sin(theta2)))/1000 + (11*sin(theta6)*(sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + cos(theta1)*cos(theta4)*sin(theta2)))/125 - (11*cos(theta6)*(cos(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta1)*sin(theta2)*sin(theta4)) + sin(theta5)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))))/125 + (33*cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))/400 - (107*sin(theta6)*(cos(theta5)*(cos(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) - cos(theta1)*sin(theta2)*sin(theta4)) + sin(theta5)*(cos(theta3)*sin(theta1) + cos(theta1)*cos(theta2)*sin(theta3))))/1000 + (48*sin(theta4)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))/125 - (33*cos(theta1)*sin(theta2)*sin(theta4))/400 + (33*cos(theta1)*cos(theta2)*cos(theta3))/400 + (48*cos(theta1)*cos(theta4)*sin(theta2))/125, (11*cos(theta6)*(cos(theta5)*(cos(theta2)*sin(theta1)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta1)*sin(theta2)) + sin(theta1)*sin(theta2)*sin(theta3)*sin(theta5)))/125 + (79*cos(theta2)*sin(theta1))/250 + (107*sin(theta6)*(cos(theta5)*(cos(theta2)*sin(theta1)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta1)*sin(theta2)) + sin(theta1)*sin(theta2)*sin(theta3)*sin(theta5)))/1000 - (107*cos(theta6)*(cos(theta2)*cos(theta4)*sin(theta1) + cos(theta3)*sin(theta1)*sin(theta2)*sin(theta4)))/1000 + (11*sin(theta6)*(cos(theta2)*cos(theta4)*sin(theta1) + cos(theta3)*sin(theta1)*sin(theta2)*sin(theta4)))/125 - (33*cos(theta3)*sin(theta1)*sin(theta2))/400 - (33*cos(theta2)*sin(theta1)*sin(theta4))/400 + (48*cos(theta2)*cos(theta4)*sin(theta1))/125 + (33*cos(theta3)*cos(theta4)*sin(theta1)*sin(theta2))/400 + (48*cos(theta3)*sin(theta1)*sin(theta2)*sin(theta4))/125, (33*cos(theta1)*cos(theta3))/400 - (33*cos(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)))/400 - (48*sin(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)))/125 - (11*cos(theta6)*(sin(theta5)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) - cos(theta4)*cos(theta5)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))))/125 - (107*sin(theta6)*(sin(theta5)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) - cos(theta4)*cos(theta5)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))))/1000 - (33*cos(theta2)*sin(theta1)*sin(theta3))/400 + (107*cos(theta6)*sin(theta4)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)))/1000 - (11*sin(theta4)*sin(theta6)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3)))/125, (107*cos(theta6)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta1)*sin(theta2)*sin(theta4)))/1000 - (11*sin(theta6)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta1)*sin(theta2)*sin(theta4)))/125 - (48*cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)))/125 + (33*sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)))/400 - (33*cos(theta4)*sin(theta1)*sin(theta2))/400 - (48*sin(theta1)*sin(theta2)*sin(theta4))/125 - (11*cos(theta5)*cos(theta6)*(sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) - cos(theta4)*sin(theta1)*sin(theta2)))/125 - (107*cos(theta5)*sin(theta6)*(sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) - cos(theta4)*sin(theta1)*sin(theta2)))/1000, - (11*cos(theta6)*(sin(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta1)*sin(theta2)*sin(theta4)) - cos(theta5)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))))/125 - (107*sin(theta6)*(sin(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta1)*sin(theta2)*sin(theta4)) - cos(theta5)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))))/1000, (107*cos(theta6)*(cos(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta1)*sin(theta2)*sin(theta4)) + sin(theta5)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))))/1000 - (107*sin(theta6)*(sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) - cos(theta4)*sin(theta1)*sin(theta2)))/1000 - (11*cos(theta6)*(sin(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) - cos(theta4)*sin(theta1)*sin(theta2)))/125 - (11*sin(theta6)*(cos(theta5)*(cos(theta4)*(cos(theta1)*sin(theta3) + cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta1)*sin(theta2)*sin(theta4)) + sin(theta5)*(cos(theta1)*cos(theta3) - cos(theta2)*sin(theta1)*sin(theta3))))/125, 0;...
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                 (33*sin(theta2)*sin(theta4))/400 - (33*cos(theta2)*cos(theta3))/400 - (48*cos(theta4)*sin(theta2))/125 - (79*sin(theta2))/250 - (11*cos(theta6)*(cos(theta5)*(sin(theta2)*sin(theta4) + cos(theta2)*cos(theta3)*cos(theta4)) - cos(theta2)*sin(theta3)*sin(theta5)))/125 - (107*sin(theta6)*(cos(theta5)*(sin(theta2)*sin(theta4) + cos(theta2)*cos(theta3)*cos(theta4)) - cos(theta2)*sin(theta3)*sin(theta5)))/1000 + (107*cos(theta6)*(cos(theta4)*sin(theta2) - cos(theta2)*cos(theta3)*sin(theta4)))/1000 - (11*sin(theta6)*(cos(theta4)*sin(theta2) - cos(theta2)*cos(theta3)*sin(theta4)))/125 + (33*cos(theta2)*cos(theta3)*cos(theta4))/400 + (48*cos(theta2)*cos(theta3)*sin(theta4))/125,                                                                                                                                                                                                                                                                                                                                                                                (33*sin(theta2)*sin(theta3))/400 + (11*cos(theta6)*(cos(theta3)*sin(theta2)*sin(theta5) + cos(theta4)*cos(theta5)*sin(theta2)*sin(theta3)))/125 + (107*sin(theta6)*(cos(theta3)*sin(theta2)*sin(theta5) + cos(theta4)*cos(theta5)*sin(theta2)*sin(theta3)))/1000 - (33*cos(theta4)*sin(theta2)*sin(theta3))/400 - (48*sin(theta2)*sin(theta3)*sin(theta4))/125 + (107*cos(theta6)*sin(theta2)*sin(theta3)*sin(theta4))/1000 - (11*sin(theta2)*sin(theta3)*sin(theta4)*sin(theta6))/125,                                                                                                                                                                                                                                                                                                                         (107*cos(theta6)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)))/1000 - (48*cos(theta2)*sin(theta4))/125 - (33*cos(theta2)*cos(theta4))/400 - (11*sin(theta6)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)))/125 - (33*cos(theta3)*sin(theta2)*sin(theta4))/400 + (11*cos(theta5)*cos(theta6)*(cos(theta2)*cos(theta4) + cos(theta3)*sin(theta2)*sin(theta4)))/125 + (107*cos(theta5)*sin(theta6)*(cos(theta2)*cos(theta4) + cos(theta3)*sin(theta2)*sin(theta4)))/1000 + (48*cos(theta3)*cos(theta4)*sin(theta2))/125,                                                                                                                                                                                         - (11*cos(theta6)*(sin(theta5)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) - cos(theta5)*sin(theta2)*sin(theta3)))/125 - (107*sin(theta6)*(sin(theta5)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) - cos(theta5)*sin(theta2)*sin(theta3)))/1000,                                                                                                                                                                                                                                                                                                 (107*cos(theta6)*(cos(theta5)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) + sin(theta2)*sin(theta3)*sin(theta5)))/1000 - (11*sin(theta6)*(cos(theta5)*(cos(theta2)*sin(theta4) - cos(theta3)*cos(theta4)*sin(theta2)) + sin(theta2)*sin(theta3)*sin(theta5)))/125 + (11*cos(theta6)*(cos(theta2)*cos(theta4) + cos(theta3)*sin(theta2)*sin(theta4)))/125 + (107*sin(theta6)*(cos(theta2)*cos(theta4) + cos(theta3)*sin(theta2)*sin(theta4)))/1000, 0];
%%%%%%%%%%%%

v_workspace=A*[dtheta1;dtheta2;dtheta3;dtheta4;dtheta5;dtheta6;dtheta7];  %
vstore=[vstore,v_workspace];

end

 
 figure(4),
 hold on,
 plot( myPostionstore(1,:), myPostionstore(2,:)) 
 


x00= [1 2.5];
y00=[1 3];

axis equal



    
    
TT=t_store(1,1:end-1);

T=2*pi;
a=4;
b=1;
rx=0.005*(6*(cos(pi/5*a*TT/T*10-pi/2))+5*cos(pi/5*b*TT/T*10-pi/2));
ry=0.005*(-6*(sin(pi/5*a*TT/T*10-pi/2))+5*sin(pi/5*b*TT/T*10-pi/2)-1);
rz=0.785104;

pd=[rx+ones(1,length(TT))*x001(i);ry++ones(1,length(TT))*y001(i);ones(1,length(TT))*0.785104];



 


errorP=pd-myPostionstore;

numberS=num2str(i);
peS=['position error',' ',numberS];



figure(5),title('manipulator 1,2,3')
hold on,
grid on;
for kk=0:0.045:0.5
 plot3([RobotPos(1),interp1(t_store(1:length(t_store)-1),Pos1store(1,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos2store(1,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos3store(1,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos4store(1,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos5store(1,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos6store(1,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos7store(1,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos8store(1,:),kk*ts,'linear')],...
       [RobotPos(2),interp1(t_store(1:length(t_store)-1),Pos1store(2,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos2store(2,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos3store(2,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos4store(2,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos5store(2,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos6store(2,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos7store(2,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos8store(2,:),kk*ts,'linear')],...
       [RobotPos(3),interp1(t_store(1:length(t_store)-1),Pos1store(3,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos2store(3,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos3store(3,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos4store(3,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos5store(3,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos6store(3,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos7store(3,:),kk*ts,'linear'),interp1(t_store(1:length(t_store)-1),Pos8store(3,:),kk*ts,'linear')]);
 axis equal
end 
 
 path1=plot3(Pos8store1(1,:),Pos8store1(2,:),Pos8store1(3,:),'color',[117/256 209/256 135/256]);
 set(path1,'linewidth',1); hold on;
 path2=plot3(Pos8store2(1,:),Pos8store2(2,:),Pos8store2(3,:),'color',[83/256 190/256 153/256]);
 set(path2,'linewidth',1); hold on;
 path3=plot3(Pos8store3(1,:),Pos8store3(2,:),Pos8store3(3,:),'color',[89/256 137/256 165/256]);
 set(path3,'linewidth',1); hold on;
 path4=plot3(Pos8store4(1,:),Pos8store4(2,:),Pos8store4(3,:),'color',[239/256 121/256 126/256]);%plot the  red trajectory
 set(path4,'linewidth',1); hold on;
 path5=plot3(Pos8store5(1,:),Pos8store5(2,:),Pos8store5(3,:),'color',[253/256 204/256 165/256]);
 set(path5,'linewidth',1); hold on;
 path6=plot3(Pos8store6(1,:),Pos8store6(2,:),Pos8store6(3,:),'color',[196/256 93/256 147/256]);%plot the  red trajectory
 set(path6,'linewidth',1); hold on;    
 legend('robot1','robot2','robot3','robot4','robot5','robot6');
 
 axis equal
  view(3)

  
speed_i=(vstore(1,:).^2+vstore(2,:).^2+vstore(3,:).^2).^(0.5);
all_speed = [all_speed, speed_i'];

error_i=(errorP(1,:).^2+errorP(2,:).^2+errorP(3,:).^2).^(0.5);
all_error_epplot = [all_error_epplot, error_i'];

  
  
end


figure(100),clf(100)
times1 = plot(t_store(1:length(t_store)-1),all_error_epplot(:,1),'color',[252/256 171/256 143/256]);
set(times1,'linewidth',2,'LineStyle','-'); hold on;
times2 = plot(t_store(1:length(t_store)-1),all_error_epplot(:,2),'color',[226/256 85/256 8/256]);
set(times2,'linewidth',2,'LineStyle','--'); hold on;
times3 = plot(t_store(1:length(t_store)-1),all_error_epplot(:,3),'color',[46/256 151/256 78/256]);
set(times3,'linewidth',2,'LineStyle','-.'); hold on;
times4 = plot(t_store(1:length(t_store)-1),all_error_epplot(:,4),'color',[46/256 126/256 187/256]);
set(times4,'linewidth',2,'LineStyle',':'); hold on;
times5 = plot(t_store(1:length(t_store)-1),all_error_epplot(:,5),'color',[114/256 98/256 172/256]);
set(times5,'linewidth',1,'LineStyle',':'); hold on;
times6 = plot(t_store(1:length(t_store)-1),all_error_epplot(:,6),'color',[95/256 95/256 95/256]);
set(times6,'linewidth',1,'LineStyle','-.'); hold on;
legend('robot1','robot2','robot3','robot4','robot5','robot6'); 
title('error')
figure(101),clf(101)
times1 = plot(t_store(1:length(t_store)-1),all_speed(:,1),'color',[252/256 171/256 143/256]);
set(times1,'linewidth',2,'LineStyle','-'); hold on;
times2 = plot(t_store(1:length(t_store)-1),all_speed(:,2),'color',[226/256 85/256 8/256]);
set(times2,'linewidth',2,'LineStyle','--'); hold on;
times3 = plot(t_store(1:length(t_store)-1),all_speed(:,3),'color',[46/256 151/256 78/256]);
set(times3,'linewidth',2,'LineStyle','-.'); hold on;
times4 = plot(t_store(1:length(t_store)-1),all_speed(:,4),'color',[46/256 126/256 187/256]);
set(times4,'linewidth',2,'LineStyle',':'); hold on;
times5 = plot(t_store(1:length(t_store)-1),all_speed(:,5),'color',[114/256 98/256 172/256]);
set(times5,'linewidth',1,'LineStyle',':'); hold on;
times6 = plot(t_store(1:length(t_store)-1),all_speed(:,6),'color',[95/256 95/256 95/256]);
set(times6,'linewidth',1,'LineStyle','-.'); hold on;
legend('robot1','robot2','robot3','robot4','robot5','robot6'); 
title('speed')



figure(8),clf(8), 
hold on;
aplot=plot(t_store(1:length(t_store)-1),Y01(1:7,:),'-','color',[252/256 171/256 143/256],'linewidth',2);
aplot2=plot(t_store(1:length(t_store)-1),Y01(8:14,:),'--','color',[226/256 85/256 8/256],'linewidth',2);
aplot3=plot(t_store(1:length(t_store)-1),Y01(15:21,:),'-.','color',[46/256 151/256 78/256],'linewidth',2);
aplot4=plot(t_store(1:length(t_store)-1),Y01(22:28,:),':','color',[46/256 126/256 187/256], 'linewidth',2);
aplot5=plot(t_store(1:length(t_store)-1),Y01(29:35,:), '--','color',[114/256 98/256 172/256], 'linewidth',1);
aplot6=plot(t_store(1:length(t_store)-1),Y01(36:42,:), ':','color',[95/256 95/256 95/256], 'linewidth',1);
set(gca,'linewidth',1.2)
title('\theta')
legend([aplot(1),aplot2(1),aplot3(1),aplot4(1),aplot5(1),aplot6(1)],'robot1','robot2','robot3','robot4','robot5','robot6');

figure(9),clf(9),   % plot of dot_theta
hold on;
eepplot1=plot(t_store(1:length(t_store)-1),Y01(61:67,:),'-','color',[252/256 171/256 143/256],'linewidth',2);
eepplot2=plot(t_store(1:length(t_store)-1),Y01(68:74,:),'--','color',[226/256 85/256 8/256],'linewidth',2);
eepplot3=plot(t_store(1:length(t_store)-1),Y01(75:81,:),'-.','color',[46/256 151/256 78/256],'linewidth',2);
eepplot4=plot(t_store(1:length(t_store)-1),Y01(82:88,:),':','color',[46/256 126/256 187/256], 'linewidth',2);
eepplot5=plot(t_store(1:length(t_store)-1),Y01(89:95,:),'--','color',[114/256 98/256 172/256], 'linewidth',1);
eepplot6=plot(t_store(1:length(t_store)-1),Y01(96:102,:),':','color',[95/256 95/256 95/256], 'linewidth',1);
set(gca,'linewidth',1.2)
title('$\dot{\theta}$','interpreter','latex')
legend([eepplot1(1),eepplot2(1),eepplot3(1),eepplot4(1),eepplot5(1),eepplot6(1)],'robot1','robot2','robot3','robot4','robot5','robot6');
figure(10),clf(10),   
hold on;
eepplot1=plot(t_store(1:length(t_store)),myDDTheta(1:7,:),'-','color',[252/256 171/256 143/256],'linewidth',2);
eepplot2=plot(t_store(1:length(t_store)),myDDTheta(8:14,:),'--','color',[226/256 85/256 8/256],'linewidth',2);
eepplot3=plot(t_store(1:length(t_store)),myDDTheta(15:21,:),'-.','color',[46/256 151/256 78/256],'linewidth',2);
eepplot4=plot(t_store(1:length(t_store)),myDDTheta(22:28,:),':','color',[46/256 126/256 187/256], 'linewidth',2);
eepplot5=plot(t_store(1:length(t_store)),myDDTheta(29:35,:),'--','color',[114/256 98/256 172/256], 'linewidth',1);
eepplot6=plot(t_store(1:length(t_store)),myDDTheta(36:42,:),'-.','color',[95/256 95/256 95/256], 'linewidth',1);
set(gca,'linewidth',1.2)
title('$\ddot{\theta}$','interpreter','latex')
legend([eepplot1(1),eepplot2(1),eepplot3(1),eepplot4(1),eepplot5(1),eepplot6(1)],'robot1','robot2','robot3','robot4','robot5','robot6');


figure(11),clf(11),   % plot of \tau
hold on;

eepplot1=plot(t_store(1:length(t_store)-1),tau_store(1:7,:),'-','color',[252/256 171/256 143/256],'linewidth',2);
eepplot2=plot(t_store(1:length(t_store)-1),tau_store(8:14,:),'--','color',[226/256 85/256 8/256],'linewidth',2);
eepplot3=plot(t_store(1:length(t_store)-1),tau_store(15:21,:),'-.','color',[46/256 151/256 78/256],'linewidth',2);
eepplot4=plot(t_store(1:length(t_store)-1),tau_store(22:28,:),':','color',[46/256 126/256 187/256], 'linewidth',2);
eepplot5=plot(t_store(1:length(t_store)-1),tau_store(29:35,:),'--','color',[114/256 98/256 172/256], 'linewidth',1);
eepplot6=plot(t_store(1:length(t_store)-1),tau_store(36:42,:),'-.','color',[95/256 95/256 95/256], 'linewidth',1);
set(gca,'linewidth',1.2)
title('$\tau$','interpreter','latex')
legend([eepplot1(1),eepplot2(1),eepplot3(1),eepplot4(1),eepplot5(1),eepplot6(1)],'robot1','robot2','robot3','robot4','robot5','robot6');

figure(12),clf(12),  
hold on;

eepplot1=plot(t_store(1:length(t_store)-1),error_Estore(1:2,:),'-','color',[252/256 171/256 143/256],'linewidth',2);
eepplot2=plot(t_store(1:length(t_store)-1),error_Estore(3:4,:),'--','color',[226/256 85/256 8/256],'linewidth',2);
eepplot3=plot(t_store(1:length(t_store)-1),error_Estore(5:6,:),'-.','color',[46/256 151/256 78/256],'linewidth',2);
eepplot4=plot(t_store(1:length(t_store)-1),error_Estore(7:8,:),':','color',[46/256 126/256 187/256], 'linewidth',2);
eepplot5=plot(t_store(1:length(t_store)-1),error_Estore(9:10,:),'--','color',[114/256 98/256 172/256], 'linewidth',1);
eepplot6=plot(t_store(1:length(t_store)-1),error_Estore(11:12,:),'-.','color',[95/256 95/256 95/256], 'linewidth',1);
set(gca,'linewidth',1.2)
legend([eepplot1(1),eepplot2(1),eepplot3(1),eepplot4(1),eepplot5(1),eepplot6(1)],'robot1','robot2','robot3','robot4','robot5','robot6');



