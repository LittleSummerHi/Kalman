clear all;
close all;
clc;
% 初始值设定
k=500;
Var_num=6;
t=1:1:k;
T=1;
xita=0.05;
xitaR=10;
xitaaf=0.4;
xitabeta=0.5;
x=zeros(1,k);
x_daoshu=zeros(1,k);
y=zeros(1,k);
y_daoshu=zeros(1,k);
z=zeros(1,k);
z_daoshu=zeros(1,k);
Vx=zeros(1,k);
Vy=zeros(1,k);
Vz=zeros(1,k);
for i=1:k
    Vx(i)=3;
    Vy(i)=4;
    Vz(i)=5;
end
x0=10;
y0=20;
z0=30;
x(1)=x0;
y(1)=y0;
z(1)=z0;
x_daoshu(1)=Vx(1);
y_daoshu(1)=Vy(1);
z_daoshu(1)=Vz(1);
for n=2:k
    x(n)=x0+Vx(n)*n*T;
    x_daoshu(n)=(x(n)-x(n-1))/T;
    y(n)=y0+Vy(n)*n*T;
    y_daoshu(n)=(y(n)-y(n-1))/T;
    z(n)=z0+Vz(n)*n*T;
    z_daoshu(n)=(z(n)-z(n-1))/T;
end
System_Noise=randn(1,k);
x=x+System_Noise;
x_daoshu=x_daoshu+System_Noise;
y=y+System_Noise;
y_daoshu=y_daoshu+System_Noise;
z=z+System_Noise;
z_daoshu=z_daoshu+System_Noise;
R=zeros(1,k);
af=zeros(1,k);
beta=zeros(1,k);
for i=1:k
    R(i)=sqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i));
    af(i)=asin(z(i)/R(i));
    beta(i)=atan(y(i)/x(i));
end
VR=randn(1,k);
Vaf=randn(1,k);
Vbeta=randn(1,k);
R=R+VR*xitaR;
af=af+Vaf*xitaaf;
beta=beta+Vbeta*xitabeta;
Z_x=zeros(1,k);
Z_y=zeros(1,k);
Z_z=zeros(1,k);
for i=1:k
     Z_x(i)=R(i)*cos(beta(i))*cos(af(i));
     Z_y(i)=R(i)*cos(beta(i))*sin(af(i));
     Z_z(i)=R(i)*sin(beta(i));
end
% 计算直角坐标Rk
xitax=zeros(1,k);
xitay=zeros(1,k);
xitaz=zeros(1,k);
xitaxy=zeros(1,k);
xitaxz=zeros(1,k);
xitayz=zeros(1,k);
for i=1:k
        xitax(i)=cos(beta(i))*cos(beta(i))*cos(af(i))*cos(af(i))*xitaR+R(i)*R(i)*sin(beta(i))*sin(beta(i))*cos(af(i))*cos(af(i))*xitabeta+R(i)*R(i)*cos(beta(i))*cos(beta(i))*sin(af(i))*sin(af(i))*xitaaf;
        xitay(i)=cos(beta(i))*cos(beta(i))*sin(af(i))*sin(af(i))*xitaR+R(i)*R(i)*sin(beta(i))*sin(beta(i))*sin(af(i))*sin(af(i))*xitabeta+R(i)*R(i)*cos(beta(i))*cos(beta(i))*cos(af(i))*cos(af(i))*xitaaf;
        xitaz(i)=sin(beta(i))*sin(beta(i))*xitaR+R(i)*R(i)*cos(beta(i))*cos(beta(i))*xitabeta;
        xitaxy(i)=cos(beta(i))*cos(beta(i))*cos(af(i))*sin(af(i))*xitaR+R(i)*R(i)*sin(beta(i))*sin(beta(i))*cos(af(i))*sin(af(i))*xitabeta-R(i)*R(i)*cos(beta(i))*cos(beta(i))*cos(af(i))*sin(af(i))*xitabeta;
        xitaxz(i)=cos(beta(i))*cos(af(i))*sin(beta(i))*xitaR-R(i)*R(i)*sin(beta(i))*cos(af(i))*cos(beta(i))*xitabeta;
        xitayz(i)=cos(beta(i))*sin(af(i))*sin(beta(i))*xitaR-R(i)*R(i)*sin(beta(i))*sin(af(i))*cos(beta(i))*xitabeta;
end
% 卡尔曼滤波部分
s0=[1000 30 2000 40 3000 50]';
Z0=[Z_x(1),Z_y(1),Z_z(1)]';
s=[s0 zeros(Var_num,k-1)];
s_pre=[s0 zeros(Var_num,k-1)];
Z=[Z0 zeros(3,k-1) ];
for i=1:k
    Z(1,i)=Z_x(i);
    Z(2,i)=Z_y(i);
    Z(3,i)=Z_z(i);
end
fi=zeros(6,6);
for i=1:6
    fi(i,i)=1;
end
fi(1,2)=T;
fi(3,4)=T;
fi(5,6)=T;
H=zeros(3,6);
H(1,1)=1;
H(2,3)=1;
H(3,5)=1;
Qk=xita*eye(3,3);
G=zeros(6,3);
G(1,1)=T*T/2;
G(2,1)=T;
G(3,2)=T*T/2;
G(4,2)=T;
G(5,3)=T*T/2;
G(6,3)=T;
P0=G*Qk*G';
P1=P0;
P2=zeros(Var_num,Var_num);

s(:,1)=s0;
Rk=[xitax(1),xitaxy(1),xitaxz(1);xitaxy(1),xitay(1),xitayz(1);xitaxz(1),xitayz(1),xitaz(1)];
s_pre(:,1)=fi*s(:,1);
P2=fi*P1*fi'+G*Qk*G';
Kk = P2*H'*inv(H*P2*H'+Rk);    
P1 = (eye(Var_num,Var_num)-Kk*H)*P2; 
s(:,1) = s_pre(:,i)+Kk*(Z(:,1)-H*s_pre(:,1));
for i=2:k
    Rk=[xitax(i),xitaxy(i),xitaxz(i);xitaxy(i),xitay(i),xitayz(i);xitaxz(i),xitayz(i),xitaz(i)];
    s_pre(:,i)=fi*s(:,i-1);
    P2=fi*P1*fi'+G*Qk*G';
    Kk = P2*H'*inv(H*P2*H'+Rk);    
    P1 = (eye(Var_num,Var_num)-Kk*H)*P2; 
    s(:,i) = s_pre(:,i)+Kk*(Z(:,i)-H*s_pre(:,i));
end 
figure(1)
plot(s(1,:),'-b.');hold on;
plot(s(3,:),'-g.');hold on;
plot(s(5,:),'-r.');hold on;
figure(2)
plot(Z(1,:),'-b.');hold on;
plot(Z(2,:),'-g.');hold on;
plot(Z(3,:),'-r.');hold on;
figure(3)
plot(s_pre(1,:),'-b.');hold on;
plot(s_pre(3,:),'-g.');hold on;
plot(s_pre(5,:),'-r.');hold on;



    









    












