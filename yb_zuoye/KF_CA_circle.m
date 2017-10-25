clear all;
close all;
clc;
% 初始值设定
k=600;
Var_num=9;
t=1:1:k;
T=1;
xita=0.05;
xitaR=0.1;
xitaaf=0.1;
xitabeta=0.1;
x=zeros(1,k);
x_daoshu=zeros(1,k);
x_daoshu2=zeros(1,k);
y=zeros(1,k);
y_daoshu=zeros(1,k);
y_daoshu2=zeros(1,k);
z=zeros(1,k);
z_daoshu=zeros(1,k);
z_daoshu2=zeros(1,k);
Vx=zeros(1,k);
Vy=zeros(1,k);
Vz=zeros(1,k);
ax=zeros(1,k);
ay=zeros(1,k);
az=zeros(1,k);
for i=1:k
    V_l(i)=3;
    V_w(i)=0.005*pi;
end
x0=10;
y0=20;
z0=30;
x(1)=x0;
y(1)=y0;
z(1)=z0;
x_daoshu(1)=0;
y_daoshu(1)=0;
z_daoshu(1)=0;
x_daoshu(2)=0;
y_daoshu(2)=0;
z_daoshu(2)=0;
for n=2:k
    x(n)=x0+V_l(n)/V_w(n)*cos(V_w(n)*n);
    x_daoshu(n)=(x(n)-x(n-1))/T;
    x_daoshu2(n)=( x_daoshu(n)- x_daoshu(n-1))/T;
    y(n)=y0+V_l(n)/V_w(n)*sin(V_w(n)*n);
    y_daoshu(n)=(y(n)-y(n-1))/T;
    y_daoshu2(n)=( y_daoshu(n)- y_daoshu(n-1))/T;
    z(n)=z0;
    z_daoshu(n)=(z(n)-z(n-1))/T;
    z_daoshu2(n)=( z_daoshu(n)- z_daoshu(n-1))/T;
end
System_Noise=randn(1,k);
x=x+System_Noise;
x_daoshu=x_daoshu+System_Noise;
x_daoshu2=x_daoshu2+System_Noise;
y=y+System_Noise;
y_daoshu=y_daoshu+System_Noise;
y_daoshu2=y_daoshu2+System_Noise;
z=z+System_Noise;
z_daoshu=z_daoshu+System_Noise;
z_daoshu2=z_daoshu2+System_Noise;
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
s0=[10 3 0 20 4 0 30 5 0]';
Z0=[Z_x(1),Z_y(1),Z_z(1)]';
s=[s0 zeros(Var_num,k-1)];
s_pre=[s0 zeros(Var_num,k-1)];
Z=[Z0 zeros(3,k-1) ];
for i=1:k
    Z(1,i)=Z_x(i);
    Z(2,i)=Z_y(i);
    Z(3,i)=Z_z(i);
end
fi=zeros(9,9);
for i=1:9
    fi(i,i)=1;
end
fi(1,2)=T;
fi(1,3)=T*T/2;
fi(2,3)=T;
fi(4,5)=T;
fi(4,6)=T*T/2;
fi(5,6)=T;
fi(7,8)=T;
fi(7,9)=T*T/2;
fi(8,9)=T;
H=zeros(3,9);
H(1,1)=1;
H(2,4)=1;
H(3,7)=1;
Qk=xita*eye(3,3);
G=zeros(9,3);
G(1,1)=T*T/2;
G(2,1)=T;
G(3,1)=1;
G(4,2)=T*T/2;
G(5,2)=T;
G(6,2)=1;
G(7,3)=T*T/2;
G(8,3)=T;
G(9,3)=1;
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

plot(Z(1,:),'-b.');hold on;
plot(s_pre(1,:),'-g.');hold on;
plot(s(1,:),'-r.');hold on;



figure(2)
plot(Z(2,:),'-b.');hold on;
plot(s_pre(3,:),'-g.');hold on;
plot(s(3,:),'-r.');hold on;



figure(3)
plot(Z(3,:),'-b.');hold on;
plot(s_pre(5,:),'-g.');hold on;
plot(s(5,:),'-r.');hold on;

    

    









    












