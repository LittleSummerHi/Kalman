%kalman filter目标跟踪
%产生目标航迹
%设定初值
%使用EKF

N=100;%连续计算N各时刻
n=3; %状态维度
I=eye(6);
T=1;% 估计采样间隔 T
s=[3000,30,4000,40,5000,50]';     % 初始状态 % 匀速直线运动速度 
sigma_u=5;           % 策动噪声标准差
sigma_r=50 ;sigma_alpha=0.05 ;sigma_beta=0.04;       % 观测噪声方差
vsigma=randn(N,1); %零均值，标准差为1的观测噪声

R=blkdiag(sigma_r^2,sigma_alpha^2,sigma_beta^2);  %%%观测噪声矩阵
Gu=[T^2/2;T];
G=blkdiag(Gu,Gu,Gu); %策动矩阵G
Fu=[1,T;0,1];
F=blkdiag(Fu,Fu,Fu);  %状态转移矩阵F

Q=eye(3)*sigma_u^2;     %%%状态噪声
%H=[1 0;1 0;1 0];
P=100*eye(6);%初始方差P0
f=@(x)(F*x);  %状态方程
%观测方程
%直角坐标转换为极坐标    
h=@(x)[sqrt(x(1)^2+x(3)^2+x(5)^2);asin(x(5)/sqrt(x(1)^2+x(3)^2+x(5)^2))*180/pi;atan(x(3)/x(1))*180/pi]; 
%极坐标转化为直角坐标
Trans=@(z)[z(1)*cos(z(3)/180*pi)*cos(z(2)/180*pi);z(1)*cos(z(3)/180*pi)*sin(z(2)/180*pi);z(1)*sin(z(3)/180*pi)];
x=s;


%保存结果
xV = zeros(6,N);  %滤波结果        
sV = zeros(n,N);  %实际结果
zV = zeros(n,N);%测量值
pV =zeros(6,N);%预测值
XV = zeros(n,N);%转化为极坐标的滤波结果
SV = zeros(n,N);  %极坐标下的实际结果
ZV = zeros(n,N);%转化为直角坐标的测量值
PV =zeros(n,N);%极坐标预测值
%滤波过程
for k=1:N
  z = h(s) + R*randn(3,1);   %极坐标下状态测量值     
  pre=f(x)+G*Q*randn(3,1);%状态预测结果
  Pre=h(pre);%极坐标下的预测值
  PV(:,k)=Pre;
  sV(:,k)=h(s);  %极坐标下真实值
  zV(:,k) =z; %状态测量值
  Z=Trans(z);
  ZV(:,k)=Z;%转化为直角坐标系的测量值12
  P_=F*P*F'+G*Q*G';         %过程方差预测
  [z1,H]=jaccsd(h,x); %计算h的雅可比矩阵
  K=P_*H'*inv(H*P_*H'+R); %卡尔曼增益
  Z1=Trans(z1);
  x=x+K*(Z-Z1);        %状态EKF估计值
  P=(I-K*H)*P_*(I-K*H)'+K*R*K';            %EKF方差
  xV(:,k) = x;          %save
  X=h(x);
  XV(:,k)=X;
  s = f(s);  %update process 
  

end
figure(1)
plot(zV(1,:),'-b.');hold on;
plot(PV(1,:),'-g.');hold on;
plot(XV(1,:),'-r.');hold on;
plot(sV(1,:),'-k.');hold on;
legend('观测','预测','滤波','真实');
figure(2)
plot(zV(2,:),'-b.');hold on;
plot(PV(2,:),'-g.');hold on;
plot(XV(2,:),'-r.');hold on;
plot(sV(2,:),'-k.');hold on;
legend('观测','预测','滤波','真实');
figure(3)
plot(zV(3,:),'-b.');hold on;
plot(PV(3,:),'-g.');hold on
plot(XV(3,:),'-r.');hold on;
plot(sV(3,:),'-k.');hold on;
legend('观测','预测','滤波','真实');