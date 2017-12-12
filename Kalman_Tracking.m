%kalman filterĿ�����
%����Ŀ�꺽��
%�趨��ֵ
%ʹ��EKF

N=100;%��������N��ʱ��
n=3; %״̬ά��
I=eye(6);
T=1;% ���Ʋ������ T
s=[3000,30,4000,40,5000,50]';     % ��ʼ״̬ % ����ֱ���˶��ٶ� 
sigma_u=5;           % �߶�������׼��
sigma_r=50 ;sigma_alpha=0.05 ;sigma_beta=0.04;       % �۲���������
vsigma=randn(N,1); %���ֵ����׼��Ϊ1�Ĺ۲�����

R=blkdiag(sigma_r^2,sigma_alpha^2,sigma_beta^2);  %%%�۲���������
Gu=[T^2/2;T];
G=blkdiag(Gu,Gu,Gu); %�߶�����G
Fu=[1,T;0,1];
F=blkdiag(Fu,Fu,Fu);  %״̬ת�ƾ���F

Q=eye(3)*sigma_u^2;     %%%״̬����
%H=[1 0;1 0;1 0];
P=100*eye(6);%��ʼ����P0
f=@(x)(F*x);  %״̬����
%�۲ⷽ��
%ֱ������ת��Ϊ������    
h=@(x)[sqrt(x(1)^2+x(3)^2+x(5)^2);asin(x(5)/sqrt(x(1)^2+x(3)^2+x(5)^2))*180/pi;atan(x(3)/x(1))*180/pi]; 
%������ת��Ϊֱ������
Trans=@(z)[z(1)*cos(z(3)/180*pi)*cos(z(2)/180*pi);z(1)*cos(z(3)/180*pi)*sin(z(2)/180*pi);z(1)*sin(z(3)/180*pi)];
x=s;


%������
xV = zeros(6,N);  %�˲����        
sV = zeros(n,N);  %ʵ�ʽ��
zV = zeros(n,N);%����ֵ
pV =zeros(6,N);%Ԥ��ֵ
XV = zeros(n,N);%ת��Ϊ��������˲����
SV = zeros(n,N);  %�������µ�ʵ�ʽ��
ZV = zeros(n,N);%ת��Ϊֱ������Ĳ���ֵ
PV =zeros(n,N);%������Ԥ��ֵ
%�˲�����
for k=1:N
  z = h(s) + R*randn(3,1);   %��������״̬����ֵ     
  pre=f(x)+G*Q*randn(3,1);%״̬Ԥ����
  Pre=h(pre);%�������µ�Ԥ��ֵ
  PV(:,k)=Pre;
  sV(:,k)=h(s);  %����������ʵֵ
  zV(:,k) =z; %״̬����ֵ
  Z=Trans(z);
  ZV(:,k)=Z;%ת��Ϊֱ������ϵ�Ĳ���ֵ12
  P_=F*P*F'+G*Q*G';         %���̷���Ԥ��
  [z1,H]=jaccsd(h,x); %����h���ſɱȾ���
  K=P_*H'*inv(H*P_*H'+R); %����������
  Z1=Trans(z1);
  x=x+K*(Z-Z1);        %״̬EKF����ֵ
  P=(I-K*H)*P_*(I-K*H)'+K*R*K';            %EKF����
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
legend('�۲�','Ԥ��','�˲�','��ʵ');
figure(2)
plot(zV(2,:),'-b.');hold on;
plot(PV(2,:),'-g.');hold on;
plot(XV(2,:),'-r.');hold on;
plot(sV(2,:),'-k.');hold on;
legend('�۲�','Ԥ��','�˲�','��ʵ');
figure(3)
plot(zV(3,:),'-b.');hold on;
plot(PV(3,:),'-g.');hold on
plot(XV(3,:),'-r.');hold on;
plot(sV(3,:),'-k.');hold on;
legend('�۲�','Ԥ��','�˲�','��ʵ');