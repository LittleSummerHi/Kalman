function [Pc,acur]=Kalman_filter(apre,Ppre,T,y,yavail)
%%%apre是前一个时刻的角度，Ppre是前一个时刻的P阵，T是时间间隔，y是观测量，yavail是是否存在当前时刻的观测量
%%%acur是当前时刻滤波得到的角度，Pc是当前时刻的P阵
wsigma=1e-4;vsigma=0.2;
I=eye(3);
F=[0 1 0;0 0 1;0 0 0];
Phi=I+F*T;%%%状态转移矩阵
G=[0 0 1]';%%%状态噪声系数矩阵
W=[0 0 wsigma];%%%状态噪声向量
Tao=([F*T]^2/6 +[F*T]/2+I)*G*T;
M=[1 0 0];

Q=wsigma^2;%%%状态噪声
R=vsigma^2;%%%观测噪声

P=Phi*Ppre*Phi'+Tao*Q*Tao';%%%最佳估计误差的方差

K=P*M'*inv(M*P*M'+R);%%%增益矩阵

if(yavail==1)
   acur=Phi*apre+K*(y-M*Phi*apre);
else
    acur=Phi*apre;
end
Pc=(I-K*M)*P;

