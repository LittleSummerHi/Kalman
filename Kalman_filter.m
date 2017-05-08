function [Pc,acur]=Kalman_filter(apre,Ppre,T,y,yavail)
%%%apre��ǰһ��ʱ�̵ĽǶȣ�Ppre��ǰһ��ʱ�̵�P��T��ʱ������y�ǹ۲�����yavail���Ƿ���ڵ�ǰʱ�̵Ĺ۲���
%%%acur�ǵ�ǰʱ���˲��õ��ĽǶȣ�Pc�ǵ�ǰʱ�̵�P��
wsigma=1e-4;vsigma=0.2;
I=eye(3);
F=[0 1 0;0 0 1;0 0 0];
Phi=I+F*T;%%%״̬ת�ƾ���
G=[0 0 1]';%%%״̬����ϵ������
W=[0 0 wsigma];%%%״̬��������
Tao=([F*T]^2/6 +[F*T]/2+I)*G*T;
M=[1 0 0];

Q=wsigma^2;%%%״̬����
R=vsigma^2;%%%�۲�����

P=Phi*Ppre*Phi'+Tao*Q*Tao';%%%��ѹ������ķ���

K=P*M'*inv(M*P*M'+R);%%%�������

if(yavail==1)
   acur=Phi*apre+K*(y-M*Phi*apre);
else
    acur=Phi*apre;
end
Pc=(I-K*M)*P;

