Z=(1:100); %�۲�ֵ    
noise=randn(1,100); %����Ϊ1�ĸ�˹����    
Z=Z+noise;    
%plot(noise) 
X=[0; 0]; %״̬    
Sigma = [1 0; 0 1]; %״̬Э�������    
F=[1 1; 0 1]; %״̬ת�ƾ���    
Q=[0.0001, 0; 0 0.0001]; %״̬ת��Э�������    
H=[1 0]; %�۲����    
R=1; %�۲���������    
figure;    
hold on;    
    
for i=1:100    
    
  X_ = F*X;    
  Sigma_ = F*Sigma*F'+Q;    
 K = Sigma_*H'/(H*Sigma_*H'+R);    
  X = X_+K*(Z(i)-H*X_);    
  Sigma = (eye(2)-K*H)*Sigma_;    
plot(X(1), X(2), '.','MarkerSize',10); %���㣬�����ʾλ�ã������ʾ�ٶ�    
end  

plot([0,100],[1,1],'r-');   
