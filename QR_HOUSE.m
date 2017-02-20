function [QT,Q,U,V,R] = QR_HOUSE(A);
%%%QR�ֽ⣬����n��������m
%%%A=Q*R;QT��Q��ת��;U,V�ֱ���QT��������������
n = size(A,1);
m = size(A,2);
t=min(n,m);
H = eye(n);
for k = 1:t
    a=A(k:n,k);
    if(length(a)>1)
    [v] = HOUSEHOLDER(a);
	Hk = eye(n);
	Hk(k:n,k:n) = eye(n+1-k) - (2/(v'*v))*(v*v');
	H = Hk*H;
    A = Hk*A;
    end
end
R=A;
QT=H;
Q=QT';
U=QT(1:t,:);
V=QT(t+1:n,:);
R=A(1:t,:);
