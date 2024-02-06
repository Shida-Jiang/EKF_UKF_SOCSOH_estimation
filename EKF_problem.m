%A simple experiment that demonstrates the problem of EKF
iterations=2;%feel free to change this number
x1=[1;0.99];
x2=[1;0.99];
x3=[1;0.99];
%true state: x1=1 x2=1
P1=[1 0;0 1];
P2=[1 0;0 1];
P3=[1 0;0 1];
Q=1e-8;
for j=1:iterations
    r=0;% Actual measurement noise is zero
    %KF
    H1=[2 1];
    S1=H1*P1*H1.'+Q;
    K1=P1*H1.'/S1;
    x1=x1+K1*(3+r-2*x1(1)-x1(2));
    P1=([1 0; 0 1]-K1*H1)*P1;
    %Traditional EKF
    H2=[2*x2(1) 1];
    S2=H2*P2*H2.'+Q; 
    K2=P2*H2.'/S2; 
    x2=x2+K2*(2+r-x2(1)*x2(1)-x2(2));
    P2=([1 0; 0 1]-K2*H2)*P2;
    %Our EKF
    H3=[2*x3(1) 1];
    S3=H3*P3*H3.'+Q; 
    K3=P3*H3.'/S3; 
    x3=x3+K3*(2+r-x3(1)*x3(1)-x3(2));
    H3_1=[2*x3(1) 1];%The extra step
    if(trace(([1 0; 0 1]-K3*H3_1)*P3)>trace(P3) || trace(([1 0; 0 1]-K3*H3_1)*P3)<0)
        P3=([1 0; 0 1]-K3*H3)*P3;
    else
        P3=([1 0; 0 1]-K3*H3_1)*P3;
    end  
end
disp(['After ',int2str(iterations),' iterations,'])
disp('The covariance matrix of KF is:')
P1
disp('The covariance matrix of traditional EKF is:')
P2
disp('The covariance matrix of our improved EKF is:')
P3
