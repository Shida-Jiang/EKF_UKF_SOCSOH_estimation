% Simulation study: ideal first-order RC model, 1C CC
% If you want to disable improvement 1,2,3, then change their value to 0.
improvemnet1=1; %new EKF
improvement2=1; %reasonable Q and R
improvement3=1; %OCV-SOH
repeat=1000; %Number of trials, change the number to 1000 if you want to reproduce the results in Table 2 in the paper; change the number to 3 if you want to reproduce Figure 4,5 in the paper
%% Generate OCV curves
P = matfile('saveP.mat').P;
%% plot the curve
SOC=0:0.01:1;
SOH=0.8:0.001:1;
figure()
hold on
for i=1:11
    plot(SOC*100,get_OCV_from_fitting(SOC,P(i*20-19,:)),'DisplayName',append('SOH = ',num2str(100*SOH(i*20-19)),"%"))
end
xlabel('SOC (%)')
ylabel('OCV (V)')
legend('location','southeast')
grid on

figure()
hold on
for i=1:11
    plot(SOC,get_dOCVdSOC_from_fitting(SOC,P(i*20-19,:)),'DisplayName',append('SOH = ',num2str(100*SOH(i*20-19)),"%"))
end
xlabel('SOC (%)')
ylabel('\partial OCV/\partial SOC (V)')
legend('location','southeast')
ylim([0 2])
grid on
%% EKF "true" data generation: CC
rng('default');
rng(100)
%setting some parameters
Qmax=2.2;
R1=0.08;
R2=0.05;
R2C_reciporal=0.008;
sigma_R1=0.008;
sigma_R2=0.005;
sigma_R2C_reciporal=0.0008;
sigma_I=5e-4;
sigma_OCV=4e-3;
sigma_V=1e-3;
%Values that are tracked
dt=1;
t=dt:dt:(3600+3600+300+300);
I=zeros(1,length(t));
U=zeros(1,length(t));
SOH_true=ones(1,length(t)+1)*0.9;
SOC_true=zeros(1,length(t)+1);
for i=1:3600
    I(i+300)=-Qmax*SOH_true(1);
    I(length(I)+1-i)=Qmax*SOH_true(1);
end
SOH_RMS=0;
SOC_RMS=0;
SOH_RMS_est=0;
SOC_RMS_est=0;
SOH_inc=0;
SOC_inc=0;
AvgTime=0;
for iii=1:repeat
    % initialize state
    Uc=0;
    SOC_true(1)=1;
    U_c=0;
    %start simulation
    for i=1:length(I)
        R1noisy=R1+sigma_R1*randn;
        R2noisy=R2+sigma_R2*randn;
        R2C_reciporalnoisy=R2C_reciporal+sigma_R2C_reciporal*randn;

        SOC_true(i+1)=SOC_true(i)+I(i)*dt/3600/(SOH_true(i)*Qmax);
        U_c=U_c*exp(-dt*R2C_reciporalnoisy)+(R2noisy-R2noisy*exp(-dt*R2C_reciporalnoisy))*I(i);
        row=101;
        U(i)=get_OCV_from_fitting(SOC_true(i+1),P(row,:))+U_c+I(i)*R1noisy+(sigma_OCV+sigma_V)*randn;
    end
    %Add noise to current
    I=I+randn(size(I))*sigma_I;
    Irms=rms(I);
    %% Our improved EKF algorithm
    % set some parameters
    % Values that are tracked
    intimeSOC=1;
    intimeSOH=1;
    SOH_est=zeros(1,length(t)+1);
    Upperbound_SOH=zeros(1,length(t)+1);
    Lowerbound_SOH=zeros(1,length(t)+1);
    SOC_est=zeros(1,length(t)+1);
    Upperbound_SOC=zeros(1,length(t)+1);
    Lowerbound_SOC=zeros(1,length(t)+1);
    %Initialization
    state0=[0.9*0.95*Qmax*3600 0 Qmax*3600*0.95];%initial state: remaining capacity, Uc, maximum capacity
    state=state0.';
    initialCovariance=[(0.05*Qmax*3600)^2 0 0;0 0 0;0 0 (0.05*Qmax*3600)^2];
    Variance = initialCovariance;
    SOC_est(1)=state(1)/state(3);
    SOH_est(1)=state(3)/Qmax/3600;
    Upperbound_SOH(1)=(state(3)+2*sqrt(Variance(3,3)))/Qmax/3600;
    Lowerbound_SOH(1)=(state(3)-2*sqrt(Variance(3,3)))/Qmax/3600;
    Upperbound_SOC(1)=state(1)/state(3)+2*sqrt(Variance(1,1)/state(3)^2+Variance(3,3)*state(1)^2/state(3)^4-2*state(1)/state(3)^3*Variance(1,3));
    Lowerbound_SOC(1)=state(1)/state(3)-2*sqrt(Variance(1,1)/state(3)^2+Variance(3,3)*state(1)^2/state(3)^4-2*state(1)/state(3)^3*Variance(1,3));
    %estimation, record the error & upper and lower bound
    %SOH & SOC estimation during discharge
    tic;
    for j=1:length(I)
        %EKF
        %predict
        if(improvement2==1)
            a11=dt^2*sigma_I^2;
            a22=dt^2*sigma_R2C_reciporal^2*(state(2)-R2*I(j))^2+(R2-R2*exp(-dt*R2C_reciporal))^2*sigma_I^2+sigma_R2^2*I(j)*I(j)*(1-exp(-dt*R2C_reciporal))^2;
            a12=dt*(R2-R2*exp(-dt*R2C_reciporal))*sigma_I^2;
            processnoise=[a11 a12 0; a12 a22 0; 0 0 0];
            measurenoise=(sigma_OCV+sigma_V)^2+sigma_R1^2*I(j)^2+sigma_I^2*R1^2;
        else
            processnoise=[dt^2*sigma_I^2 0 0; 0 (Irms*R2*0.1)^2 0; 0 0 0];
            measurenoise=(sigma_OCV+sigma_V)^2;
        end
        state=stateModel(state,dt,I(j),R2,R2C_reciporal);
        F=get_F(dt,R2C_reciporal);
        Variance=F*Variance*F.'+processnoise;
        %H matrix
        if(improvement3==1)
            row=floor((state(3)/Qmax/3600-0.8)*1000)+1;
            if(row<1)
                row=1;
            elseif(row>200)
                row=200;
            end
            dOCVdSOH = get_dOCVdSOH_from_fitting(state(1)/state(3),P(row+1,:),P(row,:),1e-3);
        else
            row=200;
            dOCVdSOH=0;
        end
        dOCVdSOC = get_dOCVdSOC_from_fitting(state(1)/state(3),P(row,:));
        H=[dOCVdSOC/state(3) 0 dOCVdSOH/(Qmax*3600)-dOCVdSOC*state(1)/state(3)/state(3)];
        %kalman gain
        K=Variance*H.'/(H*Variance*H.'+measurenoise);
        %update
        m_exp=get_OCV_from_fitting(state(1)/state(3),P(row,:))+state(2)+R1*I(j);
        state=state+K*(U(j)-m_exp);
        %H matrix update
        if(improvemnet1==1)
            row=floor((state(3)/Qmax/3600-0.8)*1000)+1;
            if(row<1)
                row=1;
            elseif(row>200)
                row=200;
            end
            dOCVdSOH = get_dOCVdSOH_from_fitting(state(1)/state(3),P(row+1,:),P(row,:),1e-3);
            dOCVdSOC = get_dOCVdSOC_from_fitting(state(1)/state(3),P(row,:));
            H2=[dOCVdSOC/state(3) 0 dOCVdSOH/(Qmax*3600)-dOCVdSOC*state(1)/state(3)/state(3)];
            temp = ([1,0,0;0,1,0;0,0,1]-K*H2)*Variance;
            if(trace(temp)<trace(Variance) && trace(temp)>=0)
                Variance=temp;
            else
                Variance=([1,0,0;0,1,0;0,0,1]-K*H)*Variance;
            end
        else
            Variance=([1,0,0;0,1,0;0,0,1]-K*H)*Variance;
        end

        
        %record the values
        SOC_est(j+1)=state(1)/state(3);
        SOH_est(j+1)=state(3)/Qmax/3600;
        Upperbound_SOH(j+1)=(state(3)+2*sqrt(Variance(3,3)))/Qmax/3600;
        Lowerbound_SOH(j+1)=(state(3)-2*sqrt(Variance(3,3)))/Qmax/3600;
        Upperbound_SOC(j+1)=state(1)/state(3)+2*sqrt(abs(Variance(1,1)/state(3)^2+Variance(3,3)*state(1)^2/state(3)^4-2*state(1)/state(3)^3*Variance(1,3)));
        Lowerbound_SOC(j+1)=real(state(1)/state(3)-2*sqrt(abs(Variance(1,1)/state(3)^2+Variance(3,3)*state(1)^2/state(3)^4-2*state(1)/state(3)^3*Variance(1,3))));
        if(Upperbound_SOC(j+1)>SOC_true(j+1) && Lowerbound_SOC(j+1)<SOC_true(j+1))
            intimeSOC=intimeSOC+1;
        end
        if(Upperbound_SOH(j+1)>SOH_true(j+1) && Lowerbound_SOH(j+1)<SOH_true(j+1))
            intimeSOH=intimeSOH+1;
        end
    end
    AvgTime=AvgTime+toc/repeat;
    SOH_RMS=SOH_RMS+rms(SOH_est-SOH_true)/repeat;
    SOC_RMS=SOC_RMS+rms(SOC_est-SOC_true)/repeat;
    SOH_RMS_est=SOH_RMS_est+rms((Upperbound_SOH-SOH_est)/2)/repeat;
    SOC_RMS_est=SOC_RMS_est+rms((Upperbound_SOC-SOC_est)/2)/repeat;
    SOH_inc=SOH_inc+intimeSOH/(length(t)+1)/repeat;
    SOC_inc=SOC_inc+intimeSOC/(length(t)+1)/repeat;
end
disp(AvgTime*1000)
%% Plotting the results
%SOC estimation result
figure(2)
hold on
plot([0 t],SOC_true,'DisplayName','True SOC')
plot([0 t],SOC_est,'DisplayName','Estimated SOC')
plot([0 t],Upperbound_SOC,'--r','DisplayName','Estimated upper bound')
plot([0 t],Lowerbound_SOC,'--r','DisplayName','Estimated lower bound')
legend()
grid on
xlabel('Time (s)')
ylabel('SOC')

figure(4)
hold on
plot([0 t],100*(SOC_est-SOC_true),'DisplayName','SOC error')
plot([0 t],100*(Upperbound_SOC-SOC_est),'--r','DisplayName','Estimated error upper bound')
plot([0 t],100*(Lowerbound_SOC-SOC_est),'--r','DisplayName','Estimated error lower bound')
legend()
grid on
xlabel('Time (s)')
ylabel('SOC error (%)')
ylim([-0.5 0.5])

%SOH estimation result
figure(3)
hold on
plot([0 t],SOH_true,'DisplayName','True SOH')
plot([0 t],SOH_est,'DisplayName','Estimated SOH')
plot([0 t],Upperbound_SOH,'--r','DisplayName','Estimated upper bound')
plot([0 t],Lowerbound_SOH,'--r','DisplayName','Estimated lower bound')
legend()
grid on
xlabel('Time (s)')
ylabel('SOH')

figure(5)
hold on
plot([0 t],100*(SOH_est-SOH_true),'DisplayName','SOH error')
plot([0 t],100*(Upperbound_SOH-SOH_est),'--r','DisplayName','Estimated error upper bound')
plot([0 t],100*(Lowerbound_SOH-SOH_est),'--r','DisplayName','Estimated error lower bound')
legend()
grid on
xlabel('Time (s)')
ylabel('SOH error (%)')
ylim([-2 2])
%% display metrics
disp(append("SOC RMS error =",num2str(100*SOC_RMS),"%"))
disp(append("SOH RMS error =",num2str(100*SOH_RMS),"%"))
disp(append("Estimated SOC RMS error =",num2str(100*SOC_RMS_est),"%"))
disp(append("Estimated SOH RMS error =",num2str(100*SOH_RMS_est),"%"))
disp(append("SOC Inclusion rate = ",num2str(SOC_inc)))
disp(append("SOH Inclusion rate = ",num2str(SOH_inc)))
%% functions
function stateNext = stateModel(state,dt,I,R2,R2C_reciporal)
    A = [1 0 0; 0 exp(-dt*R2C_reciporal) 0; 0 0 1];
    B = [dt; R2-R2*exp(-dt*R2C_reciporal); 0];
    stateNext = A*state+B*I;
end
function F_matrix=get_F(dt,R2C_reciporal)
    F_matrix=[1 0 0; 0 exp(-dt*R2C_reciporal) 0; 0 0 1];
end
function dOCVdSOH = get_dOCVdSOH_from_fitting(SOC,coeficients1,coeficients2,deltaSOH)
xspace=NaN(length(coeficients1),length(SOC));
for i=1:length(coeficients1)
    xspace(i,:)=SOC.^(length(coeficients1)-i);
end
dOCVdSOH=(coeficients1-coeficients2)*xspace/deltaSOH;
end

function dOCVdSOC = get_dOCVdSOC_from_fitting(SOC,coeficients)
xspace=zeros(length(coeficients),length(SOC));
for i=1:length(coeficients)-1
    xspace(i,:)=(length(coeficients)-i).*SOC.^(length(coeficients)-i-1);
end
dOCVdSOC=coeficients*xspace;
end
function OCV = get_OCV_from_fitting(SOC,coeficients)
xspace=NaN(length(coeficients),length(SOC));
for i=1:length(coeficients)
    xspace(i,:)=SOC.^(length(coeficients)-i);
end
OCV=coeficients*xspace;
end
