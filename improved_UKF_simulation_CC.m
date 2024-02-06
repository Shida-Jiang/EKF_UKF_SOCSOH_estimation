% Simulation study: ideal first-order RC model, 1C CC
% If you want to disable improvement 2,3, then change their value to 0.
improvement2=1; %reasonable Q and R
improvement3=1; %OCV-SOH
repeat=1000; %Number of trials, change the number to 1000 if you want to reproduce the results in Table 3 in the paper
%% Data Settings
T=25;%temeperature (10/25/40)
cycle=[0 100 200 350 500 600 700 800 900];
cells=[5 6 7 8 9 10];
OCVtype="_incremental_";
Crate='01';
%% Generate OCV curves
P = matfile('saveP.mat').P;
%% plot the curve
SOC=0:0.01:1;
SOH=0.8:0.001:1;
figure()
hold on
for i=1:11
    plot(SOC,get_OCV_from_fitting(SOC,P(i*20-19,:)),'DisplayName',append('SOH = ',num2str(SOH(i*20-19))))
end
xlabel('SOC')
ylabel('OCV (V)')
legend('location','southeast')
grid on

figure()
hold on
for i=1:11
    plot(SOC,get_dOCVdSOC_from_fitting(SOC,P(i*20-19,:)),'DisplayName',append('SOH = ',num2str(SOH(i*20-19))))
end
xlabel('SOC')
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
        U(i)=get_OCV_from_fitting(SOC_true(i+1),P(row,:))+U_c+I(i)*R1noisy+(sigma_V+sigma_OCV)*randn;
    end
    %Add noise to current
    I=I+randn(size(I))*sigma_I;
    Irms=rms(I);
    %% Our improved UKF algorithm
    % set some parameters
    % Values that are tracked
    K0=0;
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
    tic
    for j=1:length(I)
        %UKF
        %predict
        if(improvement2==1)
            a11=dt^2*sigma_I^2;
            a22=dt^2*sigma_R2C_reciporal^2*(state(2)-R2*I(j))^2+(R2-R2*exp(-dt*R2C_reciporal))^2*sigma_I^2+sigma_R2^2*I(j)*I(j)*(1-exp(-dt*R2C_reciporal))^2;
            a12=dt*(R2-R2*exp(-dt*R2C_reciporal))*sigma_I^2;
            processnoise=[a11 a12 0; a12 a22 0; 0 0 0];
            measurenoise=(sigma_V+sigma_OCV)^2+sigma_R1^2*I(j)^2+sigma_I^2*R1^2;
        else
            processnoise=[dt^2*sigma_I^2 0 0; 0 (Irms*R2*0.1)^2 0; 0 0 0];
            measurenoise=(sigma_V+sigma_OCV)^2;
        end
        L=real(sqrtm(Variance));
        states=zeros(3,3*2+1);
        states(:,1)=state;
        states(:,2)=state+sqrt(K0+3)*L(:,1);
        states(:,3)=state+sqrt(K0+3)*L(:,2);
        states(:,4)=state+sqrt(K0+3)*L(:,3);
        states(:,5)=state-sqrt(K0+3)*L(:,1);
        states(:,6)=state-sqrt(K0+3)*L(:,2);
        states(:,7)=state-sqrt(K0+3)*L(:,3);
        for i=1:7
            states(:,i)=stateModel(states(:,i),dt,I(j),R2,R2C_reciporal);
            if(i==1)
                state=states(:,i)*K0/(3+K0);
            else
                state=state+states(:,i)/(3+K0)/2;
            end
        end
        Variance=(states(:,1)-state)*(states(:,1)-state).'*K0/(3+K0)+processnoise;
        for i=2:7
            Variance=Variance+(state-states(:,i))*(state-states(:,i)).'/(3+K0)/2;
        end
        % Predict Measurement From Propagated Sigma Points
        measures=zeros(1,7);
        for i=1:7
            if(improvement3==1)
                row=floor((state(3)/Qmax/3600-0.8)*1000)+1;
                if(row<1)
                    row=1;
                elseif(row>200)
                    row=200;
                end
            else
                row=200;
            end
            measures(i)=get_OCV_from_fitting(states(1,i)/states(3,i),P(row,:))+states(2,i)+R1*I(j);
            if(i==1)
                m_exp=K0/(3+K0)*measures(i);
            else
                m_exp=m_exp+1/(3+K0)/2*measures(i);
            end
        end
        % Estimate Mean And Covariance of Predicted Measurement
        Py=K0/(3+K0)*(measures(1)-m_exp)*(measures(1)-m_exp).'+measurenoise;
        Pxy=K0/(3+K0)*(states(:,1)-state)*(measures(1)-m_exp);
        for i=2:7
            Py=Py+1/(3+K0)/2*(measures(i)-m_exp)*(measures(i)-m_exp).';
            Pxy=Pxy+1/(3+K0)/2*(states(:,i)-state)*(measures(i)-m_exp);
        end
        %kalman gain
        K=Pxy/Py;
        %update
        state=state+K*(U(j)-m_exp);
        Variance=Variance-K*Py*K.';
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
plot([0 t],SOC_est-SOC_true,'DisplayName','SOC error')
plot([0 t],Upperbound_SOC-SOC_est,'--r','DisplayName','Estimated error upper bound')
plot([0 t],Lowerbound_SOC-SOC_est,'--r','DisplayName','Estimated error lower bound')
legend()
grid on
xlabel('Time (s)')
ylabel('SOC error')
ylim([-0.002 0.002])

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
plot([0 t],SOH_est-SOH_true,'DisplayName','SOH error')
plot([0 t],Upperbound_SOH-SOH_est,'--r','DisplayName','Estimated error upper bound')
plot([0 t],Lowerbound_SOH-SOH_est,'--r','DisplayName','Estimated error lower bound')
legend()
grid on
xlabel('Time (s)')
ylabel('SOH error')
ylim([-0.02 0.02])

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
function P = OCV_function_fitting(SOH,cells,cycle,T,charge_SOH,OCVtype,Crate,degree)
charge_storage_coe=NaN(length(cells)*length(cycle),degree+1);
charge_storage_SOH=NaN(length(cells)*length(cycle),1);
discharge_storage_coe=NaN(length(cells)*length(cycle),degree+1);
discharge_storage_SOH=NaN(length(cells)*length(cycle),1);
num=0;
error=0;
if (OCVtype=="_incremental_")
    for ii=1:length(cells)
        cell=cells(ii);%5~10
        if(cell==5)
            channel='6';
        elseif cell==6
            channel='8';
        elseif cell==7
            channel='1';
        elseif cell==8
            channel='2';
        elseif cell==9
            channel='3';
        else
            channel='4';
        end
        name=append('SAMSUNG_Cell_',int2str(cell),OCVtype,Crate,'C_Channel_',channel,'_Wb_1.csv');
        for i=1:length(cycle)
            %read data
            path=append('DATA\',int2str(cycle(i)),' cycle\',int2str(T),' deg\',name);
            data=csvread(path,1,2);
            t=data(:,1);
            I=data(:,5);
            U=data(:,6);
            charge=0;
            %skip CCCV charge
            j=1;
            while I(j)>=0
                j=j+1;
            end
            %calculate charge capacity and discharge capacity
            %only record the data when I is about to be negative during
            %discharge and when I is about to be positive during charge
            Qcharge=0;
            Qdischarge=0;
            dOCV=[];
            cOCV=[];
            dSOC=[];
            cSOC=[];
            while j<=length(t)
                if(I(j)>0)
                    charge=1;
                end
                if(charge==0)
                    if(I(j)<0&&I(j-1)==0)
                        dSOC=[dSOC -Qdischarge];
                        %find the data during relaxation
                        milestone1=j-1;
                        milestone2=milestone1;
                        while I(milestone1-1)==0
                            milestone1=milestone1-1;
                        end
                        [xData, yData] = prepareCurveData( t(milestone1:milestone2)-t(milestone1), U(milestone1:milestone2));
                        %do a curve-fit to find the OCV
                        % Set up fittype and options.
                        ft = fittype( 'a*exp(b*x)+c', 'independent', 'x', 'dependent', 'y' );
                        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                        opts.Display = 'Off';
                        opts.StartPoint = [-0.01 -0.01 3.6];
                        % Fit model to data.
                        [fitresult, gof] = fit( xData, yData, ft, opts );
                        coe=coeffvalues(fitresult);
                        if(coe(3)>2.5)
                            dOCV=[dOCV coe(3)];
                        else
                            dOCV=[dOCV yData(length(yData))];
                        end
                        
                    end
                    Qdischarge=Qdischarge-I(j)*(t(j)-t(j-1))/3600;
                else
                    if(I(j)==0&&I(j-1)>0&&length(I)-j>60)
                        cSOC=[cSOC Qcharge];
                        %find the data during relaxation
                        milestone1=j;
                        milestone2=milestone1;
                        while milestone2+1<length(I)&&I(milestone2+1)==0
                            milestone2=milestone2+1;
                        end
                        [xData, yData] = prepareCurveData( t(milestone1:milestone2)-t(milestone1), U(milestone1:milestone2));
                        %do a curve-fit to find the OCV
                        % Set up fittype and options.
                        ft = fittype( 'a*exp(b*x)+c', 'independent', 'x', 'dependent', 'y' );
                        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                        opts.Display = 'Off';
                        if(Qcharge>0)
                            % Fit model to data.
                            opts.StartPoint = [0.01 -0.01 3.6];
                            [fitresult, ~] = fit( xData, yData, ft, opts );
                            coe=coeffvalues(fitresult);
                            cOCV=[cOCV coe(3)];
                        else
                            cOCV=[cOCV yData(length(yData))];
                        end     
                    end
                    Qcharge=Qcharge+I(j)*(t(j)-t(j-1))/3600;
                end
                j=j+1;
            end
            dSOC=(dSOC+Qdischarge)/Qdischarge;
            cSOC=cSOC/Qcharge;
            dSOH=Qdischarge/2.2;
            cSOH=Qcharge/2.2;
            if(i==1)
                cref=cSOH;
                dref=dSOH;
            end
            cSOH=cSOH/cref;
            dSOH=dSOH/dref;
            %store the fitting
            [xData, yData] = prepareCurveData( dSOC, dOCV );
            ft = fittype( append('poly',int2str(degree)) );
            [fitresult, gof] = fit( xData, yData, ft );
            error=error+gof.rmse;
            if(gof.rmse>0.01)
                1;
            end
            num=num+1;
            discharge_storage_coe((ii-1)*length(cycle)+i,:)=coeffvalues(fitresult);
            discharge_storage_SOH((ii-1)*length(cycle)+i)=dSOH;

            [xData, yData] = prepareCurveData( cSOC, cOCV );
            ft = fittype( append('poly',int2str(degree)) );
            [fitresult, gof] = fit( xData, yData, ft );
            error=error+gof.rmse;
            if(gof.rmse>0.01)
                1;
            end
            num=num+1;
            charge_storage_coe((ii-1)*length(cycle)+i,:)=coeffvalues(fitresult);
            charge_storage_SOH((ii-1)*length(cycle)+i)=cSOH;
        end
    end
end
%define the outputs

%OCV(SOC) is a ninth order polynomial function, the parameter matrix is
%below
P=NaN(length(SOH),degree+1);
for i=1:length(SOH)
    if(charge_SOH==1)
        distance=(SOH(i)-charge_storage_SOH).^2;
        weight=1./(distance+1e-4)./sum(1./(distance+1e-4));
        P(i,:)=weight.'*charge_storage_coe;
    else
        distance=(SOH(i)-discharge_storage_SOH).^2;
        weight=1./(distance+1e-4)./sum(1./(distance+1e-4));
        P(i,:)=weight.'*discharge_storage_coe;
    end
end
%error of the OCV curve
%error=error/num
end