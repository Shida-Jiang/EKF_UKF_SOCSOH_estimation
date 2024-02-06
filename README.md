<h3 style="text-align: center;">Codes for the paper "Weaknesses and Improvements of the Extended Kalman Filter for Battery State-of-Charge and State-of-Health Estimation", American Control Conference 2024</h3>
  
There are two data files and five scripts. A brief introduction of their functions is provided below:
* "saveP.mat" stores the OCV(SOC,SOH) function
* "Current_Profiles.mat" stores the UDDS current profile 
* "EKF_problem.mat" is a simple simulation where we compare the estimation results of KF, traditional EKF, and our improved EKF. You will see that the covariance estimation
 given by traditional EKF is very small, even though the system is unobservable. This suggests that the EKF can give inaccurate covariance estimation when the system is unobservable.
You can see more explanation to this example in our paper.
* "improved_EKF_simulation_CC.m" is the application of our EKF algorithm when using a constant current profile
* "improved_EKF_simulation_udds" is the application of our EKF algorithm when using a UDDS profile
* "improved_UKF_simulation_CC.m" is the application of the UKF algorithm when using a constant current profile
* "improved_UKF_simulation_udds.m" is the application of the UKF algorithm when using a UDDS profile
