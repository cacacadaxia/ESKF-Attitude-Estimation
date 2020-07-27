clear;
close all;
% clc;
addpath('utils');
tic
%% ==========================Input Data======================== %%
dataDir = '';
fileName = 'NAV_1'; 
Data = importdata(sprintf('%s.mat',fileName));
length = size(Data,1);
dataStructureDef;
%% ==========================Initial State======================== %%
%get init euler angle and quaternion
[qInit,roll(1),pitch(1),yaw(1)] = initializequat(Data);
%Set up the noise parameters
wn_var  = 1e-5 * ones(1,3);               % rot vel var
wbn_var = 1e-9* ones(1,3);                % gyro bias change var
an_var  = 1e-3 * ones(1,3);               % acc var
mn_var  = 1e-4 * ones(1,3);               % mag var
noiseParams.Q = diag([wn_var, wbn_var]); 
noiseParams.R = diag([an_var, mn_var]); 
% noiseParams.R = diag([an_var]);
%set uo P initial value
q_var_init = 1e-5 * ones(1,3);            % init rot var
wb_var_init = 1e-7 * ones(1,3);           % init gyro bias var
imuErrorStates{1}.P = diag([q_var_init, wb_var_init]);
%Use ground truth for the first state
imuNominalStates{1}.q       = qInit;
imuNominalStates{1}.wb      = wbn_var';
%set up the reset value of error state
detThetaReset = zeros(3,1);             
imuErrorStates{1}.det_theta = detThetaReset;
imuErrorStates{1}.det_wb    = detThetaReset; 
%MAIN LOOP
omega_b = zeros(3,1);
dt = 0.02;
quat = qInit;
PP = imuErrorStates{1}.P*dt^2;
R = noiseParams.R;
for state_k = 1:length-1
    index = state_k;
    %% ==========================STATE PROPAGATION======================== %%
    %Propagate nominal state 
    imuNominalStates{state_k+1}  = propagateNominalState(imuNominalStates{state_k},measurements{state_k},measurements{state_k+1});
    %Propagate error state and covariance
    imuErrorStates{state_k+1}    = propagateErrorStateAndCovar(imuNominalStates{state_k}.wb,imuErrorStates{state_k}, measurements{state_k}, noiseParams); 
    %% ==========================FILTER UPDATE======================== %%
    %Calculate H and detZ = y-h(det_x)
    [H,detZ] = calH(imuNominalStates{state_k+1}.q,  measurements{state_k+1});
    P = imuErrorStates{state_k+1}.P;
    P_save(:,:,state_k) = P;
    % Calculate Kalman gain
    K = (P*H') / ( H*P*H' + noiseParams.R); 
    % State correction
    det_x =  K * detZ;
    tp(:,state_k) = det_x;
    imuErrorStates{state_k+1}.det_theta = det_x(1:3);
    imuErrorStates{state_k+1}.det_wb    = det_x(4:6);
    % Covariance correction
    imuErrorStates{state_k+1}.P = P -  K *(H*P*H'+ noiseParams.R)*K';
  %% ==========================STATE CORRECTION======================== %%
    det_q = buildUpdateQuat(imuErrorStates{state_k+1}.det_theta);                                       %it seems more safety to use this to update the euation q=[1 1/2*det_theta]
    imuNominalStates{state_k+1}.q  = quatLeftComp(imuNominalStates{state_k+1}.q)*det_q;                 %joan sola_Quaternion kinematics for the error_state KF p24
    imuNominalStates{state_k+1}.q  = imuNominalStates{state_k+1}.q /norm(imuNominalStates{state_k+1}.q ); %%归一化
    imuNominalStates{state_k+1}.wb = imuNominalStates{state_k+1}.wb + imuErrorStates{state_k+1}.det_wb;
    [roll(state_k+1),pitch(state_k+1),yaw(state_k+1)] = quattoeuler(imuNominalStates{state_k+1}.q);     %transform the quaternion to euler for plot
    %% ==========================ERROR STATE RESET======================== %%
    imuErrorStates{state_k+1}.det_theta = detThetaReset;
    imuErrorStates{state_k+1}.det_wb    = detThetaReset;
    imuErrorStates{state_k+1}.P         = eye(6)*imuErrorStates{state_k+1}.P*eye(6)';
%     err_sigma(:,state_k) = sqrt(diag(imuErrorStates{state_k+1}.P));



%% 另外一种
omega_m = (measurements{index+1}.omega + measurements{index}.omega)/2;
av = (omega_m - omega_b)*dt;
 det_q = [1;0.5*av];
      quat = quatLeftComp(quat)*det_q;
     omega_b = omega_b;
         F1 = axisAngleToRotMat((measurements{index}.omega - omega_b)*dt)';
    Fx = [F1  , -eye(3)*dt;
    zeros(3)  , eye(3)];

sigma_wn = 1e-5;
sigma_wbn = 1e-9;
Theta_i = sigma_wn*dt^2*eye(3);
Omega_i = sigma_wbn*dt*eye(3);
Fi = eye(6);
Qi = blkdiag(Theta_i , Omega_i);
Q = Fi*Qi*Fi';

     P_ = Fx*PP*Fx' + Q;
     P_ = enforcePSD(P_);
          [H,detZ] = calH(quat,  measurements{index+1});
     K = P_*H'*inv(H*P_*H' + R)/2;
     err_state = K*detZ;
     
     PP = P_ - K*(H*P_*H' + R)*K';
          dt_theta = err_state(1:3);
     dt_omega_b = err_state(4:6);
     dt_q = buildUpdateQuat(dt_theta);
     quat = quatLeftComp(quat)*dt_q;
     quat = quat/norm(quat);
     omega_b = omega_b + dt_omega_b;
     
     % ------save angle-----------------------------
     
     [a1,a2,a3] = quattoeuler(quat);
     oula(index+1,:) = [a1,a2,a3]/180*pi;
     dt_theta_save(index+1,:) = err_state';
     % ----------------------------reset-------------------
     err_state = zeros(6,1);
     G = blkdiag(eye(3) - omegaMatrix(dt_theta/2) ,eye(3));
     PP = G*PP*G';
end 
% err_sigma(:,length)=err_sigma(:,end);
toc
% 
figure;
subplot(3,1,1)
plot(oula(:,1)/pi*180 - rollRef);
hold on;
plot(roll - rollRef);
plot(oula(:,1)/pi*180 - rollRef);
subplot(3,1,2)
plot(oula(:,2)/pi*180 - pitchRef);
hold on;
plot(pitch - pitchRef);
subplot(3,1,3)
plot(oula(:,3)/pi*180 - yawRef);
hold on;
plot(yaw - yawRef);
legend 博主 我的



figure;
subplot(3,1,1)
plot(oula(:,1)/pi*180 - rollRef);
subplot(3,1,2)
plot(oula(:,2)/pi*180 - pitchRef);
subplot(3,1,3)
plot(oula(:,3)/pi*180 - yawRef);
legend 1 2 
% figure;
% subplot(3,1,1)
% plot(roll - rollRef);
% subplot(3,1,2)
% plot(pitch - pitchRef);
% subplot(3,1,3)
% plot(yaw - yawRef);
% legend 1 2 

function [omega]=omegaMatrix(data)
% wx=data(1)*pi/180;
% wy=data(2)*pi/180;
% wz=data(3)*pi/180;
wx=data(1);
wy=data(2);
wz=data(3);
omega=[
    0,-wz,wy;
    wz,0,-wx;
    -wy,wx,0
    ];
end

function Ang3 = q2oula(q)
%四元数转欧拉角
x = atan2(2*(q(1)*q(2)+q(3)*q(4)),1 - 2*(q(2)^2+q(3)^2));
y = asin(2*(q(1)*q(3) - q(2)*q(4)));
z = atan2(2*(q(1)*q(4)+q(2)*q(3)),1 - 2*(q(3)^2+q(4)^2));
Ang3 = [x y z];
end