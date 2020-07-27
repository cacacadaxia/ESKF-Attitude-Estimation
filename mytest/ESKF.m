% =========================================================================
%
%                  ESKF的实现
% 
%
% =========================================================================
%
%　(C)2020-2022 China Academy of Railway Sciences 
%   版本：V1.0
%   日期：2020年 6月23日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能：1.实现ESKF算法，加深对于状态估计的理解
%        2.其中的问题：
%     1） 测量加上地磁计，可以降低误差
%     2） 注意误差量与标称量，在预测时，主要更新的是"标称状态量[q,wb]"，和"error state 的 P矩阵"。
%     3） 四元数转角度时有误差，有一个角度固定漂移，这个是怎么回事？
%     4） 如何确定陀螺的参数与协方差矩阵之间的关系？
%
%--------------------------------------------------------------------------
clear all;
close all;
addpath('../../ESKF-Attitude-Estimation-master')
addpath('../utils')
% -------------------- import data-------------------
fileName = '../NAV_1';
Data = importdata(sprintf('%s.mat',fileName));
lengthtp = size(Data,1);

time    = zeros(lengthtp,1);
roll    = zeros(lengthtp,1);
pitch   = zeros(lengthtp,1);
yaw     = zeros(lengthtp,1);
imuNominalStates = cell(1,lengthtp);
imuErrorStates   = cell(1,lengthtp);
measurements = cell(1,lengthtp);
%groundTruth
for state_k = 1:lengthtp 
    measurements{state_k}.dt    = 0.02;                      % sampling times 50Hz
    measurements{state_k}.omega = Data(state_k,27:29)';            
    measurements{state_k}.acc   = Data(state_k,9:11)';
    measurements{state_k}.mag   = Data(state_k,15:17)';
    time(state_k)=state_k*0.02;
end
rad2deg = 180/pi;
rollRef   = Data(:,30)*rad2deg;
pitchRef  = Data(:,31)*rad2deg;
yawRef    = Data(:,32)*rad2deg;
% --------------------Data analysis------------------
% ++++++++++++++++++++1.initialization++++++++++++++++
dt = measurements{1}.dt;

% Nominal state
omega_b = zeros(3,1);%%这个用到
theta = zeros(3,1);%%这个用不到
quat = zeros(4,1);
% error state initialization
dt_theta = zeros(3,1);
dt_omega_b = zeros(3,1);
err_state = [dt_theta;dt_omega_b];

% --------Refer to previous practice for initialization-----------------------------------
init_angle = [Data(1,30),Data(1,31),Data(1,32)]';
init_quat = oula2q(init_angle);
quat = init_quat';
% -------------------------2.covariance matrix ---------------------
% Uncertainty matrix
p1 = 1e-4;p2 =  1e-5;
P = blkdiag(p1,p1,p1,p2,p2,p2);%%初始化
% Process noise matrix(including control quantity and noise quantity)
sigma_wn_2 = 1e-6;
sigma_wbn_2 = 1e-9;
Theta_i = sigma_wn_2*dt^2*eye(3);
Omega_i = sigma_wbn_2*dt*eye(3);
Fi = eye(6);
Qi = blkdiag(Theta_i , Omega_i);
Q = Fi*Qi*Fi';

sigma_acc_2 = 1e-3;
sigma_mn_2 = 1e-4;
R = blkdiag(eye(3)*sigma_acc_2,eye(3)*sigma_mn_2);
for index = 1:lengthtp-1
    % --------------------------forecast------------
    omega_m = (measurements{index+1}.omega + measurements{index}.omega)/2;
    av = (omega_m - omega_b)*dt;
     det_q = [1;0.5*av];
     % Update nominal value
     quat = quatLeftComp(quat)*det_q;
     omega_b = omega_b;
     % Update P matrix
    F1 = Exp_lee((measurements{index+1}.omega - omega_b)*dt);
    F1 = F1';
    Fx = [F1  , -eye(3)*dt;
    zeros(3)  , eye(3)];
     P_ = Fx*P*Fx' + Q;
     % -----------------------observation---------------------
     % Prediction results and observations
     [H,detZ] = calH(quat,  measurements{index+1});
     % --------------------update-----------------
     K = P_*H'*inv(H*P_*H' + R)/2;
     err_state = K*detZ;
     P = P_ - K*(H*P_*H' + R)*K';
     % ----------------------update nominal state----------------------
     % Refer to previous functions，dt_theta-->quat,左乘
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
     % Why reset is needed？
     err_state = zeros(6,1);
     G = blkdiag(eye(3) - omegaMatrix(dt_theta/2) ,eye(3));
     P = G*P*G';
end

% figure;
% subplot(3,1,1)
% plot(pitchRef);
% hold on;plot(oula(:,2)/pi*180);
% subplot(3,1,2)
% plot(rollRef);
% hold on;plot(oula(:,1)/pi*180);
% subplot(3,1,3)
% plot(yawRef);
% hold on;plot(oula(:,3)/pi*180);
% legend 1 2 

 rotLim = [-5 5];
figure;
subplot(3,1,1)
plot(oula(:,1)/pi*180 - rollRef);
subplot(3,1,2)
plot(oula(:,2)/pi*180 - pitchRef);
subplot(3,1,3)
plot(oula(:,3)/pi*180 - yawRef);
% ylim(rotLim)
error = oula/pi*180  - [rollRef , pitchRef , yawRef];
fprintf('%4.3f,%4.3f,%4.3f\n',mean(error))

function R = q2R(q)
%四元数转旋转矩阵
R=[ 2*q(1).^2-1+2*q(2)^2    2*(q(2)*q(3)-q(1)*q(4)) 2*(q(2)*q(4)+q(1)*q(3));
    2*(q(2)*q(3)+q(1)*q(4)) 2*q(1)^2-1+2*q(3)^2     2*(q(3)*q(4)-q(1)*q(2));
    2*(q(2)*q(4)-q(1)*q(3)) 2*(q(3)*q(4)+q(1)*q(2)) 2*q(1)^2-1+2*q(4)^2];
R2 = R;
end

function Q_dt_theta = cal_Q_dt_theta(quat)
Q_dt_theta = 0.5* [-quat(2)  -quat(3)   -quat(4); ...
                quat(1)  -quat(4)    quat(3); ...
                quat(4)   quat(1)   -quat(2); ...
               -quat(3)   quat(2)    quat(1)]; 
end

function F = Exp_lee(in)
S = omegaMatrix(in);
    normV  = sqrt(S(1,2)^2+S(1,3)^2+S(1,3)^2);
    F = eye(3)+sin(normV)/normV*S(:,:)+...
            (1-cos(normV))/normV^2*S(:,:)^2;
end
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
function q = R2q(R)
%旋转矩阵转四元数
t=sqrt(1+R(1,1)+R(2,2)+R(3,3))/2;
q=[t (R(3,2)-R(2,3))/(4*t) (R(1,3)-R(3,1))/(4*t) (R(2,1)-R(1,2))/(4*t)];
Q1 = q;
end

function  q = oula2q(in)
x = in(1);
y = in(2);
z = in(3);
%欧拉角转四元数
q = [cos(x/2)*cos(y/2)*cos(z/2) + sin(x/2)*sin(y/2)*sin(z/2) ...
    sin(x/2)*cos(y/2)*cos(z/2) - cos(x/2)*sin(y/2)*sin(z/2) ...
    cos(x/2)*sin(y/2)*cos(z/2) + sin(x/2)*cos(y/2)*sin(z/2) ...
    cos(x/2)*cos(y/2)*sin(z/2) - sin(x/2)*sin(y/2)*cos(z/2)];

end
function Ang3 = q2oula(q)
%四元数转欧拉角
x = atan2(2*(q(1)*q(2)+q(3)*q(4)),1 - 2*(q(2)^2+q(3)^2));
y = asin(2*(q(1)*q(3) - q(2)*q(4)));
z = atan2(2*(q(1)*q(4)+q(2)*q(3)),1 - 2*(q(3)^2+q(4)^2));
Ang3 = [x y z];
end

function updateQuat = buildUpdateQuat(deltaTheta)
    deltaq = 0.5 * deltaTheta;
    updateQuat = [1; deltaq];
    updateQuat = updateQuat / norm(updateQuat);
end

function qLC = quatLeftComp(quat)
    vector = quat(2:4);
    scalar = quat(1);
    
    qLC = [  scalar ,  -vector';
             vector , scalar*eye(3) + crossMat(vector)  ];
end

function [H,detZ] = calH(q,measurements_k)
    % Normalise magnetometer measurement
    if(norm(measurements_k.mag) == 0), return; end	% 
    measurements_k.mag = measurements_k.mag / norm(measurements_k.mag);	% normalise magnitude,very important!!!!
    % Normalise accelerometer measurement
    if(norm(measurements_k.acc) == 0), return; end	% handle NaN
    measurements_k.acc  = measurements_k.acc / norm(measurements_k.acc);	% normalise accelerometer ,very important!!!!
    % Reference direction of Earth's magnetic feild
    h = quaternProd(q, quaternProd([0; measurements_k.mag], quatInv(q)));
    b = [0 norm([h(2) h(3)]) 0 h(4)];
    Ha = [2*q(3),                 	-2*q(4),                    2*q(1),                         -2*q(2)
         -2*q(2),                 	-2*q(1),                   -2*q(4),                         -2*q(3)
          0,                         4*q(2),                    4*q(3),                         0];
    Hm = [-2*b(4)*q(3),                2*b(4)*q(4),               -4*b(2)*q(3)-2*b(4)*q(1),       -4*b(2)*q(4)+2*b(4)*q(2)
          -2*b(2)*q(4)+2*b(4)*q(2),	   2*b(2)*q(3)+2*b(4)*q(1),    2*b(2)*q(2)+2*b(4)*q(4),       -2*b(2)*q(1)+2*b(4)*q(3)
           2*b(2)*q(3),                2*b(2)*q(4)-4*b(4)*q(2),	   2*b(2)*q(1)-4*b(4)*q(3),        2*b(2)*q(2)];
    Hx = [Ha, zeros(3,3)
          Hm, zeros(3,3)];
%     Hx = [Ha, zeros(3,3)];
    Q_detTheta  = [-q(2),    -q(3),      -q(4)
                    q(1),    -q(4),       q(3) 
                    q(4),     q(1),      -q(2) 
                   -q(3),     q(2),       q(1)];
    Xx = [0.5*Q_detTheta , zeros(4,3)
          zeros(3)       , eye(3)];
    H = Hx*Xx; 
    
    detZ_a = [ 2*(q(2)*q(4)  - q(1)*q(3)) + measurements_k.acc(1)
               2*(q(1)*q(2) + q(3)*q(4)) + measurements_k.acc(2)
               2*(0.5 - q(2)^2 - q(3)^2) + measurements_k.acc(3)];
%     detZ   = detZ_a;
    detZ_m =[((2*b(2)*(0.5 - q(3)^2 - q(4)^2) + 2*b(4)*(q(2)*q(4) - q(1)*q(3))) + measurements_k.mag(1))
             ((2*b(2)*(q(2)*q(3) - q(1)*q(4)) + 2*b(4)*(q(1)*q(2) + q(3)*q(4))) + measurements_k.mag(2))
             ((2*b(2)*(q(1)*q(3) + q(2)*q(4)) + 2*b(4)*(0.5 - q(2)^2 - q(3)^2)) + measurements_k.mag(3))]; 
    detZ   = [detZ_a;detZ_m];
end


function [roll,pitch,yaw] = quattoeuler(q)
rad2deg=180/pi;
T=[ 1 - 2 * (q(4) *q(4) + q(3) * q(3))  2 * (q(2) * q(3) +q(1) * q(4))         2 * (q(2) * q(4)-q(1) * q(3));
    2 * (q(2) * q(3)-q(1) * q(4))       1 - 2 * (q(4) *q(4) + q(2) * q(2))     2 * (q(3) * q(4)+q(1) * q(2));
    2 * (q(2) * q(4) +q(1) * q(3))      2 * (q(3) * q(4)-q(1) * q(2))          1 - 2 * (q(2) *q(2) + q(3) * q(3))];%cnb
roll  = atan2(T(2,3),T(3,3))*rad2deg;
pitch = asin(-T(1,3))*rad2deg;
yaw   = atan2(T(1,2),T(1,1))*rad2deg-8.3;%%这个固定偏差是什么鬼  
% yaw   = atan2(T(1,2),T(1,1))*rad2deg;%%这个固定偏差是什么鬼  
end