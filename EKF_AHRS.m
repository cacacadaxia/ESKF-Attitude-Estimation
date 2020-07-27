%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extend Kalman Filter for attitude determination
% ʹ�� NED Ϊ��������ϵ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc;
clear;
close all;   
dataDir = '';
fileName = 'NAV_1'; 
Data = importdata(sprintf('%s.mat',fileName));
addpath('utils');
tic
gx=Data(:,27); %%����������
gy=Data(:,28); 
gz=Data(:,29); 
ax=Data(:,9);    %%���ٶȼ�����
ay=Data(:,10); 
az=Data(:,11); 
mx=Data(:,15); %%��ǿ������
my=Data(:,16); 
mz=Data(:,17); 
rad_deg=180/pi;   %�ɻ���ת���ɶ� 
deg_rad=pi/180;   %�ɶ�ת���ɻ���
L = size(Data,1);                  %��������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���������͹۲�������Ҫͨ���۲⾲ֵ̬����У��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wn_var  = 1e-6 * ones(1,4);               % rot vel var
wbn_var = 1e-8* ones(1,3);                % gyro bias change var
an_var  = 1e-1 * ones(1,3);               % acc var
mn_var  = 1e-1 * ones(1,3);               % mag var
Q = diag([wn_var, wbn_var]); 
R = diag([an_var, mn_var]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EKF����(AHRS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=9.87;
pitchEKF0=atan2(ax(1),sqrt(ay(1)*ay(1)+az(1)*az(1)));
%pitchEKF0=asin(ax(1));
%rollEKF0=asin(ay(1)/(g*cos(pitchEKF0)));
rollEKF0=atan2(-ay(1),-az(1));
%rollEKF0=atan2(ay(1),sqrt(ax(1)*ax(1)+az(1)*az(1)));
%yawEKF0=atan2(my(1)*cos(rollEKF0)+mz(1)*sin(rollEKF0),mx(1)*cos(pitchEKF0)+my(1)*sin(pitchEKF0)*sin(rollEKF0)-mz(1)*sin(rollEKF0)*cos(pitchEKF0));
yawEKF0=atan2(-my(1)*cos(rollEKF0)+mz(1)*sin(rollEKF0),mx(1)*cos(pitchEKF0)+my(1)*sin(pitchEKF0)*sin(rollEKF0)+mz(1)*sin(pitchEKF0)*cos(rollEKF0))-8.3*deg_rad;
[a,b,c,d]=quatfromeuler(rollEKF0,pitchEKF0,yawEKF0);
yawEKF(1)=yawEKF0*rad_deg;
pitchEKF(1)=pitchEKF0*rad_deg;
rollEKF(1)=rollEKF0*rad_deg;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����״̬����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = zeros(L,7);                              % x is estimated state, [q0 q1 q2 q3 bx by bz], NOTICE: every state is a row vector
Pk = zeros(7,7,L);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����۲����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

acc = Data(:,9:11);               % accelerator
gyro = Data(:,27:29);              % gyro
mag =Data(:,15:17);              % magnetometer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x(1,:) = [a b c d 0 0 0];           %Ԥ���ֵ
Pk(:,:,1) = eye(7);
x_ = zeros(1,7);
Time(1)=0;
for i=2:L-1
 Time(i)=Time(i-1)+0.02;
end

Ts=0.02;
t_e=zeros(1,L-1);
for k=1:L-1
    tic;
    acc(k,:)=acc(k,:)/norm(Data(k,9:11));
    mag(k,:)=mag(k,:)/norm(Data(k,15:17));
    z(k,:) = [acc(k,:),mag(k,:)];
    % ��k+1������ֵ/Ԥ��״̬����
    x_(1) = x(k,1) + Ts/2*(-x(k,2)*(gyro(k,1)-x(k,5))-x(k,3)*(gyro(k,2)-x(k,6))-x(k,4)*(gyro(k,3)-x(k,7)));
    x_(2) = x(k,2) + Ts/2*(x(k,1)*(gyro(k,1)-x(k,5))+x(k,3)*(gyro(k,3)-x(k,7))-x(k,4)*(gyro(k,2)-x(k,6)));
    x_(3) = x(k,3) + Ts/2*(x(k,1)*(gyro(k,2)-x(k,6))-x(k,2)*(gyro(k,3)-x(k,7))+x(k,4)*(gyro(k,1)-x(k,5)));
    x_(4) = x(k,4) + Ts/2*(x(k,1)*(gyro(k,3)-x(k,7))+x(k,2)*(gyro(k,2)-x(k,6))-x(k,3)*(gyro(k,1)-x(k,5)));
    x_(5) = x(k,5);
    x_(6) = x(k,6);
    x_(7) = x(k,7);
    x_(1:4) = x_(1:4)/norm(x_(1:4));
    %��k+1����������/Ԥ��Э�������
    Ak = eye(7)+Ts/2*...
         [0                    -(gyro(k,1))  -(gyro(k,2)) -(gyro(k,3)) x(k,2)  x(k,3)  x(k,4);
          (gyro(k,1)) 0                      (gyro(k,3))  -(gyro(k,2)) -x(k,1) x(k,4)  -x(k,3);
          (gyro(k,2)) -(gyro(k,3))  0                     (gyro(k,1))  -x(k,4) -x(k,1) x(k,2);
          (gyro(k,3)) (gyro(k,2))   -(gyro(k,1))  0                     x(k,3)  -x(k,2) -x(k,1);
          0 0 0 0 0 0 0;
          0 0 0 0 0 0 0;
          0 0 0 0 0 0 0];%dfdx
    Pk_ = Ak * Pk(:,:,k) * Ak' + Q;
    %��k+1��У������/���������棬�Ѳ�������h���ſ˱Ⱦ��������������Ⱦ���H
    
    hk_1 = [-(2*(x_(2)*x_(4)-x_(1)*x_(3))) -(2*(x_(1)*x_(2)+x_(3)*x_(4))) -(1-2*x_(2)^2-2*x_(3)^2)];%T3_WtoB������[0,0,-g]
     h = quaternProd(x_(1:4), quaternProd([0 mag(k,:)], quatInv(x_(1:4)')));
     b = [0 norm([h(2) h(3)]) 0 h(4)];
    hk_2 = [(2*b(2)*(0.5 - x_(3)^2 - x_(4)^2) + 2*b(4)*(x_(2)*x_(4) - x_(1)*x_(3))),2*b(2)*(x_(2)*x_(3) - x_(1)*x_(4)) + 2*b(4)*(x_(1)*x_(2) + x_(3)*x_(4)) ,2*b(2)*(x_(1)*x_(3) +x_(2)*x_(4)) + 2*b(4)*(0.5 - x_(2)^2 - x_(3)^2)];%T1_WtoB�ų�[m,0,0]
    hk = [hk_1 hk_2];%1*6
    Hk_1 = 2*[x_(3)  -x_(4) x_(1)  -x_(2) 0 0 0;
                -x_(2) -x_(1) -x_(4) -x_(3) 0 0 0;
                -x_(1) x_(2)  x_(3)  -x_(4) 0 0 0];
    Hk_2 = 2*[-2*b(4)*x_(3),               2*b(4)*x_(4),               -4*b(2)*x_(3)-2*b(4)*x_(1),       -4*b(2)*x_(4)+2*b(4)*x_(2) 0 0 0;
               -2*b(2)*x_(4)+2*b(4)*x_(2),	2*b(2)*x_(3)+2*b(4)*x_(1),	2*b(2)*x_(2)+2*b(4)*x_(4),       -2*b(2)*x_(1)+2*b(4)*x_(3) 0 0 0;
                2*b(2)*x_(3),                2*b(2)*x_(4)-4*b(4)*x_(2),	2*b(2)*x_(1)-4*b(4)*x_(3),        2*b(2)*x_(2)  0 0 0];
    Hk = [Hk_1;Hk_2];%6*7
    Kk = Pk_ * Hk' * (Hk * Pk_ * Hk' + R)^(-1);%����������7*6
    x(k+1,:) = (x_' + Kk * (z(k,:) - hk)')';%����״̬����
    x(k+1,1:4) = x(k+1,1:4)/norm(x(k+1,1:4));
    Pk(:,:,k+1) = (eye(7) - Kk*Hk) * Pk_;%����Э�������
    % x(k,1:4) = x(k,1:4)/norm(x(k,1:4));
     t_e(k)=toc;
   w(k)=x(k,1);xx(k)= x(k,2);y(k)=x(k,3);zz(k)=x(k,4);
   rollEKF(k)  = atan2(2 * (y(k) * zz(k) + xx(k) * w(k)) , 1 - 2 * (y(k) * y(k) + xx(k) * xx(k)))*rad_deg;
   pitchEKF(k) = asin(2 * (w(k) * y(k) - xx(k) * zz(k)))*rad_deg;
   yawEKF(k)   = atan2(2 * (xx(k) * y(k) + zz(k) * w(k)) , 1 - 2 * (zz(k) * zz(k) + y(k) * y(k)))*rad_deg-8.3;   
end
toc
t_e(L)=t_e(L-1);
Time(L)=Time(L-1)+0.02;
yawEKF(L)=yawEKF(L-1);
rollEKF(L)=rollEKF(L-1);
pitchEKF(L)=pitchEKF(L-1);
save t_e;

%% ==========================plot======================== %%
rollRef = Data(:,30)*57.3;
pitchRef = Data(:,31)*57.3;
yawRef = Data(:,32)*57.3;
figure;
subplot(3,1,1)
plot(rollEKF' - rollRef);
subplot(3,1,2)
plot(pitchEKF' - pitchRef);
subplot(3,1,3)
plot(yawEKF' - yawRef);
legend 1 2 



figure('Name', 'Roll Angle');
hold on;
plot(Time,rollEKF,'k',Time,Data(:,30)*57.3,'b') ;
grid on;
title('Roll Angle');
xlabel('Time (s)');
ylabel('Roll(deg)');
legend('EKF','ground truth');
hold off;

figure('Name', 'Pitch Angle');
hold on;
plot(Time,pitchEKF,'c',Time,Data(:,31)*57.3,'g') ;
grid on;
title('Pitch Angle');
xlabel('Time (s)');
ylabel('Pitch(deg)');
legend( 'EKF','ground truth');
hold off;

figure('Name', 'Yaw Angle');
hold on;
plot(Time,yawEKF,'m',Time,Data(:,32)*57.3,'r') ;
grid on;
title('Yaw Angle');
xlabel('Time (s)');
ylabel('Yaw(deg)');
legend('EKF','ground truth');
hold off;



