function updateQuat = buildUpdateQuat(deltaTheta)
%
% Builds the update quaternion from the minimally parametrized update
% See Indirect Kalman Filter for 3D Attitude Estimation (Roumeliotis)
%
%% 从最小参数化更新建立更新四元数参见间接卡尔曼滤波器进行三维姿态估计（Roumeliotis）

    deltaq = 0.5 * deltaTheta;
    
    checkNorm = deltaq' * deltaq;
    
    if checkNorm > 1
        updateQuat = [1; deltaq];
        updateQuat = updateQuat / sqrt(1 + checkNorm);
    else
        updateQuat = [sqrt(1 - checkNorm);deltaq];
    end
    
    updateQuat = updateQuat / norm(updateQuat);
end