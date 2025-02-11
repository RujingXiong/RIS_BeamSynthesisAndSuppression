% function [H] = Generate(R1, N, frequency, RIS_patch_position, d, theta1, theta2, phi1, phi2)
%     % 将角度转换为弧度
%     theta1 = deg2rad(theta1);
%     theta2 = deg2rad(theta2);
%     phi1 = deg2rad(phi1);
%     phi2 = deg2rad(phi2);
%     
%     % 初始化矩阵H
%     H = zeros(N, N);
%     lamda = 3e8 / frequency; % 波长
%     
%     % 遍历矩阵H的每个元素
%     for i = 1:N
%         for j = 1:N
%             aaa = RIS_patch_position(:, i); % 第i个单元的位置
%             bbb = RIS_patch_position(:, j); % 第j个单元的位置
%             
%             % 定义用于积分的函数句柄，直接使用矩阵乘法而不是 dot
%             y = @(theta, phi) (R1(i)' .* R1(j)) ./ ...
%                 (d + aaa(1) .*  sin(theta) .* cos(phi) + aaa(2) .* sin(theta) .* sin(phi) + aaa(3) .* cos(theta)) ./...
%                 (d + bbb(1) .*  sin(theta) .* cos(phi) + bbb(2) .* sin(theta) .* sin(phi) + bbb(3) .* cos(theta)) .* ...
%                 exp(1i * 2 * pi * ((bbb(1) - aaa(1)) .*  sin(theta) .* cos(phi) + (bbb(2) - aaa(2)) .* sin(theta) .* sin(phi) + (bbb(3) - aaa(3)) .* cos(theta)) / lamda);
%             
%             % 进行双重积分
%             H(i,j) = integral2(y, theta1, theta2, phi1, phi2);
%         end
%     end
% end
function [H] = Generate(R1, N, frequency, RIS_patch_position, d, theta1, theta2, phi1, phi2)
    % 将角度转换为弧度
    theta1 = deg2rad(theta1);
    theta2 = deg2rad(theta2);
    phi1 = deg2rad(phi1);
    phi2 = deg2rad(phi2);
    
    % 初始化矩阵H
    H = zeros(N, N);
    lamda = 3e8 / frequency; % 波长
    
    % 遍历矩阵H的每个元素
    for i = 1:N
        for j = 1:N
            aaa = RIS_patch_position(:, i); % 第i个单元的位置
            bbb = RIS_patch_position(:, j); % 第j个单元的位置
            for k = theta1 : 0.01 : theta2
                for l = phi1 : 0.01 : phi2
                    ccc = [sin(k)*cos(l),sin(k)*sin(l),cos(k)];
                    H(i,j) = H(i,j) + 0.0001 * (R1(i)' * R1(j)) / ...
                    (d + ccc*aaa) / (d + ccc*bbb) * exp(1i * 2 * pi * (ccc*(bbb-aaa)) / lamda);
                end
            end
        end
    end
end
