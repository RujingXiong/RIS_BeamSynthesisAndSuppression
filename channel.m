function [R,R1,RIS_patch_position] = channel(Tx, Rx_position, frequency, RIS_row, RIS_col)
lamda = 3e8/frequency;
RIS_patch_position = zeros(3,RIS_row*RIS_col);
R = zeros(1,RIS_row*RIS_col);
R1 = zeros(1,RIS_row*RIS_col);
for k = 1:RIS_row
    for j =1:RIS_col
        RIS_patch_position(1,(k-1)*RIS_col+j) = lamda/2*(j-1-(RIS_col-1)/2);
        RIS_patch_position(2,(k-1)*RIS_col+j) = lamda/2*(k-1-(RIS_row-1)/2);
    end
end
aaa = length(Tx.d);
for l = 1 : aaa
    Tx_position = Tx.d(l).*[sin(Tx.theta(l))*cos(Tx.phi(l));sin(Tx.theta(l))*sin(Tx.phi(l));cos(Tx.theta(l))];
    R0 = zeros(1,RIS_row*RIS_col);
    for i = 1 : RIS_row*RIS_col
        d = norm(Tx_position-RIS_patch_position(:,i));
        R0(i) = 1e4/(d^2)*exp(1i*2*pi*d/lamda); % times 1e4 for calculation easily, The constant coefficient may need to be adjusted to avoid algorithm errors caused by MATLAB's numerical computations.
    end
    R1 = R1 + R0; % related to the Tx-RIS channels
end
D = norm(Rx_position);
for i = 1 : RIS_row*RIS_col
    d = D + (RIS_patch_position(:,i).')*Rx_position/D;
    R(i) = R1(i)/(d^2)*exp(1i*2*pi*d/lamda);
end


