function [strength] = cal(w,TX,frequency,d, RIS_row, RIS_col)
xita = -90:0.25:90;
xita = xita/180*pi;
L = length(xita);
strength = zeros(1,L);
for i = 1:L
    [R,~,~] = channel(TX, d.*[sin(xita(i));0;cos(xita(i))], frequency, RIS_row, RIS_col);
    strength(1,i) = w'*(R')*R*w;  % recover the strength while it has been amplificated in the channel matrix
end
