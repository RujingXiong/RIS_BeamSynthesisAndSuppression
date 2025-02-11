function [result, Coeff] = TestMaxMin(Coeff, omega0)
lambd = 10; % \lambda is important
deta = 1e-4;
sigma = 1e-4;
max_iter = 1e3;
Lip = .1;
K = Coeff.K;
T = Coeff.T;
N = Coeff.N;
Coeff1 = Coeff;
Coeff1.K = K + T;
Coeff1.D = ones(K+T,1);
t1 = -100000;% The value is related to the number of RIS units, if too large, it may cause f_1>0, f_1*f_2>0, leading to an algorithm error.

t2 = 1e10;
totalIter = 0;

operation1 = @(x) x + t1/N.*eye(N,N);
aaa1 = cellfun(operation1, Coeff.h, 'UniformOutput', false);
Coeff1.h = cat(1,aaa1,Coeff.H);
[omega, ~, ~, J, iter] = minimax_grad_unc(omega0, lambd, max_iter, Coeff1, Lip);
f1 = min(J) - 1/(4*lambd);

operation2 = @(x) x + t2/N.*eye(N,N);
aaa2 = cellfun(operation2, Coeff.h, 'UniformOutput', false);
Coeff1.h = cat(1,aaa2,Coeff.H);
[~, ~, ~, J, ~] = minimax_grad_unc(omega0, lambd, max_iter, Coeff1, Lip);
f2 = min(J) - 1/(4*lambd);

if f1*f2 > 0
    disp('错误')
else
    for loop = 1:100
        t  = (t1 +t2) / 2;
        operation = @(x) x + t/N.*eye(N,N);
        aaa = cellfun(operation, Coeff.h, 'UniformOutput', false);
        Coeff1.h = cat(1,aaa,Coeff.H);
        [omega, ~, ~, J, ~] = minimax_grad_unc(omega, lambd, max_iter, Coeff1, Lip);
        f = min(J) - 1/(4*lambd);
        f
        if f  >= 0
            t2 = t;
            totalIter = totalIter+iter;
        % elseif f <= 0
        elseif f <= -1/(2*lambd)
            t1 = t;
        elseif lambd > 1/(2*deta)
            break
        else
            lambd = 2*lambd;
        end
        if abs(t1-t2) < sigma
            break
        end
    end
end
result.omega = omega;
%result.J = J;
result.Iter = totalIter;