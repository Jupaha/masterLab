%% example how to use comp_eigenvalues()

% parameters for a simplyfied halogen lamp system
    R_c = 0.37;      % electrical resistance @ T_u in ohm
    L = 15e-6;       % inductance in H
    c_w = 7e-3;      % heat capacity of filament in J/K
    T_u = 295;       % ambient temperature in K
    b = 1.5e-12;     % coefficient to compute radiant power in W/K^4
    k = 0.9;         % exponent to compute electrical resistance at hot temperatures

% parameters for general equation form
    %p = zeros(1,6)
    p(1) = -R_c / (L * T_u * k);
    p(2) = k;
    p(3) = 1 / L;
    p(4) = -b / c_w;
    p(5) = T_u;
    p(6) = R_c / (c_w * T_u * k);

% define DGLs
    fct = @(t,x)[...
        p(1) * x(2) *p(2) * x(1);
        p(4) * (x(2) - p(5)) *4 + p(6) * x(2) * p(2) * x(1) *2];
    
% define ranges of x[i]
    r = [0, 100; 300, 10000];

% define steps in ranges
    n = 20;
    


%% compute eigenvalues

[eig1, eig2] =  comp_eigenvalues2(fct, r, n);
xe = gen_xe(n, r);
figure
surf(xe(1,:),xe(2,:),(eig1-eig2))
%hold on
%surf(xe(1,:),xe(2,:),eig2)
 
    if real([eig1, eig2]) < 0
        disp('System is stable')
        disp(['minimum step size FWE: ' num2str((2 / max(max(abs([eig1, eig2])))))])
    else
        disp('System is not stable')
    end

function [eig1, eig2] = comp_eigenvalues2(f, r, n)
    % compute eigenvalues of non linear DGL system"""

    syms t;
    x = sym('x',[1 length(r)]);
    sym(f(t, x));
    J = jacobian(sym(f(t, x)),x);
    A = matlabFunction(J,'Vars',x);
    
    % generate array with all combinations of x[i]
    xe = gen_xe(n, r);

    % solve A for all combinations of xe entrys and compute eigenvalues
    eig1 = zeros(n,n);
    eig2 = zeros(n,n);
    for i=1:n
        for j=1:n
            A_temp = A(xe(1,i), xe(2,j));
            eig_temp = eig(A_temp);
            eig1(i,j) = eig_temp(1);
            eig2(i,j) = eig_temp(2);
        end
    end
    end
function xs = gen_xe(n,r)
        %generate array with all combinations of x[i]
        for i=1:length(r)
            xs(i,:) = linspace(r(i,1), r(i,2), n);
        end
        end
   