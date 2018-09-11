%% Halogen ODE
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
fct = @(t,x,u)[...
    p(1) * x(2) *p(2) * x(1);
    p(4) * (x(2) - p(5)) *4 + p(6) * x(2) * p(2) * x(1) *2];
% Compute all eigenvalues in range
r = [0, 100; 300, 10000];
n = 10; 
eigV = comp_eigenvalues(fct,r,n);


%%  Plot stability region of RK4 
clearvars -except eigV h_eigmax
clear i;
[X,Y] = meshgrid(-10:0.005:10,-10:0.005:10);
m = X+i*Y;

% stability function R see: Stability of Runge-Kutta Methods
R = (m.^4)/24 + (m.^3)/6 + (m.^2)/2 + m + 1;
%R = m.^7./24000 + (161.*m.^6)./120000 + (195398609418584047.*m.^5)./21374506043965440000 + (118747255799808017.*m.^4)./2849934139195392000 + (419573637159321617.*m.^3)./2517441822955929600 + (559431516212428817.*m.^2)./1118863032424857600 + (2983634753132953651.*m)./2983634753132953600 + 1;

Rhat = abs(R);
R_bool = Rhat <= 1;
%%
contour(X,Y,R_bool,'-m')
hold on

%plot(real(eigV), imag(eigV), '*')
%%
% x = -5;
% y = 7;
% mu_min = x+y*i;
% plot(x, y, 'or', 'LineWidth',10)
% line([x 0], [y 0]);  %line
% F_bool = abs(X+i*Y)<=abs(mu_min);
% contour(X,Y,F_bool,'-r')
% a_sqr = (abs(x)-abs(real(scaled_mu)))^2;
% b_sqr = (abs(y)-abs(imag(scaled_mu)))^2;
% r = sqrt((a_sqr)+(b_sqr));
% th = 0:pi/50:2*pi;
% xunit = r * cos(th) + x;
% yunit = r * sin(th) + y;
% h = plot(xunit, yunit);
%%
% Plot axis
xL = xlim;
yL = ylim;
line([0 0], yL);  %x-axis
line(xL, [0 0]);  %y-axis
grid


function [eigV] = comp_eigenvalues(f, r, n)
    % compute eigenvalues of non linear DGL system"""

    syms t;
    x = sym('x',[1 size(r,1)]);
    Jf = jacobian(sym(f(t, x)),x);
    %Jf = matlabFunction(Jf,'Vars',x);
    xe = gen_xe(n,r);
    N = size(r,1);

    % solve A for all combinations of xe entrys and compute eigenvalues
    for i=1:size(xe,1)
        A_temp = Jf;
        for j=1:N
            A_temp = subs(A_temp, x(j), xe(i,j));
        end
        eigV(i,:) = eig(double(A_temp));
%         spec_large = eigs(double(A_temp), 1, 'largestabs');
%         spec_small = eigs(double(A_temp), 1, 'smallestabs');
%         eigV(i,:) = [spec_large spec_small];
    end
end
function xe = gen_xe(n,r)
    %generate array with all combinations of x[i]
    for i=1:size(r,1)
        xe(i,:) = linspace(r(i,1), r(i,2), n);
    end
    xe = cartprod(xe);
end 
function X = cartprod(matrixIn)
%CARTPROD Cartesian product of multiple sets.
%
%   X = CARTPROD(A,B,C,...) returns the cartesian product of the sets 
%   A,B,C, etc, where A,B,C, are numerical vectors.  
%
%   Example: A = [-1 -3 -5];   B = [10 11];   C = [0 1];
% 
%   X = cartprod(A,B,C)
%   X =
% 
%     -5    10     0
%     -3    10     0
%     -1    10     0
%     -5    11     0
%     -3    11     0
%     -1    11     0
%     -5    10     1
%     -3    10     1
%     -1    10     1
%     -5    11     1
%     -3    11     1
%     -1    11     1
%
%   This function requires IND2SUBVECT, also available (I hope) on the MathWorks 
%   File Exchange site.

varargin = num2cell(matrixIn.',1);

numSets = size(matrixIn,1);
for i = 1:numSets
    thisSet = sort(varargin{i});
    if ~isequal(prod(size(thisSet)),length(thisSet))
        error('All inputs must be vectors.')
    end
    if ~isnumeric(thisSet)
        error('All inputs must be numeric.')
    end
    if ~isequal(thisSet,unique(thisSet))
        error(['Input set' ' ' num2str(i) ' ' 'contains duplicated elements.'])
    end
    sizeThisSet(i) = length(thisSet);
    varargin{i} = thisSet;
end

X = zeros(prod(sizeThisSet),numSets);
for i = 1:size(X,1)
    
    % Envision imaginary n-d array with dimension "sizeThisSet" ...
    % = length(varargin{1}) x length(varargin{2}) x ...
    
    ixVect = ind2subVect(sizeThisSet,i);
    
    for j = 1:numSets
        X(i,j) = varargin{j}(ixVect(j));
    end
end
end
function X = ind2subVect(siz,ndx)
%IND2SUBVECT Multiple subscripts from linear index.
%   IND2SUBVECT is used to determine the equivalent subscript values
%   corresponding to a given single index into an array.
%
%   X = IND2SUBVECT(SIZ,IND) returns the matrix X = [I J] containing the
%   equivalent row and column subscripts corresponding to the index
%   matrix IND for a matrix of size SIZ.  
%
%   For N-D arrays, X = IND2SUBVECT(SIZ,IND) returns matrix X = [I J K ...]
%   containing the equivalent N-D array subscripts equivalent to IND for 
%   an array of size SIZ.
%
%   See also IND2SUB.  (IND2SUBVECT makes a one-line change to IND2SUB so as
%   to return a vector of N indices rather than retuning N individual
%   variables.)%IND2SUBVECT Multiple subscripts from linear index.
%   IND2SUBVECT is used to determine the equivalent subscript values
%   corresponding to a given single index into an array.
%
%   X = IND2SUBVECT(SIZ,IND) returns the matrix X = [I J] containing the
%   equivalent row and column subscripts corresponding to the index
%   matrix IND for a matrix of size SIZ.  
%
%   For N-D arrays, X = IND2SUBVECT(SIZ,IND) returns matrix X = [I J K ...]
%   containing the equivalent N-D array subscripts equivalent to IND for 
%   an array of size SIZ.
%
%   See also IND2SUB.  (IND2SUBVECT makes a one-line change to IND2SUB so as
%   to return a vector of N indices rather than returning N individual
%   variables.)
 

% All MathWorks' code from IND2SUB, except as noted:

n = length(siz);
k = [1 cumprod(siz(1:end-1))];
ndx = ndx - 1;
for i = n:-1:1
  X(i) = floor(ndx/k(i))+1;      % replaced "varargout{i}" with "X(i)"
  ndx = rem(ndx,k(i));
end
end

