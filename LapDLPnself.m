function A = LapDLPnself(s,ord)
% LAPDLPNSELF  self-interaction matrix for the target normal derivative of 
% 2D Laplace double-layer potential. 
% Quadrature correction scheme either use Sidi's staggered-grid scheme with
% diagonal zeta correction (exp. convergence) or use my local zeta
% correction scheme that combines with central difference formulae.
% 
% Inputs
%   s:      source struct, output of setupquad.m
%   ord:    desired ord of zeta correction. Use finite-part zeta correction
%           with central difference formula of order = ord-1.
%           If ord not given, use the staggered grid correction (Sidi 2013)
% Outputs
%   A:      potential matrix

% Wu 1/23/21

[~,A] = LapDLPmat(s,s); % native matrix
if nargin < 2           % Sidi's staggered-grid zeta correction
    N = size(A,1);
    assert(mod(N,2)==0,'using Sidi''s staggered zeta correction: number of points must be even!')
    A(1:2:end,1:2:end) = 0; % staggered grid
    A(2:2:end,2:2:end) = 0;
    A = 2*A;                % weights*2 on staggered grid
    A(diagind(A)) = -pi/4 ./s.w;
else
    % hypersingular correction
    A(diagind(A)) = -pi/6 ./s.w + s.cur.^2.*s.w/(4*pi);
    % constant term correction
    k = floor((ord-2)/2);
    cd = CDF2(k); % 2nd deriv central difference coeffs
    N = size(A,1);
    ind1 = (1:N)'.*ones(1,2*k+1);
    ind2 = mod((0:N-1)'+(-k:k),N)+1;
    ind = sub2ind([N,N],ind1,ind2);
    r2 = abs(s.x - s.x(ind2)).^2;
    rpt = s.w.*(-k:k);
    D = (r2 - rpt.^2)./rpt.^2; D(:,k+1)=0;
    C = (1-D+D.^2).*real(conj(s.nx).*s.nx(ind2)).*s.w(ind2)./s.w.^2;
    A(ind) = A(ind) + cd([end:-1:2,1:end])/(4*pi).*C;
end

function i = diagind(A)
% DIAGIND  Return indices of diagonal of square matrix
% 
% Example usage:   A = randn(3,3); A(diagind(A)) = 0;
N = size(A,1); i = sub2ind([N,N], 1:N, 1:N); i = i(:);
