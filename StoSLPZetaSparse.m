function [U,P] = StoSLPZetaSparse(s,mu,ord)
% STOSLPZETASPARSE  self-interaction matrix for the 2D Stokes single-layer
% velocity and pressure potentials. (Traction potential is non-singular.)
% Quadrature correction uses zeta correction and combines with central
% difference formulae.
% 
% Inputs
%   s:      source struct, output of setupquad.m
%   mu:     viscosity
%   ord:    desired ord of zeta correction. Use finite-part zeta correction
%           with central difference formula of order = ord-1.
% Outputs
%   U:      velocity potential matrix
%   P:      pressure potential matrix

% BW Jun,2021

N = numel(s.x);
k = floor((ord-1)/2);
if k>=0
    ind1 = (1:N)'.*ones(1,2*k+1);       % stencil i-indices
    ind2 = mod((0:N-1)'+(-k:k),N)+1;    % stencil j-indices
    
    % velocity: log-singular correction
    kr = kapur_rokhlin_sep_log(k+1);        % kapur-rokhlin separable-log wei
    UW = [kr(end:-1:2),kr].*s.w(ind2);
    UW(:,k+1) = UW(:,k+1) - log(s.w).*s.w;  % diagonal limit
    U = sparse(ind1,ind2,UW,N,N);
    U = kron((1/(4*pi*mu))*eye(2),U);       % block 2x2 tensor kernel
    
    % pressure: singular correction
    if nargout > 1
        if isfield(s,'xpp') % more accurate when 2nd deriv info available
            cd = CDF1(k);           % 1st deriv central difference coeffs
            t = 2*pi/N*(-k:k);      % local parameter
            d = -s.xp-s.xpp.*t/2;   % leading components of (x(0)-x(t))/t = -x'(0) - x''(0)/2*t + O(t^2)
            D = real(s.xp.*conj(s.xpp))./s.sp.^2.*t; % leading O(t) component of (r^2 - (sp*t)^2)/(sp*t)^2
            if k > 0
                PW = (1-D).*d.*s.sp(ind2)./s.sp.^2.*([-cd(end:-1:1),0,cd]/(2*pi));
            else
                PW = zeros(size(D));
            end
        else
            cd = CDF2(k); % 2nd deriv central difference coeffs
            d = s.x - s.x(ind2);
            rpt = s.w.*(-k:k);
            D = (abs(d).^2 - rpt.^2)./rpt.^2; D(:,k+1)=0;
            PW = (1-D).*d.*s.w(ind2)./s.w.^2.*(cd([end:-1:2,1:end])/(4*pi));
        end
        P = sparse([ind1,ind1],[ind2,ind2+N],[real(PW),imag(PW)],N,2*N);% block 1x2 kernel
    end
else
    U = 0; P = 0;
end
