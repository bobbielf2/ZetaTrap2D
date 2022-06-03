function [P,T] = StoDLPZetaSparse(s,mu,ord)
% STODLPZETASPARSE  self-interaction matrix for the 2D Stokes double-layer
% pressure and traction potentials. (Velocity potential is non-singular.)
% Quadrature correction uses zeta correction and combines with central
% difference formulae.
% 
% Inputs
%   s:      source struct, output of setupquad.m
%   mu:     viscosity
%   ord:    desired ord of zeta correction. Use finite-part zeta correction
%           with central difference formula of order = ord-1.
% Outputs
%   P:      pressure potential matrix
%   T:      traction potential matrix

% BW Jun,2021

N = numel(s.x);
k = ceil((ord-1)/2);
if k > 0    % hypersinuglar & constant term corrections
    ind1 = (1:N)'.*ones(1,2*k+1);       % stencil i-indices
    ind2 = mod((0:N-1)'+(-k:k),N)+1;    % stencil j-indices
    %ind = sub2ind([N,N],ind1,ind2);    % stencil linear indices
    cd = CDF2(k); % 2nd deriv central difference coeffs
    
    % pressure: hypersingular & singular correction
    ny = s.nx(ind2);
    d = s.x - s.x(ind2);    % displacement
    r2 = real(d.*conj(d));
    dnyirr = real(d.*conj(ny))./r2; % (x-y).ny/r^2
    dnyirr(:,k+1)=-s.cur/2; % diagonal limits, will be zero out in DLP pressure, but useful in DLP traction
    rpt = s.w.*(-k:k);
    D = (r2 - rpt.^2)./rpt.^2; D(:,k+1)=0;
    Er2cd = (mu/pi/2) * cd([end:-1:2,1:end]) .* (1-D+D.^2).*s.w(ind2)./s.w.^2; % prefactor, CDF weights, and expansion of s.w/r^2 
    PW = Er2cd.*(2*dnyirr.*d-ny);
    PW(:,k+1) = PW(:,k+1) + mu*pi/3*s.nx./s.w;  % add hypersingular part
    P = sparse([ind1,ind1],[ind2,ind2+N],[real(PW),imag(PW)],N,2*N);    % block 1x2 kernel
    
    % traction: hypersingular & constant term corrections
    if nargout > 1
        nx1 = real(s.nx); nx2 = imag(s.nx); % normal & tangent vectors on the curve
        ny1 = real(ny);   ny2 = imag(ny);
        t1  = -nx2;       t2  = nx1;
        d1  = real(d);    d2  = imag(d);
        
        % compute constants from singular terms via central difference
        dnxirr = real(d.*conj(s.nx))./r2; % (x-y).nx/r^2
        dnxirr(:,k+1)=s.cur/2;  % diagonal limits
        
        % numerator of singular terms (as tensors)
        TW = cell(2,2);
        TW{1,1} = dnxirr.*ny1.*d1 + dnyirr.*d1.*nx1;
        TW{1,2} = dnxirr.*ny1.*d2 + dnyirr.*d1.*nx2;
        TW{2,1} = dnxirr.*ny2.*d1 + dnyirr.*d2.*nx1;
        TW{2,2} = dnxirr.*ny2.*d2 + dnyirr.*d2.*nx2;
        
        % numerator of hypersingular terms
        nxny = nx1.*ny1+nx2.*ny2;      	% nx.ny
        Tt11 = d1.^2 ./r2; Tt11(:,k+1) = t1.^2;     % (d x d)/r^2
        Tt12 = d1.*d2./r2; Tt12(:,k+1) = t1.*t2;
        Tt22 = d2.^2 ./r2; Tt22(:,k+1) = t2.^2;
        TW{1,1} = TW{1,1} + nxny.*Tt11 + nx1.*ny1;
        TW{1,2} = TW{1,2} + nxny.*Tt12 + nx1.*ny2;
        TW{2,1} = TW{2,1} + nxny.*Tt12 + nx2.*ny1;
        TW{2,2} = TW{2,2} + nxny.*Tt22 + nx2.*ny2;
        
        % multiply by prefactors and CDF weights
        TW{1,1} = Er2cd.*TW{1,1};
        TW{1,2} = Er2cd.*TW{1,2};
        TW{2,1} = Er2cd.*TW{2,1};
        TW{2,2} = Er2cd.*TW{2,2};
        
        % correct hypersingular terms
        c2 = -mu*pi/3 ./s.w;
        TW{1,1}(:,k+1) = TW{1,1}(:,k+1) + c2;
        TW{2,2}(:,k+1) = TW{2,2}(:,k+1) + c2;
        
        % compute constants from regular terms
        c0 = -s.cur.^2 .*s.w * (mu/pi/4); % diagonal limit (-s.cur^2/4)*s.w*(mu/pi)
        TW{1,1}(:,k+1) = TW{1,1}(:,k+1) + c0.*(1-8*Tt11(:,k+1));
        TW{1,2}(:,k+1) = TW{1,2}(:,k+1) - 8*c0.*Tt12(:,k+1);
        TW{2,1}(:,k+1) = TW{2,1}(:,k+1) - 8*c0.*Tt12(:,k+1);
        TW{2,2}(:,k+1) = TW{2,2}(:,k+1) + c0.*(1-8*Tt22(:,k+1));
        
        % form sparse matrix
        T = [sparse(ind1,ind2,TW{1,1},N,N), sparse(ind1,ind2,TW{1,2},N,N);
             sparse(ind1,ind2,TW{2,1},N,N), sparse(ind1,ind2,TW{2,2},N,N)];
        
    end
elseif k > -1   % only hypersinuglar term correction when ord < 1
    PW = spdiags(mu*pi/3*s.nx./s.w,0,N,N);
    P = [real(PW),imag(PW)];     
    TW = spdiags(-mu*pi/3 ./s.w,0,N,N);
    T = kron(eye(2),TW);
else
    P = 0; T = 0;
end    
