function A = HelmDLPnself(ka,s,ord)
% HELMDLPNSELF  self-interaction matrix for the target normal derivative of
% 2D Helmholtz double-layer potential. Use zeta correction for log singular
% components and finite-part zeta correction based on central difference
% approx for hypersingular components.
%
% Inputs
%   s:      source struct, output of setupquad.m
%   ord:    desired ord of zeta correction. Default ord = 16.
%           Use finite-part zeta correction
%           with central difference formula of order = ord-1.
% Outputs
%   A:      potential matrix

% BW Jan 24

if nargin < 3, ord = 16; end % default order

% native matrix
d = s.x-s.x.'; r = abs(d);  % displacements mat
csry = s.nx'.*d;            % (cos phi + i sin phi).r
csrx = conj(s.nx).*d;       % (cos th + i sin th).r
cxcy = real(csry).*real(csrx) ./ (r.*r);  	% cos phi cos th
cdor = -real(csry.*csrx) ./ (r.*r.*r);   	% -cos(phi+th) / r = (nxny - 2*cxcy)/r
H0 = triu(besselh(0,triu(ka*r,1)),1);       % symmetry, 2x speedup
H1 = triu(besselh(1,triu(ka*r,1)),1);
A = (H1+H1.').*cdor + ka*cxcy.*(H0+H0.');
A = bsxfun(@times, A, (.25i*ka)*s.w(:).');    % quadr wei

% prepare local param & indices for correction
N = size(A,1);
k = floor((ord-2)/2); % use (2*k+1)-pt correction
subi = (1:N)'.*ones(1,2*k+1);       % i of local stencils
subj = mod((0:N-1)'+(-k:k),N)+1;    % periodized global j of local stencils
ind = sub2ind([N,N],subi,subj);     % linear indexing
% d = s.x-s.x(ind2); r = abs(d); 	% local displacement d = x-y, r = |d|
r = r(ind); % local r
nxny = real(conj(s.nx).*s.nx(subj));    % local dot(nx,ny)
% cxcy=real(conj(s.nx).*d).*real(conj(s.nx(ind2)).*d)./r.^2; % local cos(phi)*cos(th)
% cdor = (nxny-2*cxcy)./r;   % local -cos(phi+th)/r

% singular correction
A(diagind(A)) = - pi/6 ./s.w + s.cur.^2/(4*pi).*s.w ...    	% finite-part zeta + diag const
    + ka^2*s.w.*(1i/8+(.5+psi(1)-log(ka/2*s.w))/(4*pi)) ; 	% H1 term diag const + log zeta corr const
% constant term correction via central difference
wcd = CDF2(k);	% 2nd deriv central difference coeffs
rpt = s.w.*(-k:k);                      % intrinsic metric in tangent space |y'*t|^2
D = (r.^2 - rpt.^2)./rpt.^2; D(:,k+1)=0;% rel deviation: (extrin - intrin)/intrin
C = (1-D+D.^2).*nxny./s.w.^2;
wcd = [wcd(end:-1:2),wcd]/(4*pi).*s.w(subj); % symm cd coeff, incl prefac & wei
A(ind) = A(ind) + wcd.*C;  % constant term 1/(2pi)*g''(0)*h/2
% now correct log kernel
wkr = kapur_rokhlin_sep_log(k+1); % K-R wei for separable log-singularity
J = ka*besselj(0,ka*r).*cxcy(ind)+besselj(1,ka*r).*cdor(ind);   % besselj smooth terms
J(:,k+1) = ka/2;
A(ind) = A(ind) + ka/(2*pi)*[wkr(end:-1:2),wkr].*s.w(subj).*J;  % zeta wei, incl prefac & wei


function i = diagind(A)
% DIAGIND  Return indices of diagonal of square matrix
%
% Example usage:   A = randn(3,3); A(diagind(A)) = 0;
N = size(A,1); i = sub2ind([N,N], 1:N, 1:N); i = i(:);
