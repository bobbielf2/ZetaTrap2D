function u = HelmDLPn(ka,t,s,a,dens)
% normal derivative of Helmholtz double-layer potential interaction from 
% curve to target in 2D. No self-interaction correction.
%
% Unit test = run without arguments.
%
% Inputs
%   ka:     complex wavenumber
%   t:      target struct, t.x = target points
%   s:      source struct, output of setupquad.m
%   a:      source shift s.x -> s.x+a
%   dens:   density function. If not given (or empty), output matrix.
% Outputs
%   u:      potential value (if dens provided) or integral operator matrix
%

% Bowei Wu 12/25/2020

if nargin == 0, testHelmDLPn; return; end % unit test
if nargin >= 4 && ~isempty(a), s.x = s.x+a; end % source shift

d = bsxfun(@minus,t.x,s.x.'); % displacements matrix
r = abs(d);
csrx = conj(t.nx).*d;
csry = s.nx'.*d;
% build BIO matrix (note scaled by 4/i/ka)
A = ka * besselh(0,ka*r) .* real(csry).*real(csrx) ./ (r.*r) ...
       - besselh(1,ka*r) .* real(csry .* csrx) ./ (r.*r.*r);

if nargin < 5 || isempty(dens) %output matrix
    u = bsxfun(@times, A, .25i*ka * s.w(:).'); % speed wei & prefactor
else
    ttau = .25i*ka*dens.*s.w;
    u = A*ttau;
end


function testHelmDLPn    % unit test
% starfish curve
b = 0.3; m = 7; sh = randn();
r = @(t) 1 + b*cos(m*t+sh); rp = @(t) -b*m*sin(m*t+sh); 
rpp = @(t) -b*m^2*cos(m*t+sh);
s.Z = @(t) r(t).*exp(1i*t); s.Zp = @(t) exp(1i*t).*(1i*r(t)+rp(t)); 
s.Zpp = @(t) exp(1i*t).*(-r(t)+2i*rp(t)+rpp(t));
% setup quadrature
N = 200; s = setupquad(s, N);
% targ eval consistency check
ka = 20;                        % wavenumber
t.x = [s.x*0.5;s.x*1.5]; % targ inside & outside
t.nx = [s.nx;s.nx]; % targ normal
dens = cos(2*s.t+randn())+randn();  % rand dens func
p1 = HelmDLPn(ka,t,s)*dens; % generate full matrix then mat-vec
p2 = HelmDLPn(ka,t,s,[],dens); % mat-vec
fprintf('targ eval consistency err: %.15g\n',norm(p1-p2))


function i = diagind(A)
% DIAGIND  Return indices of diagonal of square matrix
% 
% Example usage:   A = randn(3,3); A(diagind(A)) = 0;
N = size(A,1); i = sub2ind([N,N], 1:N, 1:N); i = i(:);