function u = HelmSLPn(ka,t,s,a,dens,ord)
% normal derivative of Helmholtz single-layer potential interaction from 
% closed curve to target in 2D. Self-interaction via zeta-corrected
% trapezoidal rule [1]. 
%
% Unit test = run without arguments.
%
% Inputs
%   ka:     complex wavenumber
%   t:      target struct, t.x = target points, t.nx = target unit normal
%   s:      source struct, output of setupquad.m
%   a:      source shift s.x -> s.x+a
%   dens:   density function. If not given (or empty), output matrix.
%   ord:    order of zeta correction for self interaction (default 16)
% Outputs
%   u:      potential value (if dens provided) or integral operator matrix
%
% [1] Wu, B., & Martinsson, P.G. (2020). Zeta Correction: A New Approach 
%     to Constructing Corrected Trapezoidal Rule for Singular Integral 
%     Operators. (arxiv)

% Bowei Wu 6/7/2020

if nargin == 0, testHelmSLPn; return; end % unit test
if nargin > 3 && ~isempty(a), s.x = s.x+a; end % source shift
if nargin < 6 || isempty(ord), ord = 16; end

N = numel(s.x);
d = bsxfun(@minus,t.x,s.x.'); % displacements matrix
r = abs(d);

if (numel(s.x)==numel(t.x)) && max(abs(s.x-t.x))<1e-14 % self eval
    % build BIO matrix using symmetry
    A = triu(besselh(1,triu(ka*r,1)),1);
    A = -.25i*ka * (A + A.') .* real(conj(s.nx).*d)./r;
    
    % Use 2*K+1 local corr points,i.e. K on each side and 1 center point
    K = ceil((ord-2)/2);            % def K based on ord==2*K+2
    K = max(0,min(K,floor((N-1)/2)));
    c = kapur_rokhlin_sep_log(K+1); % K-R wei for separable log-singularity
    A(diagind(A)) = -s.cur/(4*pi);     	% diag limit
    % prepare loc corr indices & aux bessel func vals
    ind1 = repmat((1:N)',1,2*K);
    ind2 = mod(([-K:-1,1:K]) + (0:N-1)',N)+1;
    ind = sub2ind([N,N],ind1,ind2);
    rJ = r(ind);
    J = -ka*besselj(1,ka*rJ).*real(conj(s.nx).*d(ind))./rJ;
    if nargin < 5 || isempty(dens) %output matrix
        A(ind) = A(ind) + [c(K+1:-1:2),c(2:K+1)]/(2*pi).*J; % loc corr
        u = bsxfun(@times, A, s.w(:).'); % speed wei & prefactor
    else
        ttau = dens.*s.w;
        u = A*ttau; % trapz rule w\ diag corr
        circ = ttau(ind2).*J; % circulant (hankel) matrix consist of (-K:K) neighboring data for each pt
        u = u + circ * [c(K+1:-1:2),c(2:K+1)]'/(2*pi); % local correction
    end
else
    % build BIO matrix (note scaled by 2*pi)
    A = -.25i*ka * besselh(1,ka*r) .* real(conj(t.nx).*d)./r;
    if nargin < 5 || isempty(dens) %output matrix
        u = bsxfun(@times, A, s.w(:).'); % speed wei & prefactor
    else
        ttau = dens.*s.w;
        u = A*ttau;
    end
end

function testHelmSLPn   % unit test
% starfish curve
b = 0.3; m = 7; shft = randn();
r = @(t) 1 + b*cos(m*t+shft); rp = @(t) -b*m*sin(m*t+shft); rpp = @(t) -b*m^2*cos(m*t+shft);
s.Z = @(t) r(t).*exp(1i*t); s.Zp = @(t) exp(1i*t).*(1i*r(t)+rp(t)); s.Zpp = @(t) exp(1i*t).*(-r(t)+2i*rp(t)+rpp(t));
% setup quadrature
N = 200; s = setupquad(s, N);
% self eval consistency check
ka = 11+12i;                        % wavenumber
% ka = 20;                        % wavenumber
dens = cos(2*s.t+randn())+randn();  % rand dens func
p1 = HelmSLPn(ka,s,s)*dens; % generate full matrix then mat-vec
p2 = HelmSLPn(ka,s,s,[],dens); % mat-vec
fprintf('self eval consistency err: %.15g\n',norm(p1-p2))
% targ eval consistency check
t.x = [s.x*0.5;s.x*1.5]; % targ inside & outside
t.nx = [s.nx;s.nx];      % targ normal
p1 = HelmSLPn(ka,t,s)*dens; % generate full matrix then mat-vec
p2 = HelmSLPn(ka,t,s,[],dens); % mat-vec
fprintf('targ eval consistency err: %.15g\n',norm(p1-p2))

% self interaction convergence
fprintf('self interaction convergence...')
figure; subplot(1,2,1) % plot domain & target
plot(real(s.Z(0)),imag(s.Z(0)),'*k'); hold on; 
plot(s.Z(linspace(0,2*pi,200)),'k'); hold off
axis equal; legend({'target'})
mkr = {'o','^','v','+','*','s'}; % markers for plots
ORD = [2,6,10,16,32,42]; j = 1;
for ord = ORD % expected order of convergence is O(log(h)*h^(2k+1))
    a = randn()*3; sig = @(t) a*exp(cos(t+a));    
    NN = 40:40:800; I = [];
    for N = NN
        s = setupquad(s, N);  % set up quadr
        II = HelmSLPn(ka,s,s,[],sig(s.t),ord); % self eval
        I = [I,II(1)];
    end
    N = 1500; s = setupquad(s, N);
    II = HelmSLPn(ka,s,s,[],sig(s.t)); Iref = II(1);
    subplot(1,2,2)
    err = abs(I - Iref)/abs(Iref);
    semilogy(NN,err, mkr{j}); hold on; drawnow
    j=j+1;
end
legend(strcat('ord=',string(ORD)))
xlabel('N'); title('zeta-corrected trapz rule')
fprintf('done\n')

function i = diagind(A)
% DIAGIND  Return indices of diagonal of square matrix
% 
% Example usage:   A = randn(3,3); A(diagind(A)) = 0;
N = size(A,1); i = sub2ind([N,N], 1:N, 1:N); i = i(:);