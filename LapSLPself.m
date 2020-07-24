function u = LapSLPself(s,ord,dens)
% LAPSLPSELF  2D Laplace single-layer potential self interaction from curve
% to itself using the zeta-corrected trapezoidal rule [1]. 
%
% Unit test = run without arguments.
%
% Inputs
%   s:      src seg.
%   ord:    order of zeta correction (default 16)
%   dens:   density function. If not given (or empty), output matrix.
% Outputs
%   u:      potential value (if dens provided) or integral operator matrix
% 
% [1] Wu, B., & Martinsson, P.G. (2020). Zeta Correction: A New Approach 
%     to Constructing Corrected Trapezoidal Rule for Singular Integral 
%     Operators. (arxiv)

% Bowei Wu 6/4/20

if nargin == 0, testLapSLP; return; end
if nargin < 2 || isempty(ord), ord = 16; end % default 16-th order

N = numel(s.x);
d = bsxfun(@minus,s.x,s.x.'); 	% C-# displacements mat
A = -log(abs(d));               % peri log
A(diagind(A)) = -log(s.w);      % diagonal limit

% Use 2*K+1 local corr points,i.e. 1 center point plus K pts on each side
K = max(0,ceil((ord-2)/2)); % def K based on ord==2*K+2
K = min(K,floor((N-1)/2));
c = kapur_rokhlin_sep_log(K+1); % kapur-rokhlin separable-log wei

if nargin < 3 || isempty(dens) % matrix only
    %A = A + toeplitz([c,zeros(1,N-2*K-1),c(K+1:-1:2)]);
    ind1 = (1:N)'.*ones(1,2*K+1);
    ind2 = mod((-K:K) + (0:N-1)',N)+1;
    ind = sub2ind([N,N],ind1,ind2);
    A(ind) = A(ind) + [c(K+1:-1:2),c];
    u = bsxfun(@times, A, (1/2/pi)*s.w(:)');   % prefactor & wei
else % mat-vec
    tau = (1/2/pi)*dens.*s.w(:); % smooth dens*wei, incl prefactor
    u = A*tau;
    % (2K+1)-point correction
    circ = hankel(tau([N-K+1:N,1:N-K]),tau([N-K:N,1:K])); % circulant matrix consist of (-K:K) neighboring data for each pt
    u = u + circ * [c(K+1:-1:2),c(1:K+1)]'; % correction
end

function testLapSLP
% starfish domain
a = .3; w = 7; sh = randn()*3;
R = @(t) (1 + a*cos(w*t+sh)); s.Z = @(t) R(t).*exp(1i*t);
% setup quadrature
N = 200; s = setupquad(s, N);
% consistency check
dens = cos(s.t);
p1 = LapSLPself(s)*dens; % generate full matrix then mat-vec
p2 = LapSLPself(s,[],dens); % mat-vec
fprintf('consistency err: %.15g\n',norm(p1-p2))

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
    NN = [50:50:1000,1500]; I = [];
    for N = NN
        s = setupquad(s, N);  % set up quadr
        II = LapSLPself(s,ord,sig(s.t)); % self eval
        I = [I,II(1)];
    end
    subplot(1,2,2)
    NN = NN(1:end-1); err = abs(I(1:end-1) - I(end))/abs(I(end));
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