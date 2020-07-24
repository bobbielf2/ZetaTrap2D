% TEST_HELM2D_BIE
% Test zeta-corrected trapezoidal rule for Helmholtz layer potentials
% on smooth geometries by solving the Helmholtz exterior Dirichlet and
% Neumann problem via combined field integral equation formulation:
%      	(Dirichlet ansatz) u = (D   - 1i*eta*S)*tau
%     	(Neumann ansatz)   u = (D*S - 1i*eta*S)*tau {altern. u=(D*S+S)*tau}

ka = 20+0i;     % wavenumber
ord = 32;       % pick desired convergence order of singular quad

% set up source geometry (starfish domain)
a = .3; w = 7; R = @(t) (1 + a*cos(w*t)); s.Z = @(t) R(t).*exp(1i*t);
Rp = @(t) -a*w*sin(w*t); s.Zp = @(t) exp(1i*t).*(Rp(t)+1i*R(t));

% generate (random) reference solution at test points
ns = 10;                          	% num of source points
s_ps.x = .4*exp(2i*pi*rand(ns,1));	% random source location
s_ps.w = 1;                         % dummy wei
den_source = randn(ns,1);           % random source densities
t.x = 1.5*exp(2i*pi*(1:20)/20).';  	% test points
uexact = HelmSLP(ka,t,s_ps,den_source);    % ref soln at test pts

% (optional) plot reference soln (real part), point sources, & test points
nx = 150; gx = ((1:nx)/nx*2-1)*1.6;         % set up plotting grid
[xx, yy] = meshgrid(gx); zz = (xx+1i*yy); fhom = nan*(1+1i)*zz;
N = 600; xx = s.Z((0:N-1)'*(2*pi/N));       % sample nodes on curve
ii = inpolygon(real(zz),imag(zz),real(xx),imag(xx)); p.x = zz(~ii); % targ pts outside
fhom(~ii) = HelmSLP(ka,p,s_ps,den_source);  % generate ext field
figure; imagesc(gx,gx,real(fhom)); colorbar; hold on;
fill(real(xx),imag(xx),'w'); title('solution (real part)')
plot(s_ps.x,'.','MarkerSize',20,'LineWidth',2)
plot(t.x,'dk','MarkerSize',8,'LineWidth',2); axis equal off

% BVPs convergence test
fprintf('Dirichlet problem...\n')
NN = 20:40:700;
errD = zeros(size(NN));
ieta = 1i*real(ka);
j = 1;
for N = NN
  s = setupquad(s, N);                            	% set up quadr
  f = HelmSLP(ka,s,s_ps,den_source);                % Dirichlet bdry condition
  Ad = HelmDLP(ka,s,s,[],ord); As = HelmSLP(ka,s,s,[],ord);
  A = eye(N)/2 + Ad - ieta * As;                	% construct CFIE operator A=0.5+D-i*eta*S
  tau = A\f;                                        % solve CFIE for tau
  u = HelmDLP(ka,t,s,tau)-ieta*HelmSLP(ka,t,s,tau);	% eval soln u = (D-i*eta*S)*tau
  errD(j) = max(abs(u-uexact));                     % err
  fprintf('N=%d, \terr = %.6g\n',N,errD(j)); j = j+1;
end
fprintf('Neumann problem...\n')
errN = zeros(size(NN));
j = 1;
for N = NN
  s = setupquad(s, N);                              % set up quadr
  g = HelmSLPn(ka,s,s_ps,den_source);               % Neumann bdry condition
  Ad = HelmSLPn(ka,s,s,[],ord); As = HelmSLP(ka,s,s,[],ord);
  A = (0.5*ieta-0.25)*eye(N)+Ad*(Ad-ieta*eye(N));  	% construct CFIE operator A=(0.5*i*eta-0.25)+D'^2-i*eta*D'
  tau = A\g;                                        % solve CFIE for tau
  u = HelmDLP(ka,t,s,As*tau)-ieta*HelmSLP(ka,t,s,tau); % eval soln u = (D*S-i*eta*S)*tau
  errN(j) = max(abs(u-uexact));                    	% err
  fprintf('N=%d, \terr = %.6g\n',N,errN(j)); j = j+1;
end

figure; semilogy(NN,errD,'*',NN,errN,'d')
xlabel('N'); ylabel('err')
title(['Helmholtz BVP convergence, $\kappa=',num2str(ka),'$'],'interpreter','latex')
legend({'Dirichlet','Neumann'},'interpreter','latex')