% TEST_LAP2D_HYPER_BIE
% Test hypersingular zeta-corrected trapezoidal rule for Laplace layer
% potentials on smooth geometries by solving the BVPs using a direct
% approach to BIE:
%      	int Laplace ansatz: u = S*(du/dn) - D*u
%       int Calderón projection:     u =  (1/2-D)*u  +         S*(du/dn)  
%                                du/dn =       -T*u  + (1/2+D^*)*(du/dn)
%      	ext Laplace ansatz: u = D*u - S*(du/dn) + omega
%       ext Calderón projection:     u =  (1/2+D)*u  -         S*(du/dn)
%                                du/dn =        T*u  + (1/2-D^*)*(du/dn)
%
% c.f. Hsiao-Wendland 2008, Sec.1.3-1.4


clear; clc; clf
ord = 32;       % pick desired convergence order of singular quad

% set up source geometry (starfish domain)
a = .3; w = 5; R = @(t) (1 + a*cos(w*t)); s.Z = @(t) R(t).*exp(1i*t);
Rp = @(t) -a*w*sin(w*t); s.Zp = @(t) exp(1i*t).*(Rp(t)+1i*R(t));

%% Interior Laplace BVPs
% generate (random) reference solution at test points
ns = 10;                          	% num of source points
s_ps.x = 1.5*exp(2i*pi*rand(ns,1));	% random source location
s_ps.w = 1;                         % dummy wei
den_source = randn(ns,1);           % random source densities
t.x = .4*exp(2i*pi*(1:20)/20).';  	% test points
uexac=LapSLPmat(t,s_ps)*den_source; % ref soln at test pts

% (optional) plot reference soln (real part), point sources, & test points
nx = 150; gx = ((1:nx)/nx*2-1)*1.6;         % set up plotting grid
[xx, yy] = meshgrid(gx); zz = (xx+1i*yy); fhom = nan(size(zz));
N = 600; xx = s.Z((0:N-1)'*(2*pi/N));       % sample nodes on curve
ii = inpolygon(real(zz),imag(zz),real(xx),imag(xx)); p.x = zz(ii); % targ pts outside
fhom(ii) = LapSLPmat(p,s_ps)*den_source;  % generate ext field
figure(1); subplot(1,4,1); cla
imagesc(gx,gx,fhom); colorbar; hold on;
fillout(real(xx),imag(xx),[-1,1,-1,1]*1.6,'w'); 
%title('interior: $$u(x) = S\left[\frac{\partial u}{\partial n}\right](x)-D[u](x)$$','interpreter','latex')
plot(s_ps.x,'*','MarkerSize',5,'LineWidth',1)
plot(t.x,'.k','MarkerSize',10,'LineWidth',1); axis('equal','off',[-1,1,-1,1]*1.6)

% int Diri BIE: (-1/2+D^*)*(du/dn)=T*u, uniquely solvable
fprintf('Dirichlet problem...\n')
NN = 20:20:400;
errD = zeros(size(NN)); err_un = zeros(size(NN)); errD2 = zeros(size(NN));
j = 1;
for N = NN                       	% convergence test
  s = setupquad(s, N);          	% set up quadr
  [f,gexac] = LapSLPmat(s,s_ps);
  f = f*den_source;	% Dirichlet bdry condition
  gexac = gexac*den_source;% exact Neu data
  if 1 % direct approach
      At = LapDLPnself(s,ord);        	% hypersingular T op
      [~,Adt] = LapSLPmat(s,s);         % D^* operator
      A = -eye(N)/2 + Adt;           	% construct A=-0.5+D^*
      g = A\(At*f);                 	% solve for g = du/dn
      u=LapSLPmat(t,s)*g-LapDLPmat(t,s)*f; % eval soln u = S*(du/dn)-D*u
      
      % Redo using Sidi's staggered grid finite-part correction
      At2 = LapDLPnself(s);            	% Sidi's corrected T op
      g2 = A\(At2*f);                	% solve for g = du/dn
      u2=LapSLPmat(t,s)*g2-LapDLPmat(t,s)*f; % eval soln u = S*(du/dn)-D*u
      errD2(j) = max(abs(u2-uexac));   	% err at test pts
  else % indirect approach (-0.5+D)*tau=f
      Ad = LapDLPmat(s,s);
      A = -eye(N)/2 + Ad;
      tau = A\f;
      u = LapDLPmat(t,s)*tau;
      g = LapDLPnself(s,ord)*tau;
  end
  errD(j) = max(abs(u-uexac));    	% err at test pts
  err_un(j) = max(abs(gexac-g));    	% err Neu data
  fprintf('N=%d, \tpotential err = %.6g, \tNeu data err = %.6g\n', ...
            N,errD(j),err_un(j)); j = j+1;
end
%subplot(2,2,2); semilogy(NN,errD,'*',NN,err_un,'d'); legend('test pts','Neu data'); title('error')

% int Neu BIE: (1/2+D)*u=S*(du/dn), unique up to const
fprintf('Neumann problem...\n')
errN = zeros(size(NN)); err_u = zeros(size(NN));
j = 1;
for N = NN                        	% convergence test
  s = setupquad(s, N);             	% set up quadr
  [fexac,g] = LapSLPmat(s,s_ps);
  g = g*den_source;                 % Neumann bdry condition
  fexac = fexac*den_source;         % exact Diri data
  if 1 % direct
      Ad = LapDLPmat(s,s);       	% D operator
      As = LapSLPself(s,ord);       % S operator
      A = eye(N)/2+Ad+s.w.';        % construct operator A=0.5+D (add 1 to kill 1d nullspace)
      %warning off MATLAB:nearlySingularMatrix
      f = A\(As*g);                 % solve for f=u (rank-1 deficient but okay)
      u=LapSLPmat(t,s)*g-LapDLPmat(t,s)*f; % eval soln u = S*(du/dn)-D*u
  else % indirect
      Adt = LapSLPnmat(s,s);        % D^* operator
      A = eye(N)/2+Adt+s.w.';
      tau = A\g;
      u = LapSLPmat(t,s)*tau;
      f = LapSLPself(s,ord)*tau;
  end
  offset = uexac(1)-u(1);         	% recover constant
  u=u+offset; f=f+offset;           % constant shift
  errN(j) = max(abs(u-uexac));      % err at test pts
  err_u(j) = max(abs(f-fexac));     % err Diri data
  fprintf('N=%d, \terr = %.6g, \tDiri data err = %.6g\n',...
      N,errN(j),err_u(j)); j = j+1;
end

subplot(1,4,2); semilogy(NN,errD,'*-',NN,errN,'d-',NN,errD2,'v-')
xlabel('N'); ylabel('err')
yticks([1e-16,1e-8,1e0])
ylim([1e-16,1e0])
%title('Laplace interior BVP','interpreter','latex')
legend({'Dirichlet (32nd)','Neumann (32nd)','Dirichlet (spec)'},'interpreter','latex')

%% exterior Laplace BVPs
% generate (random) reference solution at test points
ns = 10;                          	% num of source points
s_ps.x = .4*exp(2i*pi*rand(ns,1));	% random source location
s_ps.w = 1;                         % dummy wei
den_source = randn(ns,1);           % random source densities
t.x = 1.5*exp(2i*pi*(1:20)/20).';  	% test points
uexac=LapSLPmat(t,s_ps)*den_source;	% ref soln at test pts

% (optional) plot reference soln (real part), point sources, & test points
nx = 150; gx = ((1:nx)/nx*2-1)*1.6;         % set up plotting grid
[xx, yy] = meshgrid(gx); zz = (xx+1i*yy); fhom = nan(size(zz));
N = 600; xx = s.Z((0:N-1)'*(2*pi/N));       % sample nodes on curve
ii = inpolygon(real(zz),imag(zz),real(xx),imag(xx)); p.x = zz(~ii); % targ pts outside
fhom(~ii) = LapSLPmat(p,s_ps)*den_source;  % generate ext field
figure(1); subplot(1,4,3); cla
imagesc(gx,gx,fhom); colorbar; hold on;
fill(real(xx),imag(xx),'w'); 
%title('exterior: $$u(x) = D[u](x) - S\left[\frac{\partial u}{\partial n}\right](x)+c$$','interpreter','latex')
plot(s_ps.x,'*','MarkerSize',5,'LineWidth',1)
plot(t.x,'.k','MarkerSize',10,'LineWidth',1); axis('equal','off',[-1,1,-1,1]*1.6)

% ext Diri BIE: (1/2+D^*)*(du/dn)=T*u, uniquely solvable after modification
% growth condition: u->E*\log|x| as |x|->inf, 
% where total interior charge E==integral(du/dn) is given
fprintf('ext Dirichlet problem...\n')
NN = 20:20:400;
errD = zeros(size(NN)); err_un = zeros(size(NN)); errD2 = zeros(size(NN));
j = 1;
for N = NN                       	% convergence test
  s = setupquad(s, N);             	% set up quadr
  [f,gexac] = LapSLPmat(s,s_ps);
  f = f*den_source;                 % Dirichlet bdry condition
  gexac=gexac*den_source;           % exact Neu data
  if 1 % direct approach
      E = -sum(den_source);           	% condition @ inf: u->E*\log|x|, E==dot(du/dn,s.w)
      At = LapDLPnself(s,ord);       	% hypersingular T op
      [~,Adt] = LapSLPmat(s,s);         % D^* operator
      A = eye(N)/2+Adt+ones(N,1)*s.w.';	% construct A=0.5+D^*+1
      g = A\(At*f+E);               	% solve for g = du/dn
      u=LapDLPmat(t,s)*f-LapSLPmat(t,s)*g; % eval soln u = D*u - S*(du/dn)
      
      % Redo using Sidi's staggered grid finite-part correction
      At2 = LapDLPnself(s);            	% Sidi's corrected T op
      g2 = A\(At2*f+E);             	% solve for tau = du/dn
      u2=LapDLPmat(t,s)*f-LapSLPmat(t,s)*g2; % eval soln u = D*u - S*(du/dn)
      errD2(j) = max(abs(u2-uexac));   	% err at test pts
  else % indirect approach (1/2+(D+1))*tau=f
      E = sum(den_source);
      f = f - E/2/pi*log(1./abs(s.x)); % no blowup at inf
      Ad = LapDLPmat(s,s);      % D operator
      A = eye(N)/2+Ad+s.w.';  % 0.5+D
      tau = A\f;
      u = LapDLPmat(t,s)*tau+s.w.'*tau + E/2/pi*log(1./abs(t.x)); % add back
      g = LapDLPnself(s,ord)*tau - E/2/pi*real(s.x.*conj(s.nx))./abs(s.x).^2;
  end
  errD(j) = max(abs(u-uexac));     	% err at test pts
  err_un(j) = max(abs(gexac-g));   	% err Neu data
  fprintf('N=%d, \tpotential err = %.6g, \tNeu data err = %.6g\n', ...
            N,errD(j),err_un(j)); j = j+1;
end
%subplot(2,2,4); semilogy(NN,errD,'*',NN,err_un,'d'); legend('test pts','Neu data'); title('error')

% ext Neu BIE: (-1/2+D)*u=S*(du/dn), uniquely solvable
fprintf('ext Neumann problem...\n')
errN = zeros(size(NN)); err_u = zeros(size(NN));
j = 1;
for N = NN                        	% convergence test
  s = setupquad(s, N);            	% set up quadr
  [fexac,g] = LapSLPmat(s,s_ps);
  g=g*den_source;                   % Neumann bdry condition
  fexac = fexac*den_source;         % exact Diri data
  if 0 % direct
      Ad = LapDLPmat(s,s);             	% D operator
      As = LapSLPself(s,ord);         	% S operator
      A = -eye(N)/2+Ad;              	% construct operator A=-0.5+D
      f = A\(As*g);                   % solve for tau=u
      u=LapDLPmat(t,s)*f-LapSLPmat(t,s)*g; % eval soln u = D*u - S*(du/dn)
  else % indirect (-1/2+D^*)*tau=g
      [~, Adt] = LapSLPmat(s,s);
      A = -eye(N)/2+Adt;
      tau = A\g;
      u = LapSLPmat(t,s)*tau;
      f = LapSLPself(s,ord)*tau;
  end
  errN(j) = max(abs(u-uexac));    	% err at test pts
  err_u(j) = max(abs(f-fexac));    	% err Diri data
  fprintf('N=%d, \terr = %.6g, \tDiri data err = %.6g\n',...
      N,errN(j),err_u(j)); j = j+1;
end

subplot(1,4,4); semilogy(NN,errD,'*-',NN,errN,'d-',NN,errD2,'v-')
xlabel('N'); ylabel('err')
yticks([1e-16,1e-8,1e0])
ylim([1e-16,1e0])
% title('Laplace exterior BVP','interpreter','latex')
% legend({'Dirichlet (32nd)','Neumann (32nd)','Dirichlet (spec)'},'interpreter','latex')
