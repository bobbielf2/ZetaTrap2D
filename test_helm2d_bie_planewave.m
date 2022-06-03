% TEST_HELM2D_BIE_PLANEWAVE
% Test zeta-corrected trapezoidal rule for Helmholtz layer potentials
% on smooth geometries by solving the problem of a planewave scattering
% from a sound-soft (Dirichlet) and sound-hard (Neumann) obstacle,
% via combined field integral equation formulation:
%      	(Dirichlet ansatz) u = (D - 1i*eta*S)*tau
%     	(Neumann ansatz)   u = (S + 1i*eta*D)*tau {altern. u=(D*S-1i*eta*S)*tau}

ka = 20+0i;     % wavenumber
ord = 32;       % pick desired convergence order of singular quad

% set up source geometry (starfish domain)
a = .2; w = 7; R = @(t) (0.7 + a*cos(w*t)); s.Z = @(t) R(t).*exp(1i*t);
Rp = @(t) -a*w*sin(w*t); s.Zp = @(t) exp(1i*t).*(Rp(t)+1i*R(t));

% set up mesh for plots
nx = 150; gx = ((1:nx)/nx*2-1)*1.6;         % set up plotting grid
[xx, yy] = meshgrid(gx); zz = (xx+1i*yy); fhom = nan*(1+1i)*zz;
N = 600; xx = s.Z((0:N-1)'*(2*pi/N));       % sample nodes on curve
ii = inpolygon(real(zz),imag(zz),real(xx),imag(xx)); p.x = zz(~ii); % targ pts outside

% test points for convergence
t.x = 1.5*exp(2i*pi*(1:20)/20).';

% BVPs convergence test
fprintf('Dirichlet problem...\n')
NN = [1000,60:20:940];
errD = zeros(size(NN(2:end)));
ieta = 1i*real(ka);
j = 0;
for N = NN
  s = setupquad(s, N);                            	% set up quadr
  f = -exp(1i*ka*real(exp(1i*pi/4)*s.x));           % Dirichlet bdry condition
  Ad = HelmDLP(ka,s,s,[],[],ord);
  As = HelmSLP(ka,s,s,[],[],ord);
  A = eye(N)/2 + Ad - ieta * As;                	% construct CFIE operator A=0.5+D-i*eta*S
  tau = A\f;                                        % solve CFIE for tau
  u = (HelmDLP(ka,t,s)-ieta*HelmSLP(ka,t,s))*tau;	% eval soln u = (D-i*eta*S)*tau
  if j == 0, uref = u; % ref. soln
  else
      errD(j) = max(abs(u-uref));                     % err
      fprintf('N=%d, \terr = %.6g\n',N,errD(j));
  end
  if N == 300
      fhom(~ii) = (HelmDLP(ka,p,s)-ieta*HelmSLP(ka,p,s))*tau;
      fhom(~ii) = fhom(~ii) + exp(1i*ka*real(exp(1i*pi/4)*p.x));
      figure; imagesc(gx,gx,real(fhom));
      axis('xy','equal','tight')
      colorbar; 
      caxis([-1,1]*max(real(fhom(~ii))))
      hold on; fill(real(xx),imag(xx),'w');
  end
  j = j+1;
end


fprintf('Neumann problem...\n')
reg = 0; % regularize? if 1 build regularized BIE; if 0 build hypersingular BIE.
errN = zeros(size(NN(2:end)));
j = 0;
for N = NN
  s = setupquad(s, N);                              % set up quadr
  g = -exp(1i*ka*real(exp(1i*pi/4)*s.x)).*(1i*ka*real(exp(1i*pi/4)*s.nx)); % Neumann bdry condition
  if reg % construct regularized CFIE operator A=(0.5*i*eta-0.25)+D'^2-i*eta*D'
      Ad = HelmSLPn(ka,s,s,[],[],ord);
      As = HelmSLP(ka,s,s,[],[],ord);
      A = (0.5*ieta-0.25)*eye(N)+Ad*(Ad-ieta*eye(N));  	
  else % construct hypersingular CFIE operator A=-0.5+D'+i*eta*T
      Ad = HelmSLPn(ka,s,s,[],[],ord);
      At = HelmDLPnself(ka,s,ord);
      A = -0.5*eye(N)+Ad+ieta*At;  	
  end
  tau = A\g;   % solve CFIE for tau
  if reg
      u = HelmDLP(ka,t,s)*(As*tau)-ieta*HelmSLP(ka,t,s)*tau; % eval soln u = (D*S-i*eta*S)*tau
  else
      u = (HelmSLP(ka,t,s)+ieta*HelmDLP(ka,t,s))*tau;        % eval soln u = (S+i*eta*D)*tau
  end
  if j == 0, uref = u; % ref. soln
  else
      errN(j) = max(abs(u-uref));                   % err
      fprintf('N=%d, \terr = %.6g\n',N,errN(j));
  end
  if N == 300
      if reg
          fhom(~ii) = HelmDLP(ka,p,s)*(As*tau)-ieta*HelmSLP(ka,p,s)*tau;
      else
          fhom(~ii) = (HelmSLP(ka,p,s)+ieta*HelmDLP(ka,p,s))*tau;
      end
      fhom(~ii) = fhom(~ii) + exp(1i*ka*real(exp(1i*pi/4)*p.x));
      figure; imagesc(gx,gx,real(fhom));
      axis('xy','equal','tight')
      colorbar; 
      caxis([-1,1]*max(real(fhom(~ii))))
      hold on; fill(real(xx),imag(xx),'w');
  end
  j = j+1;
end

NN = NN(2:end);
figure; semilogy(NN,errD,'*',NN,errN,'d')
xlabel('N'); ylabel('err')
title(['Helmholtz BVP convergence, $\kappa=',num2str(ka),'$'],'interpreter','latex')
legend({'Dirichlet','Neumann'},'interpreter','latex')