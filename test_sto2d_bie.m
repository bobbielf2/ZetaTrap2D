% TEST_STO_BIE
% Test zeta-corrected trapezoidal rule for logarithmic singularity 
% on smooth geometries by solving the Stokes exterior Dirichlet problem 
% with smooth (artificial) boundary condition 
% and using mixed potential BIE formulation

mu = 1;     % viscosity
ord = 16; 	% pick desired convergence order of singular quaduature

% set up source geometry (starfish domain)
a = .3; w = 7; R = @(t) (1 + a*cos(w*t)); s.Z = @(t) R(t).*exp(1i*t);
Rp = @(t) -a*w*sin(w*t); s.Zp = @(t) exp(1i*t).*(Rp(t)+1i*R(t));

% generate (random) reference solution at test points
ns = 10;                                % num of point-forces
pf = randn(ns*2,1);                     % random point-forces
s_pf.x = .4*exp(2i*pi*rand(1,ns)).';	% random location of the pt-forces
s_pf.w = ones(ns,1);                    % dummy wei
t.x = 1.5*exp(2i*pi*(1:20)/20).';      	% test points
uu = StoSLP(t,s_pf,mu)*pf;
uexact = uu(1:end/2) + 1i*uu(end/2+1:end);	% ref soln at test points

% (optional) plot reference soln streamlines, point forces, & test points
nx = 150; gx = ((1:nx)/nx*2-1)*1.6;             % set up plotting grid
[xx, yy] = meshgrid(gx); zz = (xx+1i*yy); fhom = nan*(1+1i)*zz;
N = 600; xx = s.Z((0:N-1)'*(2*pi/N));           % sample nodes on curve
ii = inpolygon(real(zz),imag(zz),real(xx),imag(xx)); 
p.x = zz(~ii);                                  % targ pts outside
ff = StoSLP(p,s_pf,mu)*pf;
fhom(~ii) = ff(1:end/2) + 1i*ff(end/2+1:end);   % exact soln (complex form)
figure; streamslice(gx,gx,real(fhom),imag(fhom),'noarrows');
hold on; plot(xx,'k'); plot(s_pf.x,'.','MarkerSize',10,'LineWidth',10)
quiver(real(s_pf.x),imag(s_pf.x),pf(1:ns),pf(ns+1:end),'k')
plot(t.x,'dk','MarkerSize',8); axis equal off

% BVP convergence test
NN = 20:20:700;
err = zeros(size(NN));
j = 1;
for N = NN
    s = setupquad(s, N);                        % set up quadr
    f = StoSLP(s,s_pf,mu)*pf;	% Dirichlet bdry condition
    
    % solve 2nd kind BIE
    Ad = StoDLP(s,s,mu); As = StoSLP(s,s,mu,[],ord);
    A = eye(size(Ad))/2 + Ad + As;              % mixed potential formulation
    tau = A\f;                                  % solve for density tau
    
    % evaluate velocity
    uu = StoSLP(t,s,mu,tau)+StoDLP(t,s,mu,tau);
    u = uu(1:end/2) + 1i*uu(end/2+1:end);

    err(j) = max(abs(u-uexact));
    fprintf('N=%d, \terr = %.6g\n',N,err(j))
    j = j+1;
end

figure
semilogy(NN,err,'o'); ylabel('err'); xlabel('N')
