% TEST_STO2D_ZETA
% Test high-order discretization of 2D Stokes potential self-interactions
% via zeta quadrature on smooth closed curves. 
% The involved singularities are as follows:
%   - weakly singular:  SLP velocity
%   - singular:         SLP pressure
%   - hypersingular:    DLP pressure & traction
%   - regular:          SLP traction & DLP velocity

mu = 0.7;  	% viscosity
ord = 16; 	% pick desired convergence order of singular quaduature
rng(3);     % set rand seed for replicability

% set up source geometry (starfish domain)
rph = randn()*pi;   % random phase angle
a = .3; w = 5; R = @(t) (1 + a*cos(w*t)); s.Z = @(t) R(t-rph).*exp(1i*t);
Rp = @(t) -a*w*sin(w*t); s.Zp = @(t) exp(1i*t).*(Rp(t-rph)+1i*R(t-rph));
Rpp = @(t) -a*w^2*cos(w*t); s.Zpp = @(t) exp(1i*t).*(2i*Rp(t-rph)-R(t-rph)+Rpp(t-rph));

% plot curve
s = setupquad(s, 200);  % set up quadr
subplot(1,3,1)
plot(s.x,'.','markersize',4), hold on
h = plot(real(s.x(1)),imag(s.x(1)),'o','linewidth',1.5);
hold off, axis equal
legend(h,'target')

% convergence test at the target location s.x(1), both SLP & DLP
NN = [];
err_us = []; err_ps = []; err_ts = [];
err_ud = []; err_pd = []; err_td = [];
j = 0;
for N = [1000,40:20:500]
    s = setupquad(s, N);                        % set up quadr
    tau = [exp(cos(3*s.t));exp(sin(3*s.t).^2)]; % smooth density
    
    % potential matrices via zeta quadrature
    [As,Ps,Ts] = StoSLP(s,s,mu,[],ord);
    [Ad,Pd,Td] = StoDLP(s,s,mu,[],ord);
    
    % evaluate SLP
    us = As*tau;
    ps = Ps*tau;
    ts = Ts*tau;
    
    % evaluate DLP
    ud = Ad*tau;
    pd = Pd*tau;
    td = Td*tau;
    
    % error at s.x(1)
    if j == 0   % reference values
        us_ref = us(1); ps_ref = ps(1); ts_ref = ts(1);
        ud_ref = ud(1); pd_ref = pd(1); td_ref = td(1);
    else
        NN(j) = N;
        err_us(j) = max(abs(us(1)-us_ref));
        err_ps(j) = max(abs(ps(1)-ps_ref));
        err_ts(j) = max(abs(ts(1)-ts_ref));
        fprintf('N=%d, \terr(us) = %.6g \terr(ps) = %.6g \terr(ts) = %.6g\n',N,err_us(j),err_ps(j),err_ts(j))
        err_ud(j) = max(abs(ud(1)-ud_ref));
        err_pd(j) = max(abs(pd(1)-pd_ref));
        err_td(j) = max(abs(td(1)-td_ref));
        fprintf(' \terr(ud) = %.6g \terr(pd) = %.6g \terr(td) = %.6g\n',err_ud(j),err_pd(j),err_td(j))
    end
    j = j+1;
end

subplot(1,3,2)
semilogy(NN,err_us,'o',NN,err_ps,'+',NN,err_ts,'v','linewidth',1.5)
title('Stokes SLP'); ylabel('err'); xlabel('N'); axis([0,500,1e-16,1e-1])
legend({'velocity (weak singular)','pressure (singular)','traction (regular)'})

subplot(1,3,3)
semilogy(NN,err_ud,'o',NN,err_pd,'+',NN,err_td,'v','linewidth',1.5); 
title('Stokes DLP'); ylabel('err'); xlabel('N'); axis([0,500,1e-15,1e0])
legend({'velocity (regular)','pressure (hypersingular)','traction (hypersingular)'})

%% separate convergence tests
lptype = 'd'; % s = SLP, d = DLP
NN = [];
err_u = []; err_p = []; err_t = [];
j = 0;
for N = [800,20:20:400]
    s = setupquad(s, N);                        % set up quadr
    tau = [exp(cos(3*s.t));exp(sin(3*s.t).^2)]; % smooth density
    
    % potential matrices via zeta quadrature
    if lptype == 's'
        [A,P,T] = StoSLP(s,s,mu,[],ord);
    else
        [A,P,T] = StoDLP(s,s,mu,[],ord);
    end

    % evaluate
    u = A*tau;
    p = P*tau;
    t = T*tau;
    
    % error at s.x(1)
    if j == 0   % reference values
        u_ref = u(1); p_ref = p(1); t_ref = t(1);
    else
        NN(j) = N;
        err_u(j) = max(abs(u(1)-u_ref));
        err_p(j) = max(abs(p(1)-p_ref));
        err_t(j) = max(abs(t(1)-t_ref));
        fprintf('N=%d, \terr(u) = %.6g \terr(p) = %.6g \terr(t) = %.6g\n',N,err_u(j),err_p(j),err_t(j))
    end
    j = j+1;
end

subplot(1,2,2)
semilogy(NN,err_u,'o',NN,err_p,'+',NN,err_t,'v'); 
ylabel('err'); xlabel('N')
if lptype == 's'
    title('SLP')
    legend({'velocity (weak singular)','pressure (singular)','traction (regular)'})
else
    title('DLP')
    legend({'velocity (regular)','pressure (hypersingular)','traction (hypersingular)'})
end