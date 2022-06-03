function [u,p,T] = StoDLP(t,s,mu,dens,ord)
% STODLP   Evaluate 2D Stokes double-layer velocity, pressure, and traction.
%
% [A,P,T] = StoDLP(t,s,mu,[],ord) returns dense matrices taking double-layer
%  density values on the nodes of a source curve to velocity, pressure, and
%  traction on the nodes of a target curve. Native quadrature is used; when
%  target=source, A is the self-interaction Nystrom matrix using the smooth
%  diagonal-limit formula for DLP, P and T are corrected to high-order
%  (specified by the input 'ord') using the "zeta quadrature" described in
%  [WM21] and [WM22]. 
%
% [u,p,T] = StoDLP(t,s,mu,dens,ord) evaluates the double-layer density dens,
%  returning flow velocity u, pressure p, and target-normal traction T.
%
% The normalization is as in Sec 2.3 of [HW], and [MBGV].
%
% References:
%  [HW]    "Boundary Integral Equations", G. C. Hsiao and W. L. Wendland
%          (Springer, 2008).
%
%  [MBGV]  "A fast algorithm for simulating multiphase flows through periodic
%          geometries of arbitrary shape," G. Marple, A. H. Barnett,
%          A. Gillman, and S. Veerapaneni, SIAM J. Sci. Comput., 38(5), pp.
%          B740-B772, 2016.
%
%  [WM21]  "Zeta correction: a new approach to constructing corrected
%          trapezoidal quadrature rules for singular integral operators",
%          B. Wu, P.G. Martinsson, Adv. Comput. Math., 47(3), 1–21, 2021.
% 
%  [WM22]  "A Unified Trapezoidal Quadrature Method for Singular and
%          Hypersingular Boundary Integral Operators on Curved Surfaces",
%          B. Wu and P.G. Martinsson, in preparation.
%
% Inputs: (see setupquad for source & target struct definitions)
%  s = source segment struct with s.x nodes, s.w weights on [0,2pi),
%      s.sp speed function |Z'| at the nodes, and s.tang tangent angles.
%  t = target segment struct with t.x nodes, and t.nx normals if traction needed
%  mu = viscosity
%  ord = order of correction for pressure & traction self-evaluation
%
% Outputs: (matrix case)
%  A = 2M-by-2N matrix taking density (force vector) to velocity on the
%      target curve. As always for Stokes, ordering is nodes fast,
%      components (1,2) slow, so that A has 4 large blocks A_11, A_12, etc.
%  P = M-by-2N matrix taking density to pressure (scalar) on target nodes.
%  T = 2M-by-2N matrix taking density to normal traction on target nodes.
%
% To test use STOINTDIRBVP
%
% See also: SETUPQUAD

% Barnett 6/13/16, 6/27/16
% Bowei Wu 6/30/21, added self-eval correction for pressure & traction

if numel(mu)~=1, error('mu must be a scalar'); end
if nargin < 5 || isempty(ord)
  ord = 16; % default: 16-th order zeta correction
end

if nargout==1
  u = StoDLPmat(t,s,mu,ord);
  if nargin>3 && ~isempty(dens)
    u = u * dens;
  end
elseif nargout==2
  [u p] = StoDLPmat(t,s,mu,ord);
  if nargin>3 && ~isempty(dens)
    u = u * dens;
    p = p * dens;
  end
else
  [u p T] = StoDLPmat(t,s,mu,ord);
  if nargin>3 && ~isempty(dens)
    u = u * dens;
    p = p * dens;
    T = T * dens;
  end
end
%%%%%%

function [A,P,T] = StoDLPmat(t,s,mu,ord)
% Returns native quadrature matrices, or self-evaluation matrices.
% some repmat->bsxfun 6/28/16
% added self-eval correction for pressure & traction, Bowei Wu 2021

self = sameseg(t,s);
N = numel(s.x); M = numel(t.x);
r = bsxfun(@minus,t.x,s.x.');                 % C-# displacements mat
irr = 1./(conj(r).*r);    % 1/r^2, used in all cases below
d1 = real(r); d2 = imag(r);
%rdotny = d1.*repmat(real(s.nx)', [M 1]) + d2.*repmat(imag(s.nx)', [M 1]);
rdotny = bsxfun(@times, d1, real(s.nx)') + bsxfun(@times, d2, imag(s.nx)');
rdotnir4 = rdotny.*(irr.*irr); if nargout==1, clear rdotny; end
A12 = (1/pi)*d1.*d2.*rdotnir4;  % off diag vel block
A = [(1/pi)*d1.^2.*rdotnir4,   A12;                     % Ladyzhenskaya
     A12,                      (1/pi)*d2.^2.*rdotnir4];
if sameseg(t,s)         % self-int DLP velocity diagonal correction
  c = -s.cur/2/pi;           % diagonal limit of Laplace DLP
  tx = 1i*s.nx; t1=real(tx); t2=imag(tx);     % tangent vectors on the curve
  A(sub2ind(size(A),1:N,1:N)) = c.*t1.^2;     % overwrite diags of 4 blocks
  A(sub2ind(size(A),1+N:2*N,1:N)) = c.*t1.*t2;
  A(sub2ind(size(A),1:N,1+N:2*N)) = c.*t1.*t2;
  A(sub2ind(size(A),1+N:2*N,1+N:2*N)) = c.*t2.^2;
end
A = bsxfun(@times, A, [s.w(:)' s.w(:)']);            % quadr wei

% Compute correction weights for pressure & traction self-interaction
if self
    if nargout>2
      [ZP,ZT] = StoDLPZetaSparse(s,mu,ord);
    else
      ZP = StoDLPZetaSparse(s,mu,ord);
    end
end
if nargout>1           % pressure of DLP
%  P = [(mu/pi)*(-repmat(real(s.nx)', [M 1]).*irr + 2*rdotnir4.*d1), ...
%       (mu/pi)*(-repmat(imag(s.nx)', [M 1]).*irr + 2*rdotnir4.*d2)  ];
  P = [bsxfun(@times, real(-s.nx)',irr) + 2*rdotnir4.*d1, ...
       bsxfun(@times, imag(-s.nx)',irr) + 2*rdotnir4.*d2  ];
  P = bsxfun(@times, P, (mu/pi)*[s.w(:)' s.w(:)']);            % quadr wei
  if self
    P(logical([speye(N),speye(N)])) = 0;
    P = P + ZP; % correction for self interact
  end
end
if nargout>2           % traction, my derivation of formula
  nx1 = repmat(real(t.nx), [1 N]); nx2 = repmat(imag(t.nx), [1 N]);
  rdotnx = d1.*nx1 + d2.*nx2;
  ny1 = repmat(real(s.nx)', [M 1]); ny2 = repmat(imag(s.nx)', [M 1]);
  dx = rdotnx.*irr; dy = rdotny.*irr; dxdy = dx.*dy;
  R12 = d1.*d2.*irr; R = [d1.^2.*irr, R12; R12, d2.^2.*irr];
  nydotnx = nx1.*ny1 + nx2.*ny2;
  T = R.*kron(ones(2), nydotnx.*irr - 8*dxdy) + kron(eye(2), dxdy);
  T = T + [nx1.*ny1.*irr, nx1.*ny2.*irr; nx2.*ny1.*irr, nx2.*ny2.*irr] + ...
      kron(ones(2),dx.*irr) .* [ny1.*d1, ny1.*d2; ny2.*d1, ny2.*d2] + ...
      kron(ones(2),dy.*irr) .* [d1.*nx1, d1.*nx2; d2.*nx1, d2.*nx2];
  T = bsxfun(@times, T, (mu/pi)*[s.w(:)' s.w(:)']);    % prefac & quadr wei
  if self
      T(logical(kron(ones(2),speye(N)))) = 0;
      T = T + ZT; % correction for self interact
  end
end

function a = sameseg(t,s)
% SAMEREG   Return true if two segments have numerically the same nodes

% Barnett 6/12/16
a = (numel(s.x)==numel(t.x)) && max(abs(s.x-t.x))<1e-14;