function [A,An] = LapSLPmat(t,s)
% LAPSLPMAT  2D Laplace single-layer matrix, interaction from curve to
% target. No self correction for A.
%
% Inputs
%   t:      target struct, t.x = target points
%   s:      source struct, output of setupquad.m
% Outputs
%   A:      potential matrix
%   An:     normal derivative

% Bowei Wu 1/13/21

d = bsxfun(@minus,t.x,s.x.'); 	% displacements mat
A = -log(abs(d));               % peri log
A = bsxfun(@times, A, (1/2/pi)*s.w(:).');   % prefactor & wei
if (numel(s.x)==numel(t.x)) && max(abs(s.x-t.x))<1e-14  % self eval
    A(diagind(A)) = 0; % zero out
end
if nargout > 1
    An = real(-t.nx./d);
    if (numel(s.x)==numel(t.x)) && max(abs(s.x-t.x))<1e-14  % self eval
        An(diagind(An)) = -.5*s.cur;	% diagonal limit
    end
    An = bsxfun(@times, An, (1/2/pi)*s.w(:)');  % include weights & prefac
end

function i = diagind(A)
% DIAGIND  Return indices of diagonal of square matrix
% 
% Example usage:   A = randn(3,3); A(diagind(A)) = 0;
N = size(A,1); i = sub2ind([N,N], 1:N, 1:N); i = i(:);