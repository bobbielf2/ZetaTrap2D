function [A,An] = LapDLPmat(t,s)
% LAPDLPMAT  2D Laplace double-layer matrix, interaction from curve to 
% target.
% 
% Inputs
%   t:      target struct, t.x = target points
%   s:      source struct, output of setupquad.m
% Outputs
%   A:      potential matrix
%   An:     normal derivative

% BW Jan 21

d = bsxfun(@minus,t.x,s.x.');          % displacements mat
A = real(s.nx.'./d);
if (numel(s.x)==numel(t.x)) && max(abs(s.x-t.x))<1e-14  % self eval
    A(diagind(A)) = -.5*s.cur;	% diagonal limit
end
A = bsxfun(@times, A, (1/2/pi)*s.w(:)');  % include weights & prefac
if nargout > 1
    r = abs(d);
    csry = s.nx'.*d;            % (cos phi + i sin phi).r
    csrx = conj(t.nx).*d;       % (cos th + i sin th).r
    An = -real(csry.*csrx)./((r.^2).^2);         % divide is faster than bxsfun here
    An = bsxfun(@times, An, (1/2/pi)*s.w(:)');    % prefac & quadr wei
end


function i = diagind(A)
% DIAGIND  Return indices of diagonal of square matrix
% 
% Example usage:   A = randn(3,3); A(diagind(A)) = 0;
N = size(A,1); i = sub2ind([N,N], 1:N, 1:N); i = i(:);