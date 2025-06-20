function pp = multiply_poles(p, m)

% MULTIPLY_POLES - Duplicates the elements of 'p' by the elements of 'm'.  
%
% Usage: 
%     pp = multiply_poles(p,m)
%
% Input parameters:
%     p : row vector that contains a pole only once
%     m : multiplicities related to the pole vector 'p'
%
% Output parameters:
%     pp : vector of the poles that contains the ith element of 'p' 
%          only 'm(i)' times
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

if length(p) ~= length(m)
    error('Bad poles, length of p and m must be equal!');
end

n = size(p,2);
pp = zeros(1,sum(m));
innen = 1;
for i = 1:n
    pp(1,innen:innen+m(i)-1) = p(i)*ones(1,m(i));
    innen = innen + m(i);
end
