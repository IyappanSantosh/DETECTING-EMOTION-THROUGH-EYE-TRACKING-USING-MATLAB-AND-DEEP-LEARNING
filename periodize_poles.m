function pp=periodize_poles(p, m)

% PERIODIZE_POLES - Duplicates periodically the elements of 'p' 'm' times.  
%
% Usage: 
%     pp=periodize_poles(p,m)
%
% Input parameters:
%     p : row vector that contains the poles
%     m : integer factor of duplication
%
% Output parameters:
%     pp : vector of the poles that contains 'p' sequentially
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

pp=zeros(1,m*length(p));
for i=1:1:m
    pp((i-1)*length(p)+1:i*length(p))=p;
end