function varargout = disperse(x)
if ndims(x)==2 && size(x,2)==1
    x = x';
end
if isnumeric(x) || ischar(x) || islogical(x) || isstruct(x)
    dims = 1:ndims(x)-1;
    varargout = num2cell(x,dims);
elseif iscell(x)
    if size(x,1) == 1
        varargout = x;
    else
        dims = 1:ndims(x)-1;
        varargout = num2cell(x,dims);
    end
else
    error('unknown data type');
end