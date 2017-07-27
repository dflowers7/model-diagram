function [A, B, varargout] = sxrepmat(A, B, varargin)
%SXREPMAT Resizes dimensions of arrays to match one another by singleton
%expansion.
%   [A, B, ...] = SXREPMAT( A, B, ... ) replicates A, B, and all other
%   inputted arrays along their singleton dimensions so that they match in
%   size.

assert(nargin >= nargout, 'Number of inputs cannot be less than the number of outputs.')

arrays = cell(nargin,1);
arrays{1} = A;
arrays{2} = B;
if ~isempty(varargin)
    arrays(3:end) = varargin;
end

% Determine which dimensions are singleton
ndim = cellfun(@ndims,arrays);
sizes = cell(max(ndim),1);
[sizes{:}] = cellfun(@size,arrays);
sizes = [sizes{:}]; % sizes: nargin-by-ndims array of dim sizes
issingleton = sizes == 1;

% Get new sizes
endsizes = sort(sizes,1);
endsizes = repmat(endsizes(end,:),nargin,1);

% Ensure non-singleton dimensions are the same in each array
assert(all(all(issingleton | sizes == endsizes)),'Non-singleton dimensions do not match.')

% Determine rep values
reps = zeros(size(sizes));
reps(issingleton) = endsizes(issingleton);
reps(~issingleton) = ones(sum(sum(~issingleton)),1);

reps = mat2cell(reps,ones(nargin,1),max(ndim));
arrays = cellfun(@repmat,arrays,reps,'UniformOutput',false);

A = arrays{1};
B = arrays{2};
if nargout > 2
    varargout = arrays(3:end);
end

end

