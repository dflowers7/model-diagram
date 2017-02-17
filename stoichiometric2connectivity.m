function C = stoichiometric2connectivity(S)
%C = STOICHIOMETRIC2CONNECTIVITY(S) Convert stoichiometric matrix to
%connectivity matrix
%   Converts stoichiometric matrix S to connectivity matrix C. S is an
%   nx-by-nr matrix, where nx is the number of nodes and nr is the number
%   of reactions. Element (xi,ri) indicates the stoichiometry of node xi in
%   reaction ri. C is an nx-by-nx matrix, where nx is the number of nodes.
%   Element x1i,x2i indicates the number of reactions in which node x1i is
%   a source (or reactant) and node x2i is a sink (or product).

isSreactant = S < 0;
isSproduct = S > 0;

[nx,nr] = size(S);

% Form for including reactions as nodes
C = 1*[false(nx,nx) isSreactant;
    isSproduct' false(nr,nr)];

%Multiply by 1 to convert to double
%C = (1*isSreactant)*(1*isSproduct)';

% [nx,nr] = size(S);
% 
% C = zeros(nx);
% for ri = 1:nr
%     reactants = find(S(:,ri) < 0);
%     products = find(S(:,ri) > 0);
%     nreactants = length(reactants);
%     nproducts = length(products);
%     if isempty(reactants) || isempty(products)
%         continue
%     end
%     reactants = repmat(reactants(:),1,nproducts);
%     products = repmat(products(:)',nreactants,1);
%     Cind = sub2ind([nx nx],reactants(:),products(:));
%     C(Cind) = C(Cind) + 1;
% end




end

