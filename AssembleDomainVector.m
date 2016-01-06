function b = AssembleDomainVector(geo,femrule,linearform)
% Assembles a vector b = (f,phi_i)_i from a given linearform and finite
% element rule
% INPUT: geo geometry - instance of geometry class
%        femrule - finite element method
%        linearform - function handle to linear form of user
%
% (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015


ndof = geo.nmbNodes;
nelem = geo.nmbElements;
b = zeros(ndof,1);

% loop over all elements
for i=1:nelem
    
    % get local element properties
    K = geo.getLocalElement(i);
    
    % get element entries
    bK = elementEntries( K,femrule,linearform );
    
    % get global basis indices
    globalbasis = geo.globalNodeIndices(i);
    
    % put values into load vector
    b(globalbasis) = b(globalbasis) + bK;

end

end

function y = elementEntries(K,femrule,linearform)

% get quadrature points
Gauss = Quad2DTriangle(femrule.m);
W = Gauss.getWeights;
points = Gauss.getQuadPointInfo(K);

% pre evaluate element basis functions
[~,Z] = Gauss.trafo(K);
basisfunctions = cell(3,1);
for i=1:3
    basisfunctions{i} = femrule.basis(i,points.refx,points.refy,Z);
end

% test linearform with each basis function
y = [0;0;0];
detTrafo = Gauss.detTrafo(K);

for i=1:3
    phi = basisfunctions{i};
    Y = linearform( phi, points);
    
    % Gauss integration
    y(i) = W*Y'*detTrafo;
end

end