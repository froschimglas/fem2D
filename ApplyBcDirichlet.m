function [A,b] = ApplyBcDirichlet(geo,A,b)
% applys homogeneous Dirichlet boundary conditions by
% elimination entries in system matrix and load vector.
% INPUT:    A: system matrix (resulting from bilinear form)
%           b: load vector
% OUTPUT:   A: modified system matrix
%           b: modified load vector
%
% (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015

% TODO: Apply Dirichlet Boundary Conditions only on elements given by user input 

e=geo.getBCEdges;
for bcEdge =1:size(e,1)
    elimEdge = e(bcEdge,:);
    for bcNode = 1:2
        elimNode = elimEdge(bcNode);
        A(elimNode,:)=0;
        A(:,elimNode)=0;
        A(elimNode,elimNode)=1;
        b(elimNode) = 0;
    end
end
end