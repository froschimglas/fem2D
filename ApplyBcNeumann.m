function b = ApplyBcNeumann(geo,femrule,b,g)
% applys Neumann boundary conditions by
% elimination entries in system matrix and load vector.
% INPUT:    b: load vector
% OUTPUT:   b: modified load vector
%
% (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015

% TODO: Apply Neumann Boundary Conditions only on elements given by user input 

Gauss = Quad1D(femrule.m);

e=geo.getBCEdges;
for bcEdge =1:size(e,1)
    elimEdge = e(bcEdge,:);
    A = geo.nodes( elimEdge(1),:);
    B = geo.nodes( elimEdge(2),:);
    
    detJFK = norm(B-A);
    
    fun=@(s) g( A(1)+(B(1)-A(1))*s, A(2)+(B(2)-A(2))*s  );
    
    bcBasis1 = elimEdge(1);
    bcBasis2 = elimEdge(2);
    
    fun1=@(s) fun(s).*femrule.basis(1,s);   % belongs to node A
    fun2=@(s) fun(s).*femrule.basis(2,s);   % belongs to node B
    
    val1 = integrate(Gauss,fun1,0,1 );
    val2 = integrate(Gauss,fun2,0,1 );
    
    b(bcBasis1) = b(bcBasis1) + val1*detJFK;
    b(bcBasis2) = b(bcBasis2) + val2*detJFK;
    
end
end