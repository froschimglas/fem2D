function G = Ex_Magnetic()
% Load geometry
G = Geometry(0.25);
G.shape = 'rectangle';

%G.plotEdges;

G = G.refine;
G.plot;

% Load Finite Elements on Triangles
FE2d = linearFE;
FE1d = quadraticFE;

% Assemble system matrix
%A = AssembleDomainMatrix(G,FE,@bilinearform);
A = AssembleMixedDomainMatrix(G,FE1d,FE2d,@bilinearform);

% Assemble right hand side
b = AssembleDomainVector(G,FE2d,@linearform);

% Apply boundary conditions
[A,b] = ApplyBcDirichlet(G,A,b);

% Solve linear equation
u = A\b;

% Plot solution
plotSolution(G,u,'solution');

end

% Bilinear form a(basis_i,basis_j) = (grad basis_i,grad basis_j)
function a = bilinearform( basisi,basisj, ~)

% x = points.x;
% y = points.y;

% u = points.u;
% v = points.v;
       
mu11=1;
mu22=1;


a = basisi.dx .* basisj.dx;


end

% Linear form f(phi) = (f,phi)
function f = linearform( basis, ~)

% x = points.x;
% y = points.y;

% u = points.u;
% v = points.v;
            
gamma = 1;

f = gamma.*basis.eval;

end
