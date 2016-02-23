function G= Ex_Poisson2D_homDirichlet(geo)
% Example of a Poisson equation and solution with linear finite elements
% INPUT: geo optional variable - if one already has a mesh
% OUTPUT: G geometry / mesh - especially useful if refinement ist used
%
% (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015

close all;

if nargin==0   
    % Load geometry
    G = Geometry(1);
    G.shape = 'square';
else   
    % Refine triangulation uniformly
    % G = geo.refine;
    G = geo;
end

G.plot

% Load Finite Elements on Triangles
FE = linearFE;
FE.m = 3; % default quadrature rule = 2 for linear finite elements

% Assemble right hand side
b = AssembleDomainVector(G,FE,@linearform);

% Assemble system matrix
A = AssembleDomainMatrix(G,FE,@bilinearform);

% Apply boundary conditions
[A,b] = ApplyBcDirichlet(G,A,b);

% Solve linear equation
u = A\b;

% Plot solution
plotSolution(G,u,'solution');

% Plot L-infty error
uexact=@(x,y) sin(2*pi * x).*sin(2* pi * y)/(8*pi^2);
plotError(G,u,uexact);

% Compute L2-error
L2err = exactError('L2',G,u,uexact,FE,A);
fprintf('Exact L2 error: %e\nDof: %d\n',L2err,G.nmbNodes);

% Estimate L2-error
f =@(x,y) sin(2*pi * x).*sin(2*pi * y);

[errEst,elemErr] = residualErrorEstimator(G,FE,f,u);
fprintf('Estimated H1 error: %e\nDof: %d\n',errEst,G.nmbNodes);


% refine elements which are the highest contributors
 markedElements = markElements(elemErr,'max');

 
 G = G.refine(markedElements);
 
end


% Linear form l(phi) = (f,phi)
function f = linearform( basis, points)

x = points.x;
y = points.y;

% u = points.u;
% v = points.v;
            
f = sin(2*pi * x).*sin(2*pi * y).*basis.eval;
%f = 1*basis.eval;
end

% Bilinear form a(basis_i,basis_j) = (grad basis_i,grad basis_j)
function a = bilinearform( basisi,basisj, ~)

% x = points.x;
% y = points.y;

% u = points.u;
% v = points.v;
       

a = basisi.dx .* basisj.dx + basisi.dy .* basisj.dy; % + basisj.eval.*basisi.eval
end


