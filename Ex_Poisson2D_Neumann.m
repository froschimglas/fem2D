function G= Ex_Poisson2D_Neumann(geo)
% Example of a Poisson equation and solution with linear finite elements
% INPUT: geo optional variable - if one already has a mesh
% OUTPUT: G geometry / mesh - especially useful if refinement ist used
%
% (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015

close all;

if nargin==0   
    % Load geometry
    G = Geometry(0.2);
    G.shape = 'rectangle';
else   
    % Refine triangulation uniformly
    % G = geo.refine;
    G = geo;
end

G.plot

% Load Finite Elements
FE = linearFE;

% Assemble right hand side
b = AssembleDomainVector(G,FE,@linearform);

% Assemble system matrix
A = AssembleDomainMatrix(G,FE,@bilinearform);

% Apply boundary conditions
%[A,b] = ApplyBcDirichlet(G,A,b);
% -> natural Neumann (  grad u * n = 0 )

% Solve linear equation
u = A\b;

% Plot solution
plotSolution(G,u,'solution');
 
end


% Linear form l(phi) = (f,phi)
function f = linearform( basis, points)

x = points.x;
y = points.y;
% u = points.u;
% v = points.v;
            
f = sin(pi/2 * x).*sin(pi/2 * y).*basis.eval;
%f = 1*basis.eval;
end

% Bilinear form a(phii,phij) = (phii,phij) + (grad phii,grad phij)
function a = bilinearform( basisi,basisj, ~)

% x = points.x;
% y = points.y;

% u = points.u;
% v = points.v;
       

a = basisi.eval.*basisj.eval + basisi.dx .* basisj.dx + basisi.dy .* basisj.dy;
end


