function err = exactError(errorType,geo,uh,u,FE,A)
% Computes the error e=u-uh of given solution u and numerical solution uh
% in currently only the L2-norm
% INPUT: errorType 'L2' - other norms are not yet implemented
%        geo geometry instance of Geometry class
%        uh numerical solution, i.e. vector of solution coefficients
%        FE the finite element rule used to compute uh
%        A optional variable, the precomputed stiffness/system matrix
%
% (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015

switch errorType
    case 'L2'
        %compute || u-uh||_L2(Omega)
        if nargin < 4
            % assemble stiffness matrix again
        else
            Gauss = Quad2DTriangle(FE.m);
            
            c=constantPart(geo,u,Gauss);                
            b = mixedPart(geo,FE,u,Gauss);          
            
            err = sqrt(c+uh'*A*uh - 2*b'*uh); 
        end
        
        
    otherwise
        error('This type of error is not yet implemented');
end


end


function c = constantPart(geo,u,Gauss)
c=0;
for i=1:geo.nmbElements
    K=geo.getLocalElement(i);
    ci = L2NormDomain(K,u,Gauss);
    
    c = c + ci*ci;
    
    %c= c+integrate(Gauss,f,K);
end
end

function b = mixedPart(geo,FE,u,Gauss)
b = zeros( geo.nmbNodes,1 );

for i=1:geo.nmbElements
    K=geo.getLocalElement(i);
    [P,Z] = Gauss.trafo(K);
    
    globalBasis = geo.globalNodeIndices(i);
     
    for j=1:3
        f =@(x,y) u(x,y);
        ftrafo = @(x,y) trafof(x,y,f,P,Z);
        
         
        fun =@(x,y) ftrafo(x,y).*FE.basisfunction(j,x,y);
        ymax =@(x) 1-x;
        
        globbasis = globalBasis(j);
        
        b(globbasis) = b(globbasis)+ quad2d(fun,0,1,0,ymax);        
    end
end


end

