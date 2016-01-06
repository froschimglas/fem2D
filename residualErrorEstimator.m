function [err,elemErr] = residualErrorEstimator(geo,FE,f,sol)
% Residualbased error estimator for the error e=u-uh in the L2-Norm of 
% the Poisson equation with homogeneous Dirichlet boundary conditions. 
% INPUT: 
% geo geometry-instance of Geometry class
% FE finite element rule
% f function handle to right hand side of PDE
% sol numerical solution vector of Poisson equation which was computed with
% FE finite element rule
%
% TODO: other error estimators
% FIXME: does it work correctly? Check...
%
% (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015



elemErr = getElementError(geo,FE,f);
edgeErr = getElementEdgeError(geo,FE,sol);

elemErr = elemErr.^2 + edgeErr;

err = sum( elemErr );
err = sqrt(err);



end

function elemErr = getElementError(geo,FE,f)
elemErr = zeros( geo.nmbElements,1 );
Gauss = Quad2DTriangle(FE.m);

for i=1:geo.nmbElements
    K=geo.getLocalElement(i);
    
    elemPart = L2NormDomain(K,f,Gauss);
    hK = K.diameter;
    
    elemErr(i) = hK*elemPart;
end

end
function elemErr = getElementEdgeError(geo,FE,sol)
elemErr = zeros( geo.nmbElements,1 );
for elem =1:geo.nmbElements
    
    K=geo.getLocalElement(elem);
    nodeOrder = geo.globalNodeIndices(elem);
    u = sol( nodeOrder );
    
    
    gradK = u(1)* FE.derivBasisfunction(1,0,0,1) + ...
        u(2)* FE.derivBasisfunction(2,0,0,1) + ...
        u(3)* FE.derivBasisfunction(3,0,0,1);
    
    %loop over all three edges
    for edg=1:3
        %get nodes of this edge
        N = K.edgeNodes(edg);
        
        % check if it is a boundary edge
        isNoBoundaryEdge = ~geo.isBoundaryEdge( elem,N );
        
        
        if isNoBoundaryEdge
            
            %get neighbour of this edge
            SI = edgeAttachments(geo.getTRep, nodeOrder(N));
            neighElem = setdiff(SI{1},elem);
            
            
            
            %get normal on edge
            normal = K.normal(edg);
            
            %evaluate gradient of solution on this and the neighbour element
            nodeOrderN = geo.globalNodeIndices(neighElem);
            uN = sol( nodeOrderN );
            gradN = uN(1)* FE.derivBasisfunction(1,0,0,1) + ...
                uN(2)* FE.derivBasisfunction(2,0,0,1) + ...
                uN(3)* FE.derivBasisfunction(3,0,0,1);
            
            
            jump = normal'*gradK -normal'*gradN;
            hE = norm( K.edge(edg) );
        end
        
        elemErr( elem ) = elemErr(elem) + (0.5*jump*hE)^2;
    end
end
end
