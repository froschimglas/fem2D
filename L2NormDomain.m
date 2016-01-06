function c = L2NormDomain(elem,f,Gauss)
% Computes the L2-Norm of a function f over an element using MATLABS adaptive quadrature rules
% INPUT: elem element (triangle) instance of the ELEMENT class
%        f function handle to function
%        Gauss Quadrature rule - to obtain the transformation from elem to
%        the reference triangle
% OUTPUT: c = sqrt ( int_elem |f|^2 dx )
%
% (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015



    [P,Z] = Gauss.trafo(elem);
     
    fsquared =@(x,y) f(x,y).*f(x,y);
    
    ftrafo = @(x,y) trafof(x,y,fsquared,P,Z) ;
    ymax =@(x) 1-x;
    
    c = quad2d(ftrafo,0,1,0,ymax)* Gauss.detTrafo(elem)*2;
    c = sqrt(c);
    
    %c= integrate(Gauss,f,elem);

end