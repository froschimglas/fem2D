classdef linearFE
% linearFE class linearFE
%   The linearFE class realises linear finite elements over triangles
%   TODO: also important quadratic, cubic, quadrilaterals...
%
% (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015    
    
   properties
       b1 = @(x,y) 1-x-y;
       b2 = @(x,y) x;
       b3 = @(x,y) y;
       
       db1= @(x,y) [-1;-1];
       db2= @(x,y) [ 1; 0];
       db3= @(x,y) [ 0; 1];
       
       s1 =@(s) s;
       s2 =@(s) 1-s;
       
       ds1=@(s)  ones(length(s));
       ds2=@(s) -ones(length(s));
       
       nBasis = 3;
       m = 3;
   end
   methods
       function obj = linearFE(quadrule)
           % constructor
           if nargin>1
               obj.m = quadrule;
           end
       end
       function phi = basisfunction(obj,i,x,y)
           % returns basis function handle to local basis
           switch i
               case 1
                   phi = obj.b1(x,y);
               case 2
                   phi = obj.b2(x,y);
               case 3
                   phi = obj.b3(x,y);
           end
       end
       function phi = derivBasisfunction(obj,i,x,y,deriv)
           % returns handle to derivative of basis function to local basis

           if deriv >1
               error('only first derivative implemented yet');
           end
           
           switch i
               case 1
                   phi = obj.db1(x,y);
               case 2
                   phi = obj.db2(x,y);
               case 3
                   phi = obj.db3(x,y);
           end
       end
       
       function phi = basis(obj,i,x,y,Z,derivs)
           % basis function evaluation
           phi = struct('eval',[],'dx',[],'dy',[]);           
           
           if nargin < 6
               derivs = 1;
           end
           
           phi.eval = obj.basisfunction(i,x,y);
           
           if derivs >0
               grad = obj.derivBasisfunction(i,x,y,1);                                            
               grad = Z'\grad;
               phi.dx = grad(1)*ones(size(x));
               phi.dy = grad(2)*ones(size(x));
           end
       end       
       function phi = boundarybasis(obj,i,s)
           phi = struct('eval',[],'ds',[]);
           switch i
               case 1
                   phi.eval = obj.s1(s);
                   phi.ds = obj.ds1(s);
               case 2
                   phi.eval = obj.s2(s);
                   phi.ds = obj.ds2(s);
           end
       end
   end
end