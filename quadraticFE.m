classdef quadraticFE
% linearFE class quadraticFE
%   The quadraticFE class realises quadratic finite elements over 1D
%   intervalls
%
% (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015    
    
   properties
       b1 = @(x,y) 2*(x-1).*(x-0.5);
       b2 = @(x,y) 2*x.*(x-0.5);
       b3 = @(x,y) 4*x-4*x.*x;
       
       db1= @(x,y) 4*x-3;
       db2= @(x,y) 4*x-1;
       db3= @(x,y) 4-8*x;
       
       % TODO: stimmt noch nicht!!!
       s1 =@(s) s;
       s2 =@(s) 1-s;
       
       ds1=@(s)  ones(length(s));
       ds2=@(s) -ones(length(s));
       
       nBasis = 3;
       m = 3;
   end
   methods
       function obj = quadraticFE(quadrule)
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
               grad = [grad;zeros(1,size(x,2))];
               grad = Z'\grad;
               phi.dx = grad(1,:);
               phi.dy = grad(2,:);
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