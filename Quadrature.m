classdef Quadrature
% Quadrature class Quadrature: abstract Quadrature class
%   The Quadrature class is an abstract quadrature rule
%
% (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015
    
    properties(Access = public)
         % m -> determines number of quadrature points
         m;
     end
    methods(Access = public)    
        function obj = Quadrature(quadrule)
            % constructor
            % INPUT:             
            %  quadrule:      optional quadrature parameter. default is 3.
                           
            if nargin>0            
                % quadrule must be an integer between 1 and 3
                if quadrule<0 || quadrule>3
                    error('only quadrature rules with m = 1,...,3 supported');
                else
                    obj.m = quadrule;                
                end
            else
                obj.m = 3;
            end                       
        end               
        function val = integrate(obj,f,intval)
            W = obj.getWeights;
            Y = obj.getQuadValues(f,intval);
            val = W*Y'*obj.detTrafo(intval);
        end
        function points = getQuadPointInfo(obj,intval)
            XRef = obj.getReferenceQuadPoints;
            X = obj.getQuadPoints(intval);

            points = struct('x',X(1,:),'y',X(2,:),'refx',XRef(1,:),'refy',XRef(2,:) );
            
        end
    end
    methods(Abstract)
        W = getWeights(obj);
        X = getReferenceQuadPoints(obj);
        X = getQuadPoints(obj);
        Y = getQuadValues(obj,f,intval);
    end
    methods(Static)
        [P,Z] = trafo(intval);
        dt = detTrafo(intval);
    end
end

