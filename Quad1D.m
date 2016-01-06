classdef Quad1D < Quadrature
    % Quad1D class Quad1D: Quadrature rule
    %   The Quad1D class is a quadrature rule in 1D 
    %
    % (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015
    methods(Access = public)   
        function obj = Quad1D
            %constructor
            obj = obj@Quadrature;
        end        
        function W = getWeights(obj)
            switch obj.m
                case 1
                    W = 2;
                case 2
                    W = [1,1];
                case 3
                    W = [8,5,5]/9;
                otherwise
                    error('this quadrature rule not yet implemented');
            end
        end
        function X = getReferenceQuadPoints(obj)
            switch obj.m
                case 1
                    X = 0;
                case 2
                    X = 1./[sqrt(3),-sqrt(3)];
                case 3
                    X = [0,sqrt(3/5),-sqrt(3/5)];
                otherwise
                    error('this quadrature rule not yet implemented');

            end
        end
        function X = getQuadPoints(obj,intval)
            X = obj.getReferenceQuadPoints;
            a = intval(1); b = intval(2);
            X = a + (X+1)*(b-a)/2;            
        end
        function Y = getQuadValues(obj,f,intval)
            X = obj.getQuadPoints(intval);            
            Y = f(X);
        end
    end
    methods(Static)
        function [P,Z] = trafo(intval)
            a = intval(1); 
            b = intval(2);
            P=a;
            Z = (b-a)/2 ;
        end
        function dt = detTrafo(intval)
            a = intval(1); 
            b = intval(2);
            
            dt = (b-a)/2;
        end
    end
end