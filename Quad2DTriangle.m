classdef Quad2DTriangle < Quadrature
    % Quad2DTriangle class Quad2DTriangle: Quadrature rule for triangles
    %   The Quad2D class is a quadrature rule in 2D 
    %
    % (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015
    
    methods(Access = public)   
        function obj = Quad2DTriangle(quadrule)
            % constructor            
            if nargin<1
                quadrule = 3;
            end
            
            obj = obj@Quadrature(quadrule);
            
        end        
        function W = getWeights(obj)
            switch obj.m
                case 1
                    W = 1;
                case 2
                    W = [1,1,1]/3;
                case 3
                    W = [-0.5625,0.52083333333333333333,0.52083333333333333333,0.52083333333333333333];                    
            end
        end
        function X = getReferenceQuadPoints(obj)
            switch obj.m
                case 1
                    X1 = 1/3;
                    X2 = 1/3;
                case 2
                    X1 = [0,  0.5,0.5];
                    X2 = [0.5,0,  0.5];
                case 3
                    X1 = [1/3,0.6, 0.2, 0.2];
                    X2 = [1/3,0.2,0.6, 0.2];
            end
            X = [X1;X2];
        end
        function X = getQuadPoints(obj,Triangle)
            Xref = obj.getReferenceQuadPoints;                        
            [P,Z] =  obj.trafo(Triangle);
            
            X = repmat(P,1,size(Xref,2)) + Z*Xref;            
            
        end
        function Y = getQuadValues(obj,f,Triangle)
             
            X = obj.getQuadPoints(Triangle);    
            X1 = X(1,:);
            X2 = X(2,:);
            
            Y =f(X1,X2);
        end        
    end
    methods(Static)
        function dt = detTrafo(Triangle)           
            % Trafo from Reference triangle to Triangle. Be carefule, a
            % factor 0.5 has been considered for the quadrature! If you
            % need the determinant of the Jacobian of the Trafo you need to
            % multiply dt by 2
            dt= ( Triangle.B(1)-Triangle.A(1) )*( Triangle.C(2)-Triangle.A(2) ) - ( Triangle.B(2) - Triangle.A(2) )*( Triangle.C(1)-Triangle.A(1) ) ;
            dt = dt/2;
        end
        function [P,Z] = trafo(Triangle)
            % trafo is P + Z*x
            A = Triangle.A;
            B = Triangle.B;
            C = Triangle.C;
            
            ba = B-A;
            ca = C-A;
            
            P= A';
            Z = [ ba', ca' ];
        end
    end
end