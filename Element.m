classdef Element
% Element class Element
%   The Element class represents a local element.
%
% (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015
    
    %=============================
    properties(Access = public)
        % Vertices [A,B,C]
        Vertices
    end
    methods(Access = public)        
        % Constructor
        function obj = Element(vertices)
            obj.Vertices = vertices;
            
        end
        
        % Queries
        function x = perimeter(obj)
            % returns perimeter of the triangle
            x = norm( obj.a) + norm( obj.b) + norm( obj.c );
        end
        function x = area(obj)
            % returns area of the triangle
            s = obj.perimeter / 2;
            x = s*(s-norm( obj.a))*(s-norm(obj.b))*(s-norm(obj.c));
        end
        function X = barycenter(obj)
            % returns barycenter of this element
           X = (obj.A+obj.B+obj.C)/3;
        end
        function x = diameter(obj)
            % returns diameter of the incircle of this triangle
            %
            % mesh size parameter h of this cell is incircle of cell            
            x = 2*obj.area / obj.perimeter;
        end
        
        % Nodes [A,B,C]
        function X = A(obj,i)
            % returns first vertex, if optional parameter i is given, the
            % ith component of A
            if nargin==1
                X = obj.vertex(1);
            else
                X = obj.Vertices(1,i);
            end
        end
        function X = B(obj,i)
            % returns second vertex, if optional parameter i is given, the
            % ith component of B     
            if nargin==1
                X = obj.vertex(2);
            else
                X = obj.Vertices(2,i);
            end
        end
        function X = C(obj,i)
            % returns third vertex, if optional parameter i is given, the
            % ith component of C

            if nargin==1
                X = obj.vertex(3);
            else
                X = obj.Vertices(3,i);
            end
        end
        
        function X = vertex(obj,i)
           % returns ith vertex

            if 0<i<4
                X = obj.Vertices(i,:);
            else
                error('Error: A triangle has only three vertices!');
            end                    
        end
        
        % Edges [a,b,c]
        function x = a(obj)
            % edge B-to-C
            x = obj.C-obj.B;
        end
        function x = b(obj)
            % edge C-to-A
            x = obj.A-obj.C;
        end
        function x = c(obj)
            % edge A-to-B
            x = obj.B-obj.A;
        end
        
        function x = edge(obj,i)
            % returns ith edge vector where i is local number of edge [a,b,c]
            switch i
                case 1
                    x = obj.a;
                case 2
                    x = obj.b;
                case 3
                    x = obj.c;
                otherwise
                    error('Error: A triangle has only three edges!');
            end
        end
        function i_max = longestEdge(obj)
            % returns length and index of longest edge i.e max ( a,b,c )
            d_edge = zeros(3,1);
            for i=1:3
                d_edge(i) = norm( obj.edge(i) );
            end
            [~,i_max] = max(d_edge);
        end
        function X = normal(obj,i)
            % returns the outward normal on the ith edge
            Y = obj.edge(i);
            X = [Y(2);-Y(1)]/ norm(Y);
        end
        
        % Bisection
        function [E1,VertexOrder1,E2,VertexOrder2,Bnew] =  bisect(obj,edgeIndex)
            % bisects edge given by edgeIndex
            % OUTPUT: new Elements E1 and E2
            % with Indices of Vertices VertexOrder1 and VertexOrder2, where
            % 4 is the index of the new node Bnew on edge edgeIndex.
            
            % start vertex
            si = obj.edgeStart(edgeIndex);
            Anew = obj.vertex(si);
            
            % new vertex            
            Bnew = Anew + obj.edge(edgeIndex)/2;
            
            % new end vertex
            N = obj.edgeNodes(edgeIndex);
            ci = setdiff([1,2,3],N);
            Cnew = obj.vertex(ci);
            ei = setdiff(N,si);
            V = obj.vertex(ei);
            
            % new elements            
            E1 = Element([Anew; Bnew; Cnew]);
            VertexOrder1 = [si,4,ci];
            E2 = Element([Bnew; V; Cnew]);
            VertexOrder2 = [4,ei,ci];
        end
        
        % Plot
        function plot(obj,withNumbers)
            % plots the triangle, optional parameter withNumbers
            % default=false, if true, numbers of vertices are printed as
            % well as the barycenter
            X = obj.Vertices(:,1);
            X = [X;obj.Vertices(1,1)];
                
            Y = obj.Vertices(:,2);
            Y = [Y;obj.Vertices(1,2)];
            
            plot(X,Y,'-*');
            if nargin > 1 && withNumbers
                hold on;
                
                for i=1:3
                    NC=obj.vertex(i);
                    stri = int2str(i);
                    text(NC(1),NC(2),stri,'FontSize',18);
                end
                
                NC = obj.barycenter;
                plot(NC(1),NC(2),'*r');
                hold off;
            end
        end
    end
    methods(Static)
        function N = edgeNodes(i)
            % local numbers of nodes on this edge
             N = setdiff([1,2,3],i);
        end
        function N = edgeStart(i)
            % node at the start of edge - with counter clock wise
            % numbering!
            startVs = [2,3,1];
            N = startVs(i);
        end
        function e = edgeBetween(edgeNodes)
            % determines from two nodes with local numbers the local edge
            % number of the edge between them
            i=edgeNodes(1); j=edgeNodes(2);
            A = [0, 3, 2; 3,0,1; 2,1,0];
            e = A(i,j);
        end
    end
end