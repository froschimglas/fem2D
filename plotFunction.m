function h = plotFunction(geo,fun)
% plots a function over the triangulation
% INPUT: geo - instance of geometry class
%        fun - function handle to the map over geometry
% OUTPUT: h handle to figure
%
% (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015

tri = geo.getTRep;
z = fun(geo.nodes(:,1), obj.geometry.nodes(:,2));
figure
h = trisurf(tri.Triangulation,obj.geometry.nodes(:,1),obj.geometry.nodes(:,2),z);
end

