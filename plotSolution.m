function h = plotSolution(geo,u,userTitle,fun)
% plots numerical solution over triangulation
% INPUT:geo -instance of geometry class
%       u coefficient vector of linear fem solution
%       userTitle: optional input, if given this will be the title used in
%       the plot
%       fun optional input, if given both the solution u and fun are
%       plotted in one figure
% OUTPUT: h handle to figure 
%
% (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015
if nargin<4
    tri = geo.getTRep;
    figure
    clf
    h = trisurf(tri.Triangulation,geo.nodes(:,1),geo.nodes(:,2),u);
    shading('interp');lighting flat;
    if nargin<3
        title('solution');
    else
        title(userTitle);
    end
else
    tri = obj.geometry.getTRep;
    z = fun(obj.geometry.nodes(:,1), obj.geometry.nodes(:,2));
    figure
    hold on;
    trisurf(tri.Triangulation,obj.geometry.nodes(:,1),obj.geometry.nodes(:,2),z);
    h=trisurf(tri.Triangulation,obj.geometry.nodes(:,1),obj.geometry.nodes(:,2),u);
    hold off;
end
end

