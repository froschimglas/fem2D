function h = plotError(geo,u,fun)
% plots numerical solution over triangulation
% INPUT:geo - instance of geometry class 
%       u coefficient vector of linear fem solution
%       fun - optional variable: if given, then the pointwise difference 
%       | u-fun | is plotted         
% OUTPUT: h handle to figure 
%
% (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015

if nargin<3
    tri = geo.getTRep;
    figure
    h = trisurf(tri.Triangulation,geo.nodes(:,1),geo.nodes(:,2),u);
else
    tri = geo.getTRep;
    z = fun(geo.nodes(:,1), geo.nodes(:,2));
    figure
    hold on;
    h=trisurf(tri.Triangulation,geo.nodes(:,1),geo.nodes(:,2),abs(u-z));
    title('pointwise error');
    hold off;
end
end