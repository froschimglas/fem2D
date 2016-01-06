function fval = trafof(xr,yr,f,P,Z)
% f( trafo( xref )) Transformation of a function f given over an arbitrary
% triangle to a reference triangle
% INPUT:
% xr, yr: points on the reference triangle
% f: function handle to the function f over the 'real' triangle
% P,Z: transformation parameters. The trafo from a reference triangle to a given triangle is realized via the
% formula x = P + Z*xref, where x is a point in the `real` triangle and
% xref a point in the reference triangle. P and Z must be given, for
% instance by the Gauss Quadrature rules.

% (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015



x = P(1)*ones(size(xr)) + Z(1,1)*xr + Z(1,2)*yr;
y = P(2)*ones(size(xr)) + Z(2,1)*xr + Z(2,2)*yr;

fval = f(x,y);
 

end
