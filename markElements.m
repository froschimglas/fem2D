function elements = markElements(elemErr,markMethod)
% Marks elements in a mesh for refinement by user specified methods
% INPUT: elemErr - vector of error on each element
%        markMethod - method by which elements should be marked for
%        refinement. Currently only the method max is available
%
% TODO: other mark methods, optional parameter to compute threshold
%
% (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015

[eMax,iMax] = max(elemErr);

switch markMethod
    case 'max'
        threshold = 0.75*eMax;
        elements = find( elemErr > threshold );
        
    otherwise
        error('No such marking method');
end


end