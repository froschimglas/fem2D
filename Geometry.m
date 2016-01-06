classdef Geometry
% Geometry class Geometry
%   The Geometry class partly uses distmesh to generate triangulations of
%   domains.   
%
% (c) Daniela Fusseder, Technische Universität Kaiserslautern, 2015
    
    %=============================
    properties(Access = public)
        % nodes for linear finite element mesh
        nodes;
        % element list stores index of 3 nodes
        elements
        % mesh size parameter
        h;
        % shape / geometry
        shape
        % special geometry properties
        geoProperties
    end
    properties(Access = private, Hidden=true)
        tr; % Matlab TriRep 
        BCEdgeMatrix; % matrix of nodes of edges on boundary
    end
    methods(Access = public)    
        function obj = Geometry(paramh)
            % constructor
            % INPUT:             
            %  paramh:      optional mesh size parameter. default is 0.25.
            %               A smaller parameter results in a finer mesh. It is used both 
            %               by the distgen mesh generator and the rectangle
            %               triangulation of this class. 
            if nargin>0            
                % paramh must be between 0 and 1
                if paramh<0 || paramh>1
                    error('mesh size parameter must be in (0,1]');
                else
                    obj.h = paramh;                
                end
            else
                obj.h = 0.25;
            end
            
            obj.shape = 'square';
        end
        
    end
    methods
        % SET METHODS
        function obj = set.shape(obj, newshape)            
           %reload data after input was renewed
            obj = obj.getData(newshape);
            obj.shape = newshape;
        end
        function obj = getData(obj,newshape)
            switch newshape
                case 'circle'
                    [obj.nodes,obj.elements] = getCircle(obj.h);
                case 'rectangle'
                    [obj.nodes,obj.elements] = getRectangle(obj.h);
                case 'square'
                    [obj.nodes,obj.elements] = getSquare;
                otherwise
            end
            X = obj.nodes(:,1);
            Y = obj.nodes(:,2);
            obj.tr = TriRep(obj.elements,X,Y);
            
            Edges = obj.getBCEdges;
            obj.BCEdgeMatrix = sparse(Edges(:,1),Edges(:,2),ones(size(Edges,1),1),obj.nmbNodes,obj.nmbNodes,size(Edges,1));
        end     
        
        % GET METHODS
        function n = nmbNodes(obj)
            % returns number of nodes (dofs)
            n = size(obj.nodes,1);
        end
        function n = nmbElements(obj)
            % returns number of elements (triangles, cells)
            n = size( obj.elements,1);
        end
        function bcEdges = getBCEdges( obj )
            % returns a list of edges which are on the boundary of the
            % geometry
            bcEdges =freeBoundary(obj.tr);
        end
        function TR = getTRep(obj)
            % returns the MATLAB TriRep object of this triangulation
            TR = obj.tr;
        end
        
        % ELEMENT METHODS
        function globNodes = globalNodeIndices( obj,elem )
            % returns the three indices of the nodes of this element
            % (element has index elem)
            globNodes = obj.elements(elem,:);
        end
        function coordNodes = getNodes( obj,globElementNodes )
            % returns the coordinates of the 3 nodes of element with index
            % globElementNodes, each row of the return value corresponds to one node
            coordNodes = zeros(3,2);
            for i=1:3
                coordNodes(i,:) = obj.nodes( globElementNodes(i),: );
            end
        end
        function elementVertices = getVertices(obj,elem)
            % returns coordinates of Element Nodes for elem with global
            % index elem
             globNodes = obj.globalNodeIndices( elem );
             elementVertices = obj.getNodes(globNodes);
        end
        function nodesOnEdge = getEdgeVertices(obj,edgeIndices)
            e = edges(obj.tr);
            nodesOnEdge = zeros( length(edgeIndices)*2,3 );
            
            for i=1:length(edgeIndices)            
                edgeIndex = edgeIndices(i);
                N = e(edgeIndex,:);
                
                N1.index = N(1);
                N2.index = N(2);
                
                N1.coordinates = obj.nodes(N(1),:);
                N2.coordinates = obj.nodes(N(2),:);
                
                nodesOnEdge(2*i-1,:)    = [N1.index,N1.coordinates];
                nodesOnEdge(2*i,:)      = [N2.index,N2.coordinates];                             
            end
            nodesOnEdge = sortrows(nodesOnEdge,2);
            nodesOnEdge = unique(nodesOnEdge,'rows','stable');
            nodesOnEdge = nodesOnEdge(:,2:3);
        end
        function locElem = getLocalElement(obj,elem)
            elemNodes = obj.getVertices(elem);
            locElem = Element(elemNodes);
        end
        function CP = getCenterPoint(obj,elem)
            % compute barycenter of this element
            locElem = obj.getLocalElement(elem);
            CP = locElem.barycenter;
            
        end
        function isOnBoundary = isBoundaryEdge( obj,elemIndex,localEdgeNodes )
            % check if the edge with the local edge nodes is a boundary
            % edge
            
            % get global node indices
            nodeIndices = obj.globalNodeIndices(elemIndex);
            
            E1 = localEdgeNodes(1);   % local node 1 of bisected edge
            E2 = localEdgeNodes(2);   % local node 2 of bisected edge
            
            E1_node = nodeIndices(E1);  % global node 1 of bisected edge
            E2_node = nodeIndices(E2);  % global node 2 of bisected edge
            
            isOnBoundary = true;
            if obj.BCEdgeMatrix(E1_node,E2_node) == 0 && obj.BCEdgeMatrix(E2_node,E1_node) == 0
                isOnBoundary = false;
            end
            
        end
        
        
        % REFINEMENT METHODS
        function refGeo = refine(obj,markedElements)
            if nargin<2
                % refine all
                markedElements = obj.nmbElements:-1:1;
            end
            
            refGeo = obj;
            
             while ~isempty(markedElements)
                 K = markedElements(1);
                 [refGeo,new_indices,refinedElements] = refineSimple(refGeo,K);
                 % delete refined elements from marked-element-list
                 markedElements = setdiff(markedElements,refinedElements);
                                  
                 % translate numbering of elements to new numbers
                 markedElements = new_indices(markedElements);
                 markedElements = sort(markedElements,'descend');                                  
             end
            
        end
        
        % PLOT METHODS
        function plot(obj,withNumbers)
            % plots the mesh, if input paramter withNumbers=true, the
            % element numbers are printed as well
            triplot(obj.tr);
            if nargin>1
                if withNumbers
                    hold on;                                                                                 
                                           
                    for i=1:obj.nmbElements
                        NC=obj.getCenterPoint(i);
                        stri = int2str(i);
                        text(NC(1),NC(2),stri,'FontSize',18);
                    end
                    hold off;
                end
            end
        end
        function plotEdges(obj)
            triplot(obj.tr);
            
            e = edges(obj.tr);
            
            
            
            hold on;
            
            for i=1:size(e,1)
                N= e(i,:);
                
                N1 = obj.nodes( N(1),: );
                N2 = obj.nodes( N(2),: );
                
                NC = N1 + (N2-N1)/2;
                
                stri = int2str(i);
                text(NC(1),NC(2),stri,'FontSize',15);
            end
            hold off;
        end
    end
end

%% Refinement
function [newGeo,new_indices,refinedElements] = refineSimple(G,markedElement)

refElem = markedElement;
[refinedElements,newNodes,newElements]= refineElement(G,refElem);
[newGeo,new_indices] = addElements(G,newNodes,newElements,refinedElements);

end
function [refinedElements,newNodes,newElements] =  refineElement(G,markedElement)
elem = markedElement;

% switch to local elements
elemNodes = G.getVertices(elem);
refElem = Element(elemNodes);

% refine local element
i_max = refElem.longestEdge;
[~,locE1,~,locE2,newNode] = refElem.bisect(i_max);

% check if it is on boundary
N = refElem.edgeNodes(i_max);
isOnBoundary = G.isBoundaryEdge( elem,N );

% case 1: edge is on boundary
% -> two new elements, one new knot
refinedElements = elem;
newNodes = newNode;
nodeOrder = [G.globalNodeIndices(elem),G.nmbNodes+1];
newElements = [ nodeOrder(locE1); nodeOrder(locE2)];

% case 2: edge is shared with neighbour
% -> four new elements, one new knot

% get neighbour and refine it too
if ~isOnBoundary
    SI = edgeAttachments(G.getTRep, nodeOrder(N));
    neighElem = setdiff(SI{1},elem);
    
    % in simple refinement: don't bisect largest edge of this
    % element, but the edge shared with first element!
    
    % switch to local elements
    globNodes2 = G.globalNodeIndices(neighElem);
    Nglob = nodeOrder(N);
    edgeN2 = [ find(globNodes2==Nglob(1)), find(globNodes2==Nglob(2)) ];
    elemNodes2 = G.getVertices(neighElem);
    refElem2 = Element(elemNodes2);
    
    % refine local element
    ref_edge = refElem2.edgeBetween( edgeN2 );
    [~,locE1,~,locE2,~] = refElem2.bisect(ref_edge);
    
    refinedElements = [refinedElements,neighElem];
    nodeOrder = [globNodes2, G.nmbNodes+1];
    newElements = [ newElements; ...
        nodeOrder(locE1); ...
        nodeOrder(locE2)];
end
end
function [newGeo,old_indices] = addElements(G,newNode,newElements,oldElements)
newGeo = G;
old_indices = 1:G.nmbElements;

% get old element list
prev_elements = newGeo.elements;

% delete the refined elements
prev_elements(oldElements,:)=[];

% save the new indices of the remaining elements
for i=1:length(oldElements)
    k = oldElements(i);
    old_indices(k:end) = old_indices(k:end)-1;
end

% add the new nodes
newGeo.nodes = [newGeo.nodes; newNode];

% add the new elements
newGeo.elements = [prev_elements;newElements];

% get the new triangulation trirepresentation
X = newGeo.nodes(:,1);
Y = newGeo.nodes(:,2);
newGeo.tr = TriRep(newGeo.elements,X,Y);

% update the edge matrix for boundary edges
Edges = newGeo.getBCEdges;
newGeo.BCEdgeMatrix = sparse( Edges(:,1),Edges(:,2),ones(size(Edges,1),1),newGeo.nmbNodes,newGeo.nmbNodes,size(Edges,1));
end


%% Geometry data
function [Nodes,Elements] = getSquare()
Nodes = [...
    0 0;
    1 0;
    0.5,0.5;
    0 1;
    1 1];

Elements = [ ...
    1 2 3;
    1 3 4;
    3 5 4;
    3 2 5];
end
function [Nodes,Elements] = getRectangle(h)
a = 1;
b = -1;
c = 1;
d = -1;


nX = ceil(1/h)+1;
nY = ceil(1/h)+1;

x = linspace(a,b,nX);
y = linspace(c,d,nY);


% create nodes
Nodes = zeros(nX*nY,2);
for i=1:nX
    for j=1:nY
        Nodes(i + (j-1)*nX,1)=x(i);
        Nodes(i + (j-1)*nX,2)=y(j);
    end
end

% create triangulation
Elements = zeros( (nX-1)*(nY-1)*2 ,3);
elem = 1;
for i=1:nX-1
    for j=1:nY-1
        % in rectangle (i,j)
        % there are 4 nodes
        
        N1 = (j-1)*nX + i;
        N2 = (j-1)*nX + i+1;
        N3 = j*nX + i;
        N4 = j*nX + i+1;
        
        Elements(elem,:) = [N1,N2,N3];
        Elements(elem+1,:) = [N4,N3,N2];
        elem = elem+2;
    end
end
end
function [Nodes,Elements] = getCircle(h)
fd=inline('sqrt(sum(p.^2,2))-1','p');
[Nodes,Elements]=distmesh2d(fd,@huniform,h,[-1,-1;1,1],[]);
end
