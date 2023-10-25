function [P,T,Edge] = getMesh_Gmsh(filename)

% P1 element
% input: 'filename'
% output:
%         P -- the coordinate of point
%         T -- the information of element (T(1:3,:) is the nodes index) (T(4,:) --> which domain belongs to)
%         Edge -- the nodes of boundary edges
% Delete $PhysicalNames ~ $EndPhysicalNames


%%
file = fopen(filename, 'r'); % fid means fileID, fopen:（以只读方式）打�?文件或获得有关打�?文件的信�?

% 读取xxxx.msh 的前4�?, not interested
buf = fscanf(file, '%s', 1); % fscanf: 读取文本文件中的数据�?'$MeshFormat'
buf = fscanf(file, '%s',3);
buf = fscanf(file, '%s', 1); % '$EndMeshFormat'
buf = fscanf(file, '%s', 1); % '$Nodes'


% Read the number of nodes and Coordinate information
N_Nodes = fscanf(file, '%i', 1); % number of nodes, 2575
nodes = fscanf(file, '%i %f %f %f', [4, N_Nodes]); % All the nodes,  fscanf 按列顺序填充 A。sizeA 必须为正整数或采�? [m n] 的形�?
nodes = nodes(2:3,:)'; % 只需�?(x,y)并转�?

buf = fscanf(file, '%s', 1); % '$EndNodes'
buf = fscanf(file, '%s', 1); % '$Elements'

% Read the number of element and element information
nElem = fscanf(file, '%i', 1); % number of elements including boundary edges(104), 5148=104+5044
% elements, first three rows are nodes, 4th row -- phyical region
% 5th row -- belong to which subdomain, 6th row -- possible neighboring
% subdomain if nonzero, if there are several?
elem = zeros(6,nElem); % memory allocation and initialization
elemTag = zeros(100,1);
Edge = [];
n = 0; % number of elements, only triangle elements

% 读取单元信息
for i = 1 : nElem
    elemInfo = fscanf(file, '%i %i %i', [3, 1]); % [1;1;2]
    for ii = 1 : elemInfo(3)
       elemTag(ii) = fscanf(file, '%i', 1);
    end
    % boundary edges, maybe not useful, num = 104
    if elemInfo(2) == 1
        edgeNodes = fscanf(file,'%i %i',[2,1]); 
        edgeNodes = [edgeNodes;elemTag(2)];  %edgeNodes(3):edgeTag
        Edge = [Edge edgeNodes]; % append [5; 77]
    end
    % triangle elements, num = 5044
    if elemInfo(2) == 2
        n = n + 1;
        elemNodes = fscanf(file,'%i %i %i', [3, 1]); % [196;430;331]
        elem(1:3,n) = elemNodes; % 单元顶点
        elem(4,n) = elemTag(1);  % �?属区�?
        elem(5,n) = elemTag(4);
        if elemInfo(3) > 4
            elem(6,n) = elemTag(5);
        end
    end
end
%
elems = elem(:,1:n)';
%
P = nodes';
T = elems';