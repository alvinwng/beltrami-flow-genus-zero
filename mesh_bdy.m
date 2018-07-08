function boundary=mesh_bdy(face, options)

% UNFINISHED
% Given a triangular mesh, find the boundary of it.
% Useful for finding boundaries of point clouds with local mesh.

% compute_boundary - compute the vertices on the boundary of a 3D mesh
%
%   boundary=compute_boundary(face);
%
%   Copyright (c) 2007 Gabriel Peyre

addpath('../Fast BHF/toolbox_graph/toolbox');

options.null = 0;
verb = getoptions(options, 'verb', 1);

if size(face,1)<size(face,2)
    face=face';
end

nvert=max(max(face));
nface=size(face,1);

A=sparse(nvert,nvert);
for i=1:nface
    if verb
        progressbar(i,nface);
    end
    f=face(i,:);
    A(f(1),f(2))=A(f(1),f(2))+1;
    A(f(1),f(3))=A(f(1),f(3))+1;
    A(f(3),f(2))=A(f(3),f(2))+1;
end
A=A+A';

for i=1:nvert
    u=find(A(i,:)==1);
    if ~isempty(u)
        boundary=[i u(1)];
        break;
    end
end

s=boundary(2);
i=2;
while(i<=nvert)
    u=find(A(s,:)==1);
    if length(u)~=2
        warning('problem in boundary');
    end
    if u(1)==boundary(i-1)
        s=u(2);
    else
        s=u(1);
    end
    if s~=boundary(1)
        boundary=[boundary s];
    else
        break;
    end
    i=i+1;
end
       
if i>nvert
    warning('problem in boundary');
end

boundary=boundary(length(boundary):-1:1);