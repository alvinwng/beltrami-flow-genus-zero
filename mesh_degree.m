
% Alvin Wong, November 25, 2012

% compute the degree of each vertex of a triangular mesh
% checked and should be correct

function deg=mesh_degree(pt,tri)
i1=tri(:,1);
i2=tri(:,2);
i3=tri(:,3);
edges=unique([i1 i2; i2 i1; i2 i3; i3 i2; i3 i1; i1 i3],'rows');
n=size(pt,1);
tmp=sparse(edges(:,1),ones(size(edges,1),1),ones(size(edges,1),1),n,1);
deg=full(tmp);
end