function tid=connected_to(tri,id)
% find all triangles the id-th point is connected to
[I,~]=find(tri==id);
tid=unique(I);
end