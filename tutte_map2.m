function [trg,vid,tid]=tutte_map2(pt,tri,vid)

% Alvin Wong, April 21, 2014
% Compute a Tutte embedding mapping the neighbors of vertex vid onto a unit
% disk.
% If vid is not specified, choose a vertex with the maximum number of
% neighbors.
% Return the Tutte map and the vertex chosen.

if nargin==2
    deg=mesh_degree(pt,tri);
    vid=find(deg==max(deg),1);      % a vertex with maximal degree
end

tid=connected_to(tri,vid);      % all triangles joining vid

[tri2,b]=tidy_tri(tri(tid,:));

vid_connected=mesh_bdy(tri2)';

if tri(tid(1),1)==vid
    vid_a=tri(tid(1),2);
    vid_b=tri(tid(1),3);
elseif tri(tid(1),2)==vid
    vid_a=tri(tid(1),3);
    vid_b=tri(tid(1),1);
else
    vid_a=tri(tid(1),1);
    vid_b=tri(tid(1),2);
end
% vertices a, b are in counter-clockwise order w.r.t. vid

tmp1=find(vid_connected==vid_a);
tmp2=mod(tmp1+1,length(vid_connected));
if tmp2==0
    tmp2=length(vid_connected);
end
if vid_connected(tmp2)==vid_b
    % boundary positively oriented w.r.t. vid
    vid_connected=vid_connected(end:-1:1);
else
    % boundary negatively oriented w.r.t. vid
end

% we want vid_connected to be negatively oriented w.r.t. vid

k=length(vid_connected);

trg=tutte_map1(pt,tri,b(vid_connected),[cos((0:(k-1))/k*2*pi)' sin((0:(k-1))/k*2*pi)']);

end

function [tri2,b]=tidy_tri(tri)
k=size(tri,1);
[b,m,n]=unique(tri);
tri2=reshape(n,k,3);
end