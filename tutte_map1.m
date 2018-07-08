function trg=tutte_map1(pt,tri,id_bdy,pt_bdy)

% Alvin Wong, April 20, 2014
% Compute a Tutte embedding with boundary points specified.

n=size(pt,1);
deg=mesh_degree(pt,tri)+1;

tmp=[tri(:,1) tri(:,2); tri(:,2) tri(:,1);
     tri(:,2) tri(:,3); tri(:,3) tri(:,2);
     tri(:,3) tri(:,1); tri(:,1) tri(:,3)];
tmp=unique(tmp,'rows');
I=[(1:n)';tmp(:,1)];
J=[(1:n)';tmp(:,2)];
S=[1-1./deg;-1./deg(tmp(:,1))];
M=sparse(I,J,S,n,n);
x=solve_fixed_rel(M,zeros(n,1),[id_bdy pt_bdy(:,1)]);
y=solve_fixed_rel(M,zeros(n,1),[id_bdy pt_bdy(:,2)]);

trg=[x,y];

end