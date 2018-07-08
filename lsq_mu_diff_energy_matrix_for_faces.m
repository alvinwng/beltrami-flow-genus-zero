function M=lsq_mu_diff_energy_matrix_for_faces(pt,tri,A)
% Alvin Wong, April 10, 2013
% Given a triangular mesh, return the energy matrix for the least square
% energy of the Beltrami differential defined on faces.
% verified correct
if nargin==2
    A=mesh_area(pt,tri);
end
m=size(tri,1);
M=sparse((1:2*m)',(1:2*m)',[A;A]);
end