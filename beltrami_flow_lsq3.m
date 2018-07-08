function [u_dot,v_dot]=beltrami_flow_lsq3(pt,tri,outward,u_vec,v_vec,nu,A,varargin)

% Alvin Wong, February 15, 2013
% Compute the Beltrami flow on a 3D mesh with respect to Beltrami
% differential nu.

n=size(pt,1);
% A=mesh_area(pt,tri);
% least square energy for discrete Beltrami flow
[M,C]=beltrami_flow_lsq3_energy_grad(pt,tri,outward,u_vec,v_vec,real(nu),imag(nu),A);
% result=solve_fixed_rel(M,-C,varargin{:});

l=length(varargin);
myvarargin={};
for i=1:l
    k=size(varargin{i},2);    % even number
    tmp_id=reshape(1:k,2,k/2);
    tmp_id=[tmp_id(2,:); tmp_id(1,:)];
    tmp_id=reshape(tmp_id,k,1);
    
    tmp=varargin{i};
    myvarargin{i}=[tmp(:,1) repmat(-1,size(tmp,1),1) tmp(:,2:end)];
    myvarargin{i}(:,1:k)=myvarargin{i}(:,tmp_id);
end
result=optim_quad_constrained(M,-C,myvarargin{:});

u_dot=result(1:n);
v_dot=result(n+1:2*n);

end