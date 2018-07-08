% Alvin Wong, revised on June 2, 2018

tot_steps=500;
TOL=10^(-5);

[pt,~,tri]=read_obj('./data/fandisk.obj');

% remove duplicate points
[pt,~,n]=unique(pt,'rows');
tri=n(tri);

tic

tri=[tri(:,1) tri(:,3) tri(:,2)];
n=size(pt,1);
m=size(tri,1);

[trg,vid,tid]=tutte_map2(pt,tri,30);
% figure; p1=view_mesh(trg,tri); % initial map

trg=trg*6;
trg=stereo_proj_inv(trg(:,1),trg(:,2));
trg(vid,:)=[0 0 -1];

v=[-0.1536,-0.1763,0.9723];
scale=15;

trg=rotate_space(trg,v);
[tmpx,tmpy]=stereo_proj(trg);
tmpx=tmpx*scale;
tmpy=tmpy*scale;
trg=stereo_proj_inv(tmpx,tmpy);

id_0=find(abs(trg(:,3)-1)==min(abs(trg(:,3)-1)));
id_1=find(abs(trg(:,2)-1)==min(abs(trg(:,2)-1)));
id_inf=find(abs(trg(:,3)+1)==min(abs(trg(:,3)+1)));
id_01inf=[id_0 id_1 id_inf];

mu=zeros(m,1);

A=ones(m,1)/m;

error_sup=zeros(tot_steps,1);
error_l2=zeros(tot_steps,1);
tol=zeros(tot_steps,1);

outward=compute_sphere_face_outward(trg,tri);
mu_now=compute_beltrami_differential(pt,trg,tri,outward);

M=lsq_mu_diff_energy_matrix_for_faces(pt,tri,A);
result.energy0=[real(mu_now);imag(mu_now)]'*M*[real(mu_now);imag(mu_now)];

started_positive_step=0;

for i=1:tot_steps
    disp(i);
    % compute u_vec and v_vec (on tangent plane)
    [r_vec,s_vec]=compute_sphere_uv_vec(trg);
    
    fz=compute_fz_mesh3(pt,trg,tri,outward);
    
    nu=(mu-mu_now).*fz.^2./abs(fz).^2./(1-abs(mu_now).^2);
    
    [r_dot,s_dot]=beltrami_flow_lsq3(trg,tri,outward,r_vec,s_vec,nu,A,[id_01inf' zeros(3,1); n+id_01inf' zeros(3,1);]);
    
    % determine the optimal step size from x_dot and y_dot
    [Et,Ett,t_max]=deriv_mu_energy_3d(trg,tri,r_vec,s_vec,r_dot,s_dot,M,mu_now,outward,pt);
    if -Et/Ett>1
        t_opt=min(1,min(t_max)*0.5);
        started_positive_step=1;
    elseif -Et/Ett>0
        t_opt=min(-Et/Ett,min(t_max)*0.5);
        started_positive_step=1;
    else
        if started_positive_step==0 || i<20
        	t_opt=min(0.5,min(t_max)*0.5);
        else
            t_opt=min(-Et/Ett,min(t_max)*0.5);
        end
    end
    
    trg=trg+t_opt*(r_vec.*[r_dot r_dot r_dot]+s_vec.*[s_dot s_dot s_dot]);
    
    tmp=sqrt(trg(:,1).^2+trg(:,2).^2+trg(:,3).^2);
    tmp=trg./[tmp tmp tmp];
    tmp_shift=tmp-trg;
    trg=tmp;
    
    % store the optimal step size for iteration i
    result.t_opt(i)=t_opt;
    % compute tolerence
    tmp=t_opt*(r_vec.*[r_dot r_dot r_dot]+s_vec.*[s_dot s_dot s_dot])+tmp_shift;
    result.tol(i)=max(sqrt(tmp(:,1).^2+tmp(:,2).^2+tmp(:,3).^2));
    % compute errors
    outward=compute_sphere_face_outward(trg,tri);
    mu_now=compute_beltrami_differential(pt,trg,tri,outward);
    result.sup_mu(i)=max(abs(mu-mu_now));
    result.energy(i)=[real(mu_now);imag(mu_now)]'*M*[real(mu_now);imag(mu_now)];
    
    % display restuls for iteration i
    disp('Optimal stepsize:');
    disp(result.t_opt(i));
    disp('Tolerence:');
    disp(result.tol(i));
    disp('Error in sup-norm:');
    disp(result.sup_mu(i));
    disp('Error in L2-norm:');
    disp(result.energy(i));
    
    if result.tol(i)<TOL
        break;
    end
end
toc

p1=view_mesh(pt,tri);
figure;
p2=view_mesh(trg,tri);
