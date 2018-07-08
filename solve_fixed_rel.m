function result=solve_fixed_rel(M,RHS,varargin)

% Alvin Wong, January 9, 2013
% Solve linear equations M*x=RHS with constraints given in varargin.

n=size(M,1);
id_all=1:n;
id_dep=[];

for i=1:length(varargin)
    rel=varargin{i};
    m=size(rel,1);
    p=size(rel,2)/2-1;
    for j=1:m
        id1=rel(j,1);
        for k=1:p
            id2=rel(j,1+2*k-1);
            c=rel(j,1+2*k);
            M(id2,:)=M(id2,:)+c*M(id1,:);
            M(:,id2)=M(:,id2)+c*M(:,id1);
            RHS(id2)=RHS(id2)+c*RHS(id1);
        end
        RHS=RHS-rel(j,1+2*p+1)*M(:,id1);
    end
    id_dep=[id_dep;rel(:,1)];
end

id_indep=setdiff(id_all,id_dep);

x=zeros(n,1);
x(id_indep)=M(id_indep,id_indep)\RHS(id_indep);

for i=1:length(varargin)
    rel=varargin{i};
    m=size(rel,1);
    p=size(rel,2)/2-1;
    for j=1:m
        id1=rel(j,1);
        for k=1:p
            id2=rel(j,1+2*k-1);
            c=rel(j,1+2*k);
            x(id1)=x(id1)+c*x(id2);
        end
        x(id1)=x(id1)+rel(j,1+2*p+1);
    end
end

result=x;

end