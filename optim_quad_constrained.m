function result=optim_quad_constrained(M,RHS,varargin)

% Alvin Wong, March 8, 2013
% Optimization of a quadratic function with linear constraints.
% Assume the gradient of the quadratic energy is given.
% Constraints are given in varargin as arguments of the form
% [c1 id1 c2 id2 ... ck idk c],
% meaning that c1*x(id1)+c2*x(id2)+...+ck*x(idk)+c=0.
% The problem is solved using the method of Lagrange multipliers.

n=size(M,1);
l=0;
N=sparse(n);
RHS2=zeros(n,1);
for i=1:length(varargin)
    rel=varargin{i};
    m=size(rel,1);
    p=(size(rel,2)-1)/2;
    for j=1:m
        % one more constraint
        l=l+1;
        for k=0:p-1
            N(n+l,rel(j,2*k+2))=rel(j,2*k+1);
            N(rel(j,2*k+2),n+l)=rel(j,2*k+1);
            RHS2(n+l)=rel(j,end);
        end
        
    end
end
RHS2(1:n)=RHS;
N(1:n,1:n)=M;

%tic
% r=symrcm(N);
% tmp=N(r,r)\RHS2(r);
% tmp(r)=tmp;
%toc

% tic
% p=symamd(N);
% tmp2=N(r,r)\RHS2(r);
% toc

% tic
tmp=N\RHS2;
% condest(N)
% toc

% n=length(tmp);
% tmp2=diag(N);
% tmp2(tmp2<0.000001)=0.000001;
% %tmp2=1./tmp2;
% %tmp2=sqrt(tmp2);
% tmp2=sparse(1:n,1:n,tmp2);
% N2=tmp2*N;
% tmp3=N2\RHS2;

result=tmp(1:n);
end
























