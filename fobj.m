%% frame optimization -- objective function
 function Z=fobj(x)

% The well-known spring design problem
% This program is only for 10 bar 2D truss problem

% Input for 52 bar 2D truss
% Input for 24 bar 2D truss

n=18; % no of nodes
ne=24; % no of elements

ndof=2; % no of DOF's
nen=2; % no of nodes for each element
nee=nen*ndof;
SectionsList= ((0.0254)^2)*[257,238,81,73,65.7,275,251,236,213,193];
c=1;
ar=SectionsList(ceil(x)); % area of the elements in sq.
% for i=1:8
%     for j=1:3
% ar(c)=ar1(i);
% c=c+1;
%     end
% end

fixed_dof=[1;2;3;4]; % constrained DOF's
L=[1,2,3,4,3,5,6,5,7,8,7,9,10,9,11,12,11,13,14,13,15,16,15,17;3,4,5,6,4,7,8,6,9,10,8,11,12,10,13,14,12,15,16,14,17,18,16,18]; % element connecting matrix
coord=[0,3.4,0,3.4,0,3.4,0,3.4,0,3.4,0,3.4,0,3.4,0,3.4,0,3.4;0,0,3.4,3.4,6.8,6.8,10.2,10.2,13.6,13.6,17,17,20.4,20.4,23.8,23.8,27.2,27.2];
% coordinate vector for the 6 nodes in m
load= 1000*[1.2,0,0,0,2.42,0,0,0,3.63,0,0,0,4.839,0,0,0,6.054,0,0,0,7.264,0,0,0,8.743,0,0,0,12.529,0,0,0]; %load vector at 10th dof and 12th dof i.e.downwards at 5th & 6th nodes respectively in N
E=2.07e11; % youngs modulus in Pa
den=7860*9.81; % density of material in Kg/m3

%% Calculation
%getting the ID matrix from the information given
ID=zeros(ndof,n);
r=size(ID);
c=r(1,2);
r=r(1,1);
number=1;
for c1=1:c
for r1=1:r
ID(r1,c1)=number;
number=number+1;
end
end

fixed_dof=sort(fixed_dof);

fixed1=fixed_dof;
node=1:1:n*ndof;
free_dof=node;
for i=1:length(fixed_dof)
free_dof(fixed1(i))=[];
fixed1=fixed1-1;
end

clear fixed1;
free_dof=sort(free_dof);
free_dof=free_dof';
for e=1:ne
for i=1:nen
for a=1:ndof
p=ndof*(i-1)+a;
LM(p,e)=ID(a,L(i,e));
end
end
end
K=zeros(n*ndof);
lambda_vec=[];
h_vec=[];
%assembley
for e=1:ne %coordinates matrix of each element
localcoord=[coord(:,L(1,e)) coord(:,L(2,e))];
h=0;
for count=1:ndof %length of member
temp=(localcoord(count,2)-localcoord(count,1))^2;
h=h+temp;
end
h=sqrt(h);
h_vec=[h_vec;h];
lambda=[];
for count=1:ndof %direction cosines in lambda matrix like cosx
lambda=[lambda;(localcoord(count,2)-localcoord(count,1))/h];

%first one is lambda x ans second one is lambda y and so on
end
lambda_vec=[lambda_vec lambda];
A=lambda*lambda';
k=(E*ar(e))/h*[A -A;-A A];
for p=1:nee
P=LM(p,e);
for q=1:nee
Q=LM(q,e);

K(P,Q)=K(P,Q)+k(p,q); % GLobal Stiffness matrix
end
end
end
lambda_vec=lambda_vec';
K1=K;
%applying boundary conditions
for counter=1:length(fixed_dof)
K1(fixed_dof(counter),:)=[];
K1(:,fixed_dof(counter))=[];
fixed_dof=fixed_dof-1;
end
load=load';
d=inv(K1)*load;
disp_vec=zeros(n*ndof,1);
for count=1:length(free_dof)
disp_vec(free_dof(count))=d(count);
end

%% Weight
TW=0;
for i=1:ne
W=den*ar(i)*h_vec(i);
TW=TW+W; % weight in Kg
end

%% post processing
D_big=[];
count=1;
for i=1:n
D=[];
for j=1:ndof
t=disp_vec(count);
count=count+1;
D=[D;t];
end
D_big=[D_big D];

end
axial_vec=[];
force_vec=[];
for i=1:ne
d1=lambda_vec(i,:)*D_big(:,L(1,i));
d2=lambda_vec(i,:)*D_big(:,L(2,i));
axial=(E/h_vec(i))*(d2-d1);
axial_vec=[axial_vec;axial];
force=ar(i)*axial_vec(i);
force_vec=[force_vec;force];
end
axial_vec;
%
% fprintf('***********************************************************************************************\n');
Z = (TW+sum(10^6*max(abs(disp_vec)-(5*10^-3)*ones(36,1),zeros(36,1)))+sum(10^6*max(axial_vec-(150*10^6)*ones(24,1),zeros(24,1))))/9.81;