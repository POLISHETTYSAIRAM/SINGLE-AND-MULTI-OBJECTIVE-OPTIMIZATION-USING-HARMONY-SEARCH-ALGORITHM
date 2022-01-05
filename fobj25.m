

%% 25 bar optimization -- objective function
function Z=fobj25(x)
% The well-known spring design problem
% This program is only for 10 bar 2D truss problem 
% Input for 25 bar 3D truss
n=10; % no of nodes
ne=25; % no of elements
ndof=3; % no of DOF's
nen=2; % no of nodes for each element
nee=nen*ndof;
SectionsList=[1.62 1.8 1.99 2.13 2.38 2.62 2.63 2.88 2.93 3.09 3.13 3.38 3.47 3.55 3.63 3.84 3.87 3.88 4.18 4.22 4.49 4.59 4.8 4.97 5.12 5.74 7.22 7.97 11.5 13.5 13.9 14.2 15.5 16 16.9 18.8 19.9 22 22.9 26.5 30 33.5];
ar=SectionsList(ceil(x)); % area of the elements in sq.m
fixed_dof=[19;20;21;22;23;24;25;26;27;28;29;30]; % constrained DOF's
L=[1 1 2 1 2 2 2 1 1 3 4 3 5 3 6 4 5 4 3 5 6 6 3 4 5;2 4 3 5 6 4 5 3 6 6 5 4 6 10 7 9 8 7 8 10 9 10 7 8 9]; % element connecting matrix
coord=[-40 40 -40 40 40 -40 -100 100 100 -100;0 0 40 40 -40 -40 100 100 -100 -100;200 200 100 100 100 100 0 0 0 0]; % coordinate vector for the 6 nodes in m
load=[1;-10;-10;0;-10;-10;0.5;0;0;0;0;0;0;0;0;0.6;0;0];  
% load vector at 10th dof and 12th dof i.e.downwards at 5th & 6th nodes respectively in N
E=10000; % youngs modulus in Pa
den=0.1; % density of material in Kg/m3
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
d=K1\load;
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
% fprintf('***********************************************************************************************\n');
% fprintf('\n');
% fprintf('RESULTS\n');
% fprintf('*******\n');
% fprintf('\n');
% fprintf('\n');
% fprintf('DISPLACEMENT AT THE NODES \n');
% fprintf('**********************************\n');
% for i=1:n
% fprintf('\n');
% fprintf('NODE'); disp(i);
% disp(D_big(:,i));
% fprintf('**********************************\n');
% fprintf('\n');
% end
Z = (TW)+sum(10^6*max(abs(disp_vec)-(0.35)*ones(30,1),zeros(30,1)))+sum(40*max(axial_vec-(150*10^6)*ones(25,1),zeros(25,1)));
% fprintf('***********************************************************************************************\n');
% fprintf('\n');
% fprintf('AXIAL STRESSES IN EACH OF THE MEMEBERS\n');
% disp(axial_vec);
% fprintf('***********************************************************************************************\n');
% figure(1);
% qw=500*D_big+coord;
% for i=1:ne
% X=[qw(1,L(1,i)) qw(1,L(2,i))];
% Y=[qw(2,L(1,i)) qw(2,L(2,i))];
% if ndof==3
% Z=[qw(3,L(1,i)) qw(3,L(2,i))];
% line(X,Y,Z);
% else
% line(X,Y);
% end
% hold on
% end
% hold on
% for i=1:ne
% X=[coord(1,L(1,i)) coord(1,L(2,i))];
% Y=[coord(2,L(1,i)) coord(2,L(2,i))];
% if ndof==3
% Z=[coord(3,L(1,i)) coord(3,L(2,i))];
% line(X,Y,Z,'color',[1 0 0]);
% else
% line(X,Y,'color',[1 0 0]);
% end
% hold on
% end
% gtext('RED - UNDEFORMED,BLUE - DEFORMED');
%z=z+getnonlinear(u);
% function Z=getnonlinear(u)
% Z=0;
% % Penalty constant
% lam=10^15;
% 
% % % Inequality constraints
% % g(1)=1-u(2)^3*u(3)/(71785*u(1)^4);
% % gtmp=(4*u(2)^2-u(1)*u(2))/(12566*(u(2)*u(1)^3-u(1)^4));
% % g(2)= gtmp+1/(5108*u(1)^2)-1;
% % g(3)=1-140.45*u(1)/(u(2)^2*u(3));
% % g(4)=(u(1)+u(2))/1.5-1;
% 
% % No equality constraint in this problem, so empty;
% geq=[];
% 
% % Apply inequality constraints
% for k=1:length(g),
% Z=Z+ lam*g(k)^2*getH(g(k));
% end
% % Apply equality constraints
% for k=1:length(geq),
% Z=Z+lam*geq(k)^2*getHeq(geq(k));
% end
% 
% % Test if inequalities hold
% % Index function H(g) for inequalities
% function H=getH(g)
% if g<=0,
% H=0;
% else
% H=1;
% end
% % Index function for equalities
% function H=getHeq(geq)
% if geq==0,
% H=0;
% else
% H=1;
% end
% ----------------- end ------------------------------
end








