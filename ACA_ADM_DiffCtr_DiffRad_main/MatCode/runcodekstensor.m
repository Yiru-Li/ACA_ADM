function [Efield,x,te2p2,p2,te2te]=runcodeksTensor(te2p,p,conductivity,rs,js,ro,FEMord)
%te2p is 3 by nte
%p is 3 by np
%conductivity nte by 1
%rs is 3 by nsource
%js is 3 by nsource
%ro is 3 by nobs
%FEMord is the order of the FEM
%step 1 geometric preprocessing
% tic
[te2p2,p2]=femgenmesh_c(te2p,p,FEMord);
[te2te]=gente2te(te2p);
% T1=toc;
% disp(['Time Step 1: ',num2str(T1),' s']);
%% Step 2 assemble FEM matrix
% tic
A=femassembleTensor(te2p2,p2,conductivity,FEMord);
% T2=toc;
% disp(['Time Step 2: ',num2str(T2),' s']);
%% Step 3 generate right hand side of equation
% tic
[rhs]=femgenrhskstensor(te2p2,p2,conductivity,rs,js,FEMord);
% T3=toc;
% disp(['Time Step 3: ',num2str(T3),' s']);
%% Step 4 delete one unknown and equation and define preconditioner
A=A(1:end-1,1:end-1);
PRECON=sparse(1:numel(A(:,1)),1:numel(A(:,1)),sqrt(full(diag(A))));
rhs=rhs(1:end-1);
% T4=toc;
% disp(['Time Step 4: ',num2str(T4),' s']);
%% Step 5 solve system of equations
% tic
[x,flag,relres,iter,resvec,resveccg]=minres(A,rhs,10^-10,10000,PRECON,PRECON);
x(end+1)=0;
% T5=toc;
% disp(['Time Step 5: ',num2str(T5),' s']);
%% Step 6 evaluate field at desired locations
% tic
Efield=FEMinterpolatorks(te2te,te2p2,p2,rs,js,x,ro,FEMord);
% T6=toc;
% disp(['Time Step 6: ',num2str(T6),' s']);


