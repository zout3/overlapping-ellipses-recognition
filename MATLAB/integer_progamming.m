function ind = integer_progamming(distMat, lambda)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the optimal combination of ellipses
% distMat: distance matrix between P contour segments and Q ellipses
% lambda: regularization term
% ind: index of the selected optimal ellipses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = size(distMat,1);
N = size(distMat,2);
f = [reshape(distMat, [M*N,1]);lambda*ones(N,1)];

A = zeros(2*N, M*N+N);
b = zeros(2*N, 1);
for i=1:N
    A(i,(1:M)+M*(i-1)) = -1;
    A(i,M*N+i) = 1;
    
    A(N+i,(1:M)+M*(i-1)) = 1;
    A(N+i,M*N+i) = -M;
end

Aeq = zeros(M, M*N+N);
beq = ones(M, 1);
for i=1:M
    Aeq(i,(0:(N-1))*M+i) = 1;
end

lb = zeros(M*N+N,1);
ub = ones(M*N+N,1);
intcon = 1:(M*N+N);

options = optimoptions('intlinprog','Display','off');
x = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options);
ind = find(round(x((M*N+1):end)));
end


