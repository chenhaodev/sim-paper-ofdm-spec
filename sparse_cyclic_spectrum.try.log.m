%%Failed Trial
%==========

% sparse reconstruction 2
% cannot diff signal / noise
Pre = Phi'*Sy*Phi;
D = (dftmtx(cs.N))^(-1);
for ii = 1:cs.N
	[hat(:,ii), ~] = cosamp(Pre(:,ii), D, cs.sparse, cs.iter);
end
for jj = 1:cs.N
	[W(jj,:), ~] = cosamp(hat(jj,:)', D, cs.sparse, cs.iter);
end
W = fftshift(W);
Spec_f_cs = abs(W);

% sparse reconstruction 2.1
% cannot diff signal / noise
%Pre = Phi'*Sy*Phi;
D = (dftmtx(cs.N))^(-1);
A = Phi*D;
for ii = 1:cs.M
	[hat(:,ii), ~] = cosamp(Sy(:,ii), A, cs.sparse, cs.iter);
end
for jj = 1:cs.N
	[W(jj,:), ~] = cosamp(hat(jj,:)', A, cs.sparse, cs.iter);
end
W = fftshift(W);
Spec_f_cs = abs(W);

% sparse reconstruction 3, fails
% recon fail
S = DRD;
R = inv(D) S inv(D)
Rz = ARA'
Rz = A inv(D) S inv(D) A' ;
Rz * P = A inv(D) S, where P = pinv (inv(D) * A);

% sparse recon 4, vec{} and kron{}, fails since memory full.

% sparse recon 5, alternatively cosamp:
% cannot diff signal / noise
for ii = 1:cs.M
    A= Phi*inv(D);
    for iter = 1:2
        if iter == 1    
            b = Sy(:,ii); 
        else
            b = Sy(ii,:)';
        end
        cvx_begin
            variable x(N);
            minimize(norm(x,1));
            A*x == b;
        cvx_end
        if iter == 1    
            hat1(:,ii) = x';
        else
            hat2(ii,:) = x;
        end 
    end
	%[hat1(:,ii), ~] = cosamp(Sy(:,ii), Phi*inv(D), 32, 100);
	%[hat2(ii,:), ~] = cosamp(Sy(ii,:)', Phi*inv(D), 32, 100);
end
W1 = fftshift(hat1);
W2 = fftshift(hat2);
W = W1 * W2;
Spec_f_cs = abs(W);
