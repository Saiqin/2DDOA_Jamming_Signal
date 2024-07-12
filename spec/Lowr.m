function B=Lowr(A,K)
%��K����
[U,S,V] = svd(A,'econ');
[~,ii]=sort(abs(diag(S)),'descend');
ii=ii(1:K);
B=U(:,ii)*S(ii,ii)*V(:,ii)';
end