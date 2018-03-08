function diagEN=blockones(C,N)
% this funcion will produce a block diag ones matrix
A=ones(N);
diagEN=zeros(C*N);
for ii=1:C
    mv=(ii-1)*N;
    diagEN(1+mv:N+mv,1+mv:N+mv)=A;
end

end %end of the function file