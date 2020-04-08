function F=buildF(NDimen,NY,NX,gpu)

%For fastDownsamplingFFT Fourier coeffs (only ups and downs)
F=cell(1,NDimen);

flM=mod(NY,2);
A=(NY-1)/2;
A=bsxfun(@plus,[ceil(A);NX-floor(A)],flM); %compute non-zero postions
for m=1:NDimen  
    F{m}=dynInd(single(complex(dftmtx(NX(m)))),[1:A(1,m) A(2,m):NX(m)],1); %fourier coefs
    if gpu>0;F{m}=gpuArray(F{m});end
    F{m}=F{m}/sqrt(NX(m)); %normalization
end