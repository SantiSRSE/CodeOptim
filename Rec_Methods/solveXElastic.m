function [x, yEnd,rsnew]=solveXElastic(rec,gpu)

%SOLVEX   Reconstructs an image using CG SENSE with multishot alignment
%   X=SOLVEX(X,Y,M,T,S,SCONJ,PRECOND,A,KGRID,RKGRID,NX,GPU,BLOCKSIZE) 
%   computes the best X for a given T
%   X is the image to be reconstructed
%   Y is the measured data

%   rec reconstruction struct

%   T are the transform parameters and Tinv its inverse
%   GPU is a flag that determines whether to use gpu (1) or cpu (0) 
%   computation
%   BLOCKSIZE indicates the number of shots to be processed in a chunk
%   It returns:
%   X, the reconstructed image
%

verbose=0;
x=rec.x;
%filter parameters
NX=size(x);NX(end+1:3)=1;
h=single(buildFilter(size(x),'FractionalFiniteDiscreteIsoNorm',1,gpu,rec.ParX.RegFiltOrd)); %Filtro
%figure;imshow(fftshift(abs(h)),[])
%pause
g_map=1;

if gpu>0 %go to gpu
    x=gpuArray(x);rec.A=gpuArray(rec.A);rec.S=gpuArray(rec.S); rec.y=gpuArray(rec.y);
    rec.M=gpuArray(rec.M);rec.ParX.Precond=gpuArray(rec.ParX.Precond); h=gpuArray(h);
    rec.Tinv=gpuArray(rec.Tinv); rec.T=gpuArray(rec.T);
else
    x=gather(x); rec.S=gather(rec.S); rec.y=gather(rec.y);
    rec.ParX.Precond=gather(rec.ParX.Precond); h=gather(h);
    rec.Tinv=gather(rec.Tinv); rec.T=gather(rec.T);
end
rec.M=gather(rec.M); rec.A=gather(rec.A);

if gpu>0
    yEnd=gpuArray(single(zeros(NX))); 
else
    yEnd=single(zeros(NX)); 
end

%EhE solution
%yEnd = single(Encoding_forw_back( y,T,Tinv,rec,gpu,0 ));
yEnd=single(EncodingElastic(rec,rec.y,rec.Tinv,0));


y=yEnd;
y=rec.M.*y; %masking

Ap=applyCG(x);
r=y-Ap; %residuals
clear y 
z=rec.ParX.Precond.*r; %preconditioner (initial aprox to inverse)

p=z; 
rsold=sum(sum(sum(conj(z).*r)));

fina=0;
%Iterations
n=1;
while n<rec.ParX.nXTotal %Comjugate gradient
    Ap=applyCG(p);
    al=conj(rsold)/sum(sum(sum(conj(p).*Ap)));
    xup=al*p;        
    x=x+xup;
%    if verbose, figure(1), imshow(abs(x(:,:,30,:,1)),[]),pause(0.01), end
    
    xup=real(xup.*conj(xup));
    xup=max(max(max(xup)));
    if verbose, xup, end
    if xup<rec.ParX.toler && n>=rec.ParX.nX
        break
    end
    r=r-al*Ap;
    %z=r;
    z=rec.ParX.Precond.*r;
    rsnew=sum(sum(sum(conj(z).*r)));
    be=rsnew/rsold;
    p=z+be*p;
    rsold=rsnew; 
    if verbose, rsnew, end
    if sqrt(abs(rsnew))<1e-3
        break;
    end
    n=n+1;
end

if gpu==1
    x=gather(x);
end

%%
function x=applyCG(x)
    
    xOrig=x;
    if gpu
        xEnd=gpuArray(single(zeros(NX(1:3))));
    else 
        xEnd=single(zeros(NX(1:3)));
    end   
    
    %enconding
    %[ y2 ] = Encoding_forw_back( x,T,Tinv,rec,gpu,1 );
    y2=EncodingElastic(rec,x,rec.T,1);
    %decoding
    %[ xEnd ] = Encoding_forw_back( y2,T,Tinv,rec,gpu,0 );
    xEnd=EncodingElastic(rec,y2,rec.Tinv,0);
    
    x=xEnd;    
    x=x.*rec.M;
    %Regularization
    %x=x+lambda*xOrig;   %REgularizacion Tikhonov
    %x=x+0.1*xOrig;   %REgularizacion Tikhonov lambda optima
    x=x+rec.ParX.lambda.*g_map.*filtering(xOrig,conj(h).*h);   

end

end
