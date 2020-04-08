function F=IRWLSMaterial(kv,dfL,MTs,lims,toler,pexp)

%IRWLSMATERIAL   Estimates the material deformation gradient tensor 
%   F=IRWLSMATERIAL(KV,DFL,MTS,LIMS,TOLER,PEXP) estimates the material
%   deformation gradient tensor form the wrapped gradient of the local
%   phase
%   KV are the set of wave vectors (2nd dimension) for different 
%   orientations (1st dimension) 
%   DFL is the wrapped gradient of the local phase
%   MTS is the MR-C mask in the MR-T space
%   LIMS is the limitant region for the myocardium
%   TOLER is the tolerance threshold for convergence of the IRWLS method
%   PEXP is the exponent of the fitting criteria (2 for LS / 1 for LAD)
%   Defaults to 1
%   It returns:
%   F, the estimated material deformation gradient tensor
%

if ~exist('pexp','var')
    pexp=1;
end

tstart=tic;
dfL=double(dfL);
MTs=MTs(lims(1,1):lims(2,1),lims(1,2):lims(2,2))>0.5;

iter=1;
tolant=0;
tolact=1e6;
dfL=permute(dfL,[1 2 4 5 3]);

M=size(dfL);
W=double(zeros([M(1:3) M(5) M(5)]));
F=double(zeros([M(1:3) 2 2]));
Fant=F;
F(:,:,:,1,1)=1;
F(:,:,:,2,2)=1;
Mdir=repmat(double(kv),[1 1 M(1:3)]);
MinvBis=double(zeros([M(5) M(4) M(1:3)]));
dfLBis=shiftdim(dfL,3);
Mdir=permute(Mdir,[3 4 5 1 2]);

while tolact>toler && ( (pexp~=2 && M(5)>3) || iter~=2) % && abs(tolant-tolact)>toler 
    for i=1:M(5)    
        W(:,:,:,i,i)=((Mdir(:,:,:,i,1)-sum(dfL(:,:,:,:,i).*F(:,:,:,:,1),4)).^2+(Mdir(:,:,:,i,2)-sum(dfL(:,:,:,:,i).*F(:,:,:,:,2),4)).^2).^((pexp-2)/2);
    end      
    W(W>1e9)=1e9;
    WBis=shiftdim(W,3);
    for m=1:M(1)        
        for n=1:M(2)            
            if MTs(m,n)>0
                for o=1:M(3)                
                    Waux=WBis(:,:,m,n,o);
                    Y=dfLBis(:,:,m,n,o);
                    kvpWaux=Y*Waux;
                    kvpWauxkv=kvpWaux*Y';
                    MinvAux=mldivide(kvpWauxkv,kvpWaux);      
                    MinvBis(:,:,m,n,o)=MinvAux';            
                end
            end
        end
    end
    
    Minv=permute(MinvBis,[3 4 5 1 2]);
    tolant=tolact;
    for s=1:2
        for r=1:2
            F(:,:,:,s,r)=sum(Minv(:,:,:,:,s).*Mdir(:,:,:,:,r),4);
        end
    end
    tolact=max(abs(Fant(:)-F(:)));   
    Fant=F;
    iter=iter+1;
end
F=single(F);

telapsed=toc(tstart);
fprintf('IRWLS %d orientations time: %f\n',M(5),telapsed)

%[U S V]=svd(kv,0);  
%InversionMatrix=V*(S^-1)*U';%Equivalent to: ((kv'*kv)^-1)*kv';