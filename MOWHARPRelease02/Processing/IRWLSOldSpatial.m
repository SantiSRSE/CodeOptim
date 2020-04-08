function f=IRWLSspatial(kv,dfL,pexp,toler)

tstart=tic;
%Hay que usar limitantes para hacerlo mas rapido->quiza mejor en la
%llamada e incluso quiza mejor de manera vectorial...

%X->I*2
%dfL->NR*NC*I*NT*2
N=size(dfL);
iter=1;
tolant=0;
tolact=1e6;
if N(4)~=1
    dpL=zeros([N(1:2) N(4) N(3) N(5)]);
    for s=1:N(4)
        for t=1:N(3)       
            dpL(:,:,s,t,:)=dfL(:,:,t,s,:);
        end
    end
else
    dfL=squeeze(dfL);
    dpL=dfL;
end
M=size(dpL);
if N(4)~=1
    W=zeros([M(1:4) M(4)]);
    f=zeros([M(1:3) 2 2]);
else
    W=zeros([M(1:3) M(3)]);
    f=zeros([M(1:2) 2 2]);
end
fant=f;  
if N(4)~=1
    f(:,:,:,1,1)=1;
    f(:,:,:,2,2)=1;
    Mdir=repmat(kv',[1 1 M(1:3)]);
    MinvBis=zeros([M(4) M(5) M(1:3)]);
else
    f(:,:,1,1)=1;
    f(:,:,2,2)=1;
    Mdir=repmat(kv',[1 1 M(1:2)]);
    MinvBis=zeros([M(3) M(4) M(1:2)]);
end
Mdir=shiftdim(Mdir,2);
    
while tolact>toler && ( (pexp~=2 && N(3)>3) || iter~=2) % && abs(tolant-tolact)>toler 
    for i=1:N(3)        
        if N(4)~=1
            W(:,:,:,i,i)=((dpL(:,:,:,i,1)-sum(Mdir(:,:,:,:,i).*f(:,:,:,:,1),4)).^2+(dpL(:,:,:,i,2)-sum(Mdir(:,:,:,:,i).*f(:,:,:,:,2),4)).^2).^((pexp-2)/2);
        else
            W(:,:,i,i)=((dpL(:,:,i,1)-sum(Mdir(:,:,:,i).*f(:,:,:,1),3)).^2+(dpL(:,:,i,2)-sum(Mdir(:,:,:,i).*f(:,:,:,2),3)).^2).^((pexp-2)/2);
        end
    end      
    W(W>1e9)=1e9;
    if N(4)~=1
        WBis=shiftdim(W,3);
        for m=1:M(1)        
            for n=1:M(2) 
                for o=1:M(3)                
                    Waux=WBis(:,:,m,n,o);
                    kvpWaux=kv'*Waux;
                    kvpWauxkv=kvpWaux*kv;
                    MinvAux=mldivide(kvpWauxkv,kvpWaux);      
                    MinvBis(:,:,m,n,o)=MinvAux';            
                end
            end
        end
    else
        WBis=shiftdim(W,2);
        for m=1:M(1)        
            for n=1:M(2)               
                Waux=WBis(:,:,m,n);
                kvpWaux=kv'*Waux;
                kvpWauxkv=kvpWaux*kv;
                MinvAux=mldivide(kvpWauxkv,kvpWaux);      
                MinvBis(:,:,m,n)=MinvAux';                            
            end
        end
    end        
    Minv=shiftdim(MinvBis,2);
    tolant=tolact;
    if N(4)~=1
        for s=1:2
            for r=1:2
                f(:,:,:,s,r)=sum(Minv(:,:,:,:,s).*dpL(:,:,:,:,r),4);
            end
        end
        tolact=max(max(max(max(max(abs(fant-f))))));
    else
        for s=1:2
            for r=1:2
                f(:,:,s,r)=sum(Minv(:,:,:,s).*dpL(:,:,:,r),3);
            end
        end
        tolact=max(max(max(max(abs(fant-f)))));
    end
    fant=f;
    iter=iter+1;
    [iter tolact]
end
M=size(f);
if N(4)==1
    f=reshape(f,[M(1:2) 1 M(3:4)]);
end

telapsed=toc(tstart);
tirwls=sprintf('Tiempo IRWLS %d orientaciones: %f',N(3),telapsed)

%[U S V]=svd(kv,0);  
%MatrizInversion=V*(S^-1)*U';%Equivalente a: ((kv'*kv)^-1)*kv';