function F=IRWLSmaterial(kv,dfL,MTs,lims,pexp,toler)

tstart=tic;
MTs=MTs(lims(1,1):lims(2,1),lims(1,2):lims(2,2))>0.5;
%Hay que usar calculo vectorial para acelerar...

%kv->I*2
%dfL->NR*NC*I*NT*2
N=size(dfL);
iter=1;
tolant=0;
tolact=1e6;
if N(4)~=1
    dpL=zeros([N(1:2) N(4) N(5) N(3)]);
    for s=1:N(3)
        for t=1:N(4)      
            for u=1:N(5)
                dpL(:,:,t,u,s)=dfL(:,:,s,t,u);
            end
        end
    end
else
    dfL=squeeze(dfL);
    dpL=zeros([N(1:2) N(5) N(3)]);
    for s=1:N(3)    
        for u=1:N(5)
            dpL(:,:,u,s)=dfL(:,:,s,u);
        end
    end
end
M=size(dpL);
if N(4)~=1
    W=zeros([M(1:3) M(5) M(5)]);
    F=zeros([M(1:3) 2 2]);
else
    W=zeros([M(1:2) M(4) M(4)]);
    F=zeros([M(1:2) 2 2]);
end
Fant=F;
if N(4)~=1
    F(:,:,:,1,1)=1;
    F(:,:,:,2,2)=1;
    Mdir=repmat(kv,[1 1 M(1:3)]);
    MinvBis=zeros([M(5) M(4) M(1:3)]);
    dpLBis=shiftdim(dpL,3);
else
    F(:,:,1,1)=1;
    F(:,:,2,2)=1;
    Mdir=repmat(kv,[1 1 M(1:2)]);
    MinvBis=zeros([M(4) M(3) M(1:2)]);
    dpLBis=shiftdim(dpL,2);
end
Mdir=shiftdim(Mdir,2);

while tolact>toler && ( (pexp~=2 && N(3)>3) || iter~=2) % && abs(tolant-tolact)>toler 
    for i=1:N(3)        
        if N(4)~=1
            W(:,:,:,i,i)=((Mdir(:,:,:,i,1)-sum(dpL(:,:,:,:,i).*F(:,:,:,:,1),4)).^2+(Mdir(:,:,:,i,2)-sum(dpL(:,:,:,:,i).*F(:,:,:,:,2),4)).^2).^((pexp-2)/2);
        else
            W(:,:,i,i)=((Mdir(:,:,i,1)-sum(dpL(:,:,:,i).*F(:,:,:,1),3)).^2+(Mdir(:,:,i,2)-sum(dpL(:,:,:,i).*F(:,:,:,2),3)).^2).^((pexp-2)/2);
        end
    end      
    W(W>1e9)=1e9;
    if N(4)~=1
        WBis=shiftdim(W,3);
        for m=1:M(1)        
            for n=1:M(2)            
                if MTs(m,n)>0
                    for o=1:M(3)                
                        Waux=WBis(:,:,m,n,o);
                        Y=dpLBis(:,:,m,n,o);
                        kvpWaux=Y*Waux;
                        kvpWauxkv=kvpWaux*Y';
                        MinvAux=mldivide(kvpWauxkv,kvpWaux);      
                        MinvBis(:,:,m,n,o)=MinvAux';            
                    end
                end
            end
        end
    else
        WBis=shiftdim(W,2);
        for m=1:M(1)        
            for n=1:M(2)
                if MTs(m,n)>0
                    Waux=WBis(:,:,m,n);
                    Y=dpLBis(:,:,m,n);
                    kvpWaux=Y*Waux;
                    kvpWauxkv=kvpWaux*Y';
                    MinvAux=mldivide(kvpWauxkv,kvpWaux);      
                    MinvBis(:,:,m,n)=MinvAux';
                end
            end
        end
    end
    Minv=shiftdim(MinvBis,2);
    tolant=tolact;
    if N(4)~=1
        for s=1:2
            for r=1:2
                F(:,:,:,s,r)=sum(Minv(:,:,:,:,s).*Mdir(:,:,:,:,r),4);
            end
        end
        tolact=max(max(max(max(max(abs(Fant-F))))));
    else
        for s=1:2
            for r=1:2
                F(:,:,s,r)=sum(Minv(:,:,:,s).*Mdir(:,:,:,r),3);
            end
        end
        tolact=max(max(max(max(abs(Fant-F)))));
    end    
    Fant=F;
    iter=iter+1;
    %[iter tolact]
end

%figure(11)
%for n=1:2   
%    subtightplot(2,1,n,[0 0])
%    imshow(F(:,:,n,1),[-2 2])
%end
%set(gcf, 'Position', get(0,'Screensize'))
%%export_fig(sprintf('/home/lcorgra/Escritorio/Articulos/02EnProceso/MRM14/figs/Diagrams/diagr2_im6.png'));  
%export_fig(sprintf('/home/lcg13/Articulos/02EnProceso/MRM14/figs/Diagrams/diagr2_im6.png'));  

%figure(12)
%for n=1:2   
%    subtightplot(2,1,n,[0 0])
%    imshow(F(:,:,n,2),[-2 2])
%end
%set(gcf, 'Position', get(0,'Screensize'))
%%export_fig(sprintf('/home/lcorgra/Escritorio/Articulos/02EnProceso/MRM14/figs/Diagrams/diagr2_im7.png')); 
%export_fig(sprintf('/home/lcg13/Articulos/02EnProceso/MRM14/figs/Diagrams/diagr2_im7.png')); 

M=size(F);
if N(4)==1
    F=reshape(F,[M(1:2) 1 M(3:4)]);
end
    
telapsed=toc(tstart);
tirwls=sprintf('Tiempo IRWLS %d orientaciones: %f',N(3),telapsed)

%[U S V]=svd(kv,0);  
%MatrizInversion=V*(S^-1)*U';%Equivalente a: ((kv'*kv)^-1)*kv';