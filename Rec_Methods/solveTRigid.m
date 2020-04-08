function [TBfin,xTfin,rec,EE]=solveT_rig(x,y,rec,rig,ts,term,datos,T0,PuntosC,gpu)
%Code for rigid alignment 
%first version


NX=size(x);NX(end+1:3)=1;
NY=rec.NY;
NS=rec.NS;
NT=size(rig.T);


%Second order derivatives
a=[1 2 3 1 1 2 1 2 3 1 2 3 1 2 3 4 4 5 4 5 6;1 2 3 2 3 3 4 4 4 5 5 5 6 6 6 5 6 6 4 5 6];
NHe=size(a,2); 



dHe=single(zeros([NHe NT(5)]));
dH=single(zeros([NT(6) NT(5)]));

Eprev=single(zeros(NT(1:5)));
E=single(zeros(NT(1:5)));
    
if gpu>0
    Eprev=gpuArray(Eprev);E=gpuArray(E);x=gpuArray(x);y=gpuArray(y);rec.M=gpuArray(rec.M);rec.A=gpuArray(rec.A);rec.S=gpuArray(rec.S);
end


%Iterations
for n=1:rec.nT    
    xPrev=x;    
    %Update the weights            
    rec.w(rec.flagw==2)=rec.w(rec.flagw==2)*multB; %x1 no update
    rec.w(rec.flagw==1)=rec.w(rec.flagw==1)/multA;

    for s=1:NRun                                  
        [et,etg,eth]=precomputationsSincRigidTransform(rig.kGrid,rig.kkGrid,rig.rkGrid,rig.T(:,:,:,:,vS{s},:),1,2);
        if gpu
            et{1}=gpuArray(et{1});
            for m=1:3
                etg{1}{m}=gpuArray(etg{1}{m});
                for r=2:3
                    et{r}{m}=gpuArray(et{r}{m});etg{r}{m}=gpuArray(etg{r}{m});eth{r}{m}=gpuArray(eth{r}{m});
                end
            end
            for m=1:6
                eth{1}{m}=gpuArray(eth{1}{m});
            end
        end
        [xT,xB]=sincRigidTransform(x,et,1,gpu); 
        xT=forwardApplication(xT,rec,ts,Rigid,gpu);
        %xT=forwardApplication(xT,Dimen,FullEncode);
        xT=bsxfun(@minus,xT,y);
        xT=bsxfun(@times,xT,rec.A(:,:,:,:,vS{s}));
        xTconj=conj(xT);
        Eprev(1,1,1,1,vS{s})=sum(sum(sum(sum(real(xT.*xTconj),1),2),3),4);        

        [G,GB,GC]=sincRigidTransformGradient(xB,et,etg,gpu);
        for m=1:NT(6)
            G{m}=forwardApplication(G{m},rec,ts,Rigid,gpu);
            %G{m}=forwardApplication(G{m},Dimen,FullEncode);
            G{m}=bsxfun(@times,G{m},rec.A(:,:,:,:,vS{s}));
            Gconj{m}=conj(G{m});        
        end
        
        for m=1:NHe               
            GG=sincRigidTransformHessian(xB,GB,GC,et,etg,eth,m,gpu); 
            GG=forwardApplication(GG,rec,ts,Rigid,gpu);
            %GG=forwardApplication(GG,Dimen,FullEncode);
            GG=bsxfun(@times,GG,A(:,:,:,:,vS{s}));
            GG=real(GG.*xTconj);
            GG=GG+real(G{a(1,m)}.*Gconj{a(2,m)});
            if ~gpu
                dHe(m,vS{s})=permute(sum(sum(sum(sum(GG,1),2),3),4),[1 5 3 2 4]);
            else
                dHe(m,vS{s})=gather(permute(sum(sum(sum(sum(GG,1),2),3),4),[1 5 3 2 4]));
            end
        end
        
        for m=1:NT(6)             
            G{m}=real(bsxfun(@times,G{m},xTconj));                                            
            if ~gpu
                dH(m,vS{s})=permute(sum(sum(sum(sum(G{m},1),2),3),4),[1 5 3 2 4]);                       
            else
                dH(m,vS{s})=gather(permute(sum(sum(sum(sum(G{m},1),2),3),4),[1 5 3 2 4]));
            end
        end
    end
    
    MHe=single(eye(NT(6)));      
    for s=1:NT(5)
        for k=1:NHe
            if a(1,k)==a(2,k)
                MHe(a(1,k),a(2,k))=dHe(k,s)+rec.w(1,1,1,1,s);                            
            else
               MHe(a(1,k),a(2,k))=dHe(k,s);
               MHe(a(2,k),a(1,k))=dHe(k,s);
            end              
        end    
        dH(:,s)=single(double(MHe)\double(dH(:,s)));
    end    
    Tup=permute(dH,[3 4 5 6 2 1]);
    Tup=rig.T-Tup;    
    %Restrict to allowed ranges to prevent overshooting of Newton's method
    Tang=Tup(:,:,:,:,:,4:6);
    while ~isempty(Tang(Tang>pi))
       Tang(Tang>pi)=Tang(Tang>pi)-2*pi;
    end
    while ~isempty(Tang(Tang<-pi))
       Tang(Tang<-pi)=Tang(Tang<-pi)+2*pi;
    end
    Tup(:,:,:,:,:,4:6)=Tang;
    for m=1:3
        Ttra=Tup(:,:,:,:,:,m);
        if NX(m)>1
            while ~isempty(Ttra(Ttra>rGrid{m}(end)))
               Ttra(Ttra>rGrid{m}(end))=Ttra(Ttra>rGrid{m}(end))-NX(m);
            end
            while ~isempty(Ttra(Ttra<rGrid{m}(1)))
               Ttra(Ttra<rGrid{m}(1))=Ttra(Ttra<rGrid{m}(1))+NX(m);
            end
        end
        Tup(:,:,:,:,:,m)=Ttra;
    end   
    
    
    
    for s=1:NRun
        et=precomputationsSincRigidTransform(rig.kGrid,rig.kkGrid,rig.rkGrid,Tup(:,:,:,:,vS{s},:),1,0);
        if gpu
            et{1}=gpuArray(et{1});
            for m=1:3
                for r=2:3
                    et{r}{m}=gpuArray(et{r}{m});
                end
            end
        end
        xT=sincRigidTransform(x,et,1,gpu);     
        xT=forwardApplication(xT,rec,ts,Rigid,gpu);
        %xT=forwardApplication(xT,Dimen,FullEncode);                
        xT=bsxfun(@minus,xT,y);        
        xT=bsxfun(@times,xT,rec.A(:,:,:,:,vS{s}));
        xTconj=conj(xT);
        E(1,1,1,1,vS{s})=sum(sum(sum(sum(real(xT.*xTconj),1),2),3),4);
    end
    
    
    
    rec.flagw(E<Eprev)=2;
    rec.flagw(E>=Eprev)=1;
    for s1=1:size(Tup,6)
        TauxA=rig.T(:,:,:,:,:,s1);
        TauxB=Tup(:,:,:,:,:,s1);
        TauxA(rec.flagw==2)=TauxB(rec.flagw==2);
        rig.T(:,:,:,:,:,s1)=TauxA;
    end   
        
       
    %This would diminish drifting, as suggested in ec. (28) in the paper
    if meanT
        Tmed=mean(rig.T,5);
        etD=precomputationsSincRigidTransform(rig.kGrid,rig.kkGrid,rig.rkGrid,Tmed,1,0);
        if gpu
            etD{1}=gpuArray(etD{1});
            for m=1:3
                for ns=2:3
                    etD{ns}{m}=gpuArray(etD{ns}{m});
                end
            end
        end    
        x=sincRigidTransform(xPrev,etD,1,gpu);   
        x=rec.M.*x;
        rig.T=bsxfun(@minus,rig.T,Tmed);    
    end  
 
end



end
