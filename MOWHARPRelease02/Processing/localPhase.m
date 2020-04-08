function [fL,limsb]=localPhase(ImT,phaseEstMethod,BW,sT,kT,orT,lims,Tred,NW,gpu)

%LOCALPHASE   Estimates the local phase 
%   [FL,LIMSB]=LOCALPHASE(BW,ST,KT,ORT,IMT,PHASEESTMETHOD,LIMS,TRED,NW,GPU)
%   estimates the local phase by applying extended HARP techniques
%   IMT is the the MR-T image (orientations in the 3rd dimension and
%   cardiac phases in the 4th dimension)
%   PHASEESTMETHOD is the phase estimation method to be applied. One of the
%   following: 'FT', 'WFT'
%   BW is the normalized bandwidth
%   ST is the MR-T resolution
%   KT is the spacing of the tags
%   ORT is the orientation of the tags (in radians)
%   LIMS is the limitant region for the myocardium
%   TRED are the rounded estimated translations (in pixels) to align the
%   MR-T images at different orientations
%   NW is the window width for the WFT
%   GPU is a flag that determines whether to use gpu (1) or cpu (0) 
%   computation
%   It returns:
%   FL, the estimated local phase
%   LIMSB, the updated limitant region for the myocardium
%

tstart=tic;

Coef=2
NT=size(ImT);
NT(end+1:4)=1;
rRange=0.15; %0.4 single-peak
aRange=pi/6;
mp=[1 Coef 1 Coef]; mo=[0 0 pi/2 pi/2];

maxTred=single(zeros(1,2));
for s=1:2                
    maxTred(s)=round(max(abs(Tred(s,:)/sT(s))));
    lims(1,s)=lims(1,s)-1-maxTred(s);    
    lims(2,s)=lims(2,s)+1+maxTred(s);    
end

ak=single(zeros(NT(3),2));
for s=1:NT(3)
    ak(s,:)=[-sin(orT(s)) cos(orT(s))];%Orientation of the spectral peak
end

%Apodization
alpha=2;
WAP=single(gausswin(NT(1),alpha));
WAP=WAP*WAP';
ImT=bsxfun(@times,ImT,WAP);

if strcmp(phaseEstMethod,'FT')
    dp=0;%Spectral peak is taken as the prescribed pattern frequency
    or=NT(1:2)/2+1;%Center of the k-space in pixels
    rk=NT(1:2).*sT/kT(1);%Radius of the spectral peak
    pk=bsxfun(@times,rk,ak);%k-space position of the spectral peak
    ik=bsxfun(@plus,or,pk);%k-space position of the spectral peak in pixels 
    
    rpeak=[((1-rRange)*rk(1)).^2 ((1+rRange)*rk(1)).^2];%Range for ridge searching in the radial direction    
    Xpeak=peakROI(or,NT,orT,rpeak,aRange);%ROI to search the peak
    fL(:,:,:,:,1)=fft2(ImT);       
    [peak,rtmax]=ridgeDetector(fL,Xpeak,dp,ik,rk,BW);
    v=single(zeros([NT 2]));   
    for t=1:NT(4)
        for s=1:NT(3)
            [v(:,:,s,t,1),v(:,:,s,t,2)]=ndgrid([1-peak{1}(s,t):NT(1)-peak{1}(s,t)].^2,[1-peak{2}(s,t):NT(2)-peak{2}(s,t)].^2);
        end
    end
    X=double(sum(v,5)<=rtmax);
    for t=1:NT(4)
        for s=1:NT(3)
            X(round(peak{1}(s,t)),round(peak{2}(s,t)),s,t)=1;
        end
    end    
    X=ifftshift(ifftshift(X,1),2); %Centered spectrum
    fL(:,:,:,:,1)=ifft2(X.*fL(:,:,:,:,1));
elseif strcmp(phaseEstMethod,'MP-FT')
    fL1=single(zeros([NT.*[1 1 length(mp) 1] 2]));
    dp=0;%Spectral peak is taken as the prescribed pattern frequency
    for npeak=1:length(mp) 
    rRange=0.4/abs(mp(npeak));    
    or=NT(1:2)/2+1;%Center of the k-space in pixels
    rk=NT(1:2).*sT*abs(mp(npeak))/kT(1);%Radius of the spectral peak
    pk=bsxfun(@times,rk,ak)*sign(mp(npeak));%k-space position of the spectral peak
    ik=bsxfun(@plus,or,pk);%k-space position of the spectral peak in pixels 
    
    rpeak=[((1-rRange)*rk(1)).^2 ((1+rRange)*rk(1)).^2];%Range for ridge searching in the radial direction
    Xpeak=peakROI(or,NT,orT+pi*(mp(npeak)<0)+mo(npeak),rpeak,aRange);%ROI to search the peak
    fL(:,:,:,:,1)=fft2(ImT);       
    [peak,rtmax]=ridgeDetector(fL,Xpeak,dp,ik,rk,BW/abs(mp(npeak)));
    v=single(zeros([NT 2]));   
    for t=1:NT(4)
        for s=1:NT(3)
            [v(:,:,s,t,1),v(:,:,s,t,2)]=ndgrid([1-peak{1}(s,t):NT(1)-peak{1}(s,t)].^2,[1-peak{2}(s,t):NT(2)-peak{2}(s,t)].^2);
        end
    end
    X=double(sum(v,5)<=rtmax);
    for t=1:NT(4)
        for s=1:NT(3)
            X(round(peak{1}(s,t)),round(peak{2}(s,t)),s,t)=1;
        end
    end    
    X=ifftshift(ifftshift(X,1),2); %Centered spectrum
    fL1(:,:,(npeak-1)*NT(3)+1:npeak*NT(3),:,1)=ifft2(X.*fL(:,:,:,:,1));
    fL(:,:,:,:,1)=ifft2(X.*fL(:,:,:,:,1));
    
    end
    clear fL; fL=fL1;
elseif strcmp(phaseEstMethod,'SFT')
    dp=0;%Spectral peak is taken as the prescribed pattern frequency
    or=NT(1:2)/2+1;%Center of the k-space in pixels
    rk=NT(1).*sT(1)/kT(1);%Radius of the spectral peak
    pk=[rk 0];%k-space position of the spectral peak    
    ik=bsxfun(@plus,or,pk);%k-space position of the spectral peak in pixels
    ik=repmat(ik,[NT(3) 1]);
    
    rpeak=[((1-rRange)*rk(1)).^2 ((1+rRange)*rk(1)).^2];%Range for ridge searching in the radial direction           
    Xpeak=peakROI(or,NT,pi/2*ones(1,NT(3)),rpeak,aRange,0);%ROI to search the peak   
    
    [kGrid,kkGrid,rGrid,rkGrid]=generateGrids([NT(1:2) 1]);
    TGT=single(zeros(1,1,1,1,1,6));   
    orTRot=orT-pi/2;
    ind_flip=(orTRot>=pi/2);    
    orTRot(orTRot>=pi/2)=orTRot(orTRot>=pi/2)-pi;
    for s=1:NT(3)
        TGT(4)=-orTRot(s);
        et=precomputationsSincRigidTransform(kGrid,[],rkGrid,TGT,1,0);
        if ~ind_flip(s)
            ImT(:,:,s)=ImT(end:-1:1,end:-1:1,s);
        end
        ImT(:,:,s)=sincRigidTransform(ImT(:,:,s),et,1,0);
    end      
    fL(:,:,:,:,1)=fftGPU(ImT,1,0);
    [peak,rtmax]=ridgeDetector(fL,Xpeak,dp,ik,rk,BW,1);
    v=single(zeros([NT 2]));    
    for t=1:NT(4)
        for s=1:NT(3)
            [v(:,:,s,t,1),v(:,:,s,t,2)]=ndgrid([1-peak{1}(s,t):NT(1)-peak{1}(s,t)].^2,[1-peak{2}(s,t):NT(2)-peak{2}(s,t)].^2);
        end
    end
    X=double(sum(v(:,:,:,:,1),5)<=rtmax);
    for t=1:NT(4)
        for s=1:NT(3)
            X(round(peak{1}(s,t)),round(peak{2}(s,t)),s,t)=1;
        end
    end    
    %figure(1);for s=1:44;subtightplot(4,11,s,[0 0]);imshow(X(:,:,s,1,1),[]);end
    X=ifftshift(X,1); %Centered spectrum
    %pause
    fL(:,:,:,:,1)=ifftGPU(X.*fL(:,:,:,:,1),1,0);   
    for s=1:NT(3)
        TGT(4)=orTRot(s);        
        et=precomputationsSincRigidTransform(kGrid,[],rkGrid,TGT,1,0);        
        fL(:,:,s,:,1)=sincRigidTransform(fL(:,:,s,:,1),et,1,0);
        if ~ind_flip(s)
            fL(:,:,s,:,1)=fL(end:-1:1,end:-1:1,s,:,1);
        end
    end
elseif strcmp(phaseEstMethod,'WFT')        
    dp=1;%Spectral peak is estimated on a data basis
    %dp=0;
    alpha=2.5;
    Nori=NW/2;
    Namp=NW-1;
    mDi=cell(1,2);mImO=cell(1,2);mImF=cell(1,2);MWi=single(zeros(1,2));mIn=cell(1,2);mInOr=cell(1,2);q=cell(1,2);
    for s=1:2        
        mDi{s}=lims(1,s)-NW/2+1:lims(2,s)+NW/2;
        mImO{s}=mDi{s}-Nori;
        mImF{s}=mImO{s}+Namp;
        MWi(s)=length(mDi{s});
        mIn{s}=lims(1,s):lims(2,s);
        mInOr{s}=mIn{s}-lims(1,s);
        q{s}=NW-(1:NW)+1;
    end
    
    or=[NW NW]/2+1;%Center of the k-space in pixels
    rk=[NW NW].*sT/kT(1);%Radius of the spectral peak
    pk=bsxfun(@times,rk,ak);%k-space position of the spectral peak
    ik=bsxfun(@plus,or,pk);%k-space position of the spectral peak in pixels     
    rpeak=[((1-rRange)*rk(1)).^2 ((1+rRange)*rk(1)).^2];%Range for ridge searching in the radial direction
    
    W=single(gausswin(NW,alpha));
    W=W*W';    
    fL=single(zeros([NT 2]));
      
    NA=[NW NW NT(3)];
    %v=single(zeros([NW NW MWi 2]));
    Xpeak=peakROI(or,NA,orT,rpeak,aRange);
    TF=single(zeros([NW NW MWi]));
    TFaux=single(zeros([NW NT(2) MWi(1) 1]));
    fLaux=cell(1,2);
    fLaux{1}=single(zeros([NT(1) MWi(2) 1 NW]));
    fLaux{2}=single(zeros(NT(1:2)));    
    if gpu
        fLaux{1}=gpuArray(fLaux{1});fLaux{2}=gpuArray(fLaux{2});Xpeak=gpuArray(Xpeak);TFaux=gpuArray(TFaux);TF=gpuArray(TF);ImT=gpuArray(ImT);W=gpuArray(W);fL=gpuArray(fL);
    end        
    
    for p=1:NT(4)
        for s=1:NT(3)
            for m=1:MWi(1)
                TFaux(:,:,m)=ImT(mImO{1}(m):mImF{1}(m),:,s,p);
            end
            for n=1:MWi(2)
                TF(:,:,:,n)=TFaux(:,mImO{2}(n):mImF{2}(n),:);
            end
            TFS=bsxfun(@times,TF,W);            
            for m=1:2
                TFS=fftGPU(TFS,m,gpu);
            end
            %Takes certain time, maybe possible to accelerate
            [peak rtmax]=ridgeDetector(TFS,Xpeak(:,:,s),dp,ik(s,:),rk,BW);            
            if ~dp
                for m=1:2
                    peak{m}=repmat(peak{m},[size(TFS,3) 1]);
                end                
            end                       
            grip=cell(1,2);
            for m=1:2
                peak{m}=permute(peak{m},[3 4 1 2]);
                grip{m}=single(1:NW);
                if gpu
                    grip{m}=gpuArray(grip{m});
                end
            end
            grip{1}=grip{1}';
            for m=1:2
                grip{m}=bsxfun(@minus,grip{m},peak{m}).^2;
            end
            X=single(bsxfun(@plus,grip{1},grip{2})<=rtmax);  
            X=ifftshift(ifftshift(X,1),2);
            TFS=TFS.*X;        
            for m=1:2
                TFS=ifftGPU(TFS,m,gpu);
            end
            TFS=bsxfun(@times,W,TFS);                 
            TFS=shiftdim(TFS,2);   
            
            fLaux{1}(:)=0;
            for m=1:NW
                fLaux{1}(mIn{1},:,1,:)=fLaux{1}(mIn{1},:,1,:)+TFS(mInOr{1}+m,:,q{1}(m),:);
            end
            fLaux{2}(:)=0;
            for n=1:NW                
                fLaux{2}(:,mIn{2})=fLaux{2}(:,mIn{2})+fLaux{1}(:,mInOr{2}+n,1,q{2}(n));%*W(q(s),r(t),o);
            end
            fL(:,:,s,p,1)=fLaux{2};
        end                
    end
    if gpu
        fL=gather(fL);
    end
elseif strcmp(phaseEstMethod,'AWFT')        
    dp=1;%Spectral peak is estimated on a data basis
    %dp=0;
    Nori=NW/2;
    Namp=NW-1;
    mDi=cell(1,2);mImO=cell(1,2);mImF=cell(1,2);MWi=single(zeros(1,2));mIn=cell(1,2);mInOr=cell(1,2);q=cell(1,2);
    for s=1:2        
        mDi{s}=lims(1,s)-NW/2+1:lims(2,s)+NW/2;
        mImO{s}=mDi{s}-Nori;
        mImF{s}=mImO{s}+Namp;
        MWi(s)=length(mDi{s});
        mIn{s}=lims(1,s):lims(2,s);
        mInOr{s}=mIn{s}-lims(1,s);
        q{s}=NW-(1:NW)+1;
    end
    
    or=[NW NW]/2+1;%Center of the k-space in pixels
    rk=[NW NW].*sT/kT(1);%Radius of the spectral peak
    pk=bsxfun(@times,rk,ak);%k-space position of the spectral peak
    ik=bsxfun(@plus,or,pk);%k-space position of the spectral peak in pixels     
    rpeak=[((1-rRange)*rk(1)).^2 ((1+rRange)*rk(1)).^2];%Range for ridge searching in the radial direction
    
    fL=single(zeros([NT 2]));
      
    NA=[NW NW NT(3)];
    %v=single(zeros([NW NW MWi 2]));
    Xpeak=peakROI(or,NA,orT,rpeak,aRange);
    TF=single(zeros([NW NW MWi]));
    TFaux=single(zeros([NW NT(2) MWi(1) 1]));
    fLaux=cell(1,2);
    fLaux{1}=single(zeros([NT(1) MWi(2) 1 NW]));
    fLaux{2}=single(zeros(NT(1:2)));    
    if gpu
        fLaux{1}=gpuArray(fLaux{1});fLaux{2}=gpuArray(fLaux{2});Xpeak=gpuArray(Xpeak);TFaux=gpuArray(TFaux);TF=gpuArray(TF);ImT=gpuArray(ImT);W=gpuArray(W);fL=gpuArray(fL);
    end        
    
    for p=1:NT(4)
        for s=1:NT(3)
            W=window_or(NW,kT(s)*1.1,90-orT(s)*180/pi); %alpha=(NW-1)/(2*sigma)
            for m=1:MWi(1)
                TFaux(:,:,m)=ImT(mImO{1}(m):mImF{1}(m),:,s,p);
            end
            for n=1:MWi(2)
                TF(:,:,:,n)=TFaux(:,mImO{2}(n):mImF{2}(n),:);
            end
            TFS=bsxfun(@times,TF,W);            
            for m=1:2
                TFS=fftGPU(TFS,m,gpu);
            end
            %Takes certain time, maybe possible to accelerate
            [peak rtmax]=ridgeDetector(TFS,Xpeak(:,:,s),dp,ik(s,:),rk,BW);            
            if ~dp
                for m=1:2
                    peak{m}=repmat(peak{m},[size(TFS,3) 1]);
                end                
            end                       
            grip=cell(1,2);
            for m=1:2
                peak{m}=permute(peak{m},[3 4 1 2]);
                grip{m}=single(1:NW);
                if gpu
                    grip{m}=gpuArray(grip{m});
                end
            end
            grip{1}=grip{1}';
            for m=1:2
                grip{m}=bsxfun(@minus,grip{m},peak{m}).^2;
            end
            %X2=single(bsxfun(@plus,grip{1},grip{2})<=rtmax);  
            X=single(elliptic_filter( NW, peak, orT(s), rtmax,rk(1)^2, 0.9 )); %eccentricity
            %imshow([X2(:,:,10,10).*log(abs(fftshift(fftshift(TFS(:,:,10,10),1),2))) X(:,:,10,10).*log(abs(fftshift(fftshift(TFS(:,:,10,10),1),2)))],[]), pause
            %imshow([X2(:,:,10,10) X(:,:,10,10)],[]), pause
            X=ifftshift(ifftshift(X,1),2);
            TFS=TFS.*X;        
            for m=1:2
                TFS=ifftGPU(TFS,m,gpu);
            end
            TFS=bsxfun(@times,W,TFS);                 
            TFS=shiftdim(TFS,2);   
            
            fLaux{1}(:)=0;
            for m=1:NW
                fLaux{1}(mIn{1},:,1,:)=fLaux{1}(mIn{1},:,1,:)+TFS(mInOr{1}+m,:,q{1}(m),:);
            end
            fLaux{2}(:)=0;
            for n=1:NW                
                fLaux{2}(:,mIn{2})=fLaux{2}(:,mIn{2})+fLaux{1}(:,mInOr{2}+n,1,q{2}(n));%*W(q(s),r(t),o);
            end
            fL(:,:,s,p,1)=fLaux{2};
        end                
    end
    if gpu
        fL=gather(fL);
    end
elseif strcmp(phaseEstMethod,'MP-WFT') 
    orT1=orT; fL1=single(zeros([NT.*[1 1 length(mp) 1] 2]));
    dp=1;%Spectral peak is estimated on a data basis
    for npeak=1:length(mp)     
     
    rRange=0.4/abs(mp(npeak));
            
    Nori=NW/2;
    Namp=NW-1;
    mDi=cell(1,2);mImO=cell(1,2);mImF=cell(1,2);MWi=single(zeros(1,2));mIn=cell(1,2);mInOr=cell(1,2);q=cell(1,2);
    for s=1:2        
        mDi{s}=lims(1,s)-NW/2+1:lims(2,s)+NW/2;
        mImO{s}=mDi{s}-Nori;
        mImF{s}=mImO{s}+Namp;
        MWi(s)=length(mDi{s});
        mIn{s}=lims(1,s):lims(2,s);
        mInOr{s}=mIn{s}-lims(1,s);
        q{s}=NW-(1:NW)+1;
    end
    
    or=[NW NW]/2+1;%Center of the k-space in pixels
    rk=[NW NW].*sT*abs(mp(npeak))/kT(1);%Radius of the spectral peak
    pk=bsxfun(@times,rk,ak)*sign(mp(npeak));%k-space position of the spectral peak
    ik=bsxfun(@plus,or,pk);%k-space position of the spectral peak in pixels     
    rpeak=[((1-rRange)*rk(1)).^2 ((1+rRange)*rk(1)).^2];%Range for ridge searching in the radial direction
    
    fL=single(zeros([NT 2]));
      
    orT=orT1 +pi*(mp(npeak)<0)+mo(npeak);
    NA=[NW NW NT(3)];
    %v=single(zeros([NW NW MWi 2]));
    Xpeak=peakROI(or,NA,orT,rpeak,aRange);
    TF=single(zeros([NW NW MWi]));
    TFaux=single(zeros([NW NT(2) MWi(1) 1]));
    fLaux=cell(1,2);
    fLaux{1}=single(zeros([NT(1) MWi(2) 1 NW]));
    fLaux{2}=single(zeros(NT(1:2)));    
    if gpu
        fLaux{1}=gpuArray(fLaux{1});fLaux{2}=gpuArray(fLaux{2});Xpeak=gpuArray(Xpeak);TFaux=gpuArray(TFaux);TF=gpuArray(TF);ImT=gpuArray(ImT);W=gpuArray(W);fL=gpuArray(fL);
    end        
    
    W=single(gausswin(NW,2)); W=W*W';
    for p=1:NT(4)
        for s=1:NT(3)
            %W=window_or(NW,kT(s)*1.1,90-orT1(s)*180/pi); %alpha=(NW-1)/(2*sigma)
            for m=1:MWi(1)
                TFaux(:,:,m)=ImT(mImO{1}(m):mImF{1}(m),:,s,p);
            end
            for n=1:MWi(2)
                TF(:,:,:,n)=TFaux(:,mImO{2}(n):mImF{2}(n),:);
            end
            TFS=bsxfun(@times,TF,W);            
            for m=1:2
                TFS=fftGPU(TFS,m,gpu);
            end
            %Takes certain time, maybe possible to accelerate
            [peak rtmax]=ridgeDetector(TFS,Xpeak(:,:,s),dp,ik(s,:),rk,(1.3^-(abs(mp(npeak))-1))*BW/abs(mp(npeak)));   
            if ~dp
                for m=1:2
                    peak{m}=repmat(peak{m},[size(TFS,3) 1]);
                end                
            end                       
            grip=cell(1,2);
            for m=1:2
                peak{m}=permute(peak{m},[3 4 1 2]);
                grip{m}=single(1:NW);
                if gpu
                    grip{m}=gpuArray(grip{m});
                end
            end
            grip{1}=grip{1}';
            for m=1:2
                grip{m}=bsxfun(@minus,grip{m},peak{m}).^2;
            end
            X=single(bsxfun(@plus,grip{1},grip{2})<=rtmax);  
            %X=single(elliptic_filter( NW, peak, orT(s), rtmax,rk(1)^2, 0.9 )); %eccentricity
            %imshow( imresize([X(:,:,10,10).*log(abs(fftshift(fftshift(TFS(:,:,10,10),1),2))) log(abs(fftshift(fftshift(TFS(:,:,10,10),1),2)))],4),[]), pause
            %imshow([X2(:,:,10,10) X(:,:,10,10)],[]), pause
            X=ifftshift(ifftshift(X,1),2);
            TFS=TFS.*X;        
            for m=1:2
                TFS=ifftGPU(TFS,m,gpu); 
            end
            TFS=bsxfun(@times,W,TFS);                 
            TFS=shiftdim(TFS,2);   
            
            fLaux{1}(:)=0;
            for m=1:NW
                fLaux{1}(mIn{1},:,1,:)=fLaux{1}(mIn{1},:,1,:)+TFS(mInOr{1}+m,:,q{1}(m),:);
            end
            fLaux{2}(:)=0;
            for n=1:NW                
                fLaux{2}(:,mIn{2})=fLaux{2}(:,mIn{2})+fLaux{1}(:,mInOr{2}+n,1,q{2}(n));%*W(q(s),r(t),o);
            end
            fL(:,:,s,p,1)=fLaux{2};
            fL1(:,:,(npeak-1)*NT(3)+s,p,1)=fLaux{2};
        end                
    end
    end
    fL=fL1;
    if gpu
        fL=gather(fL); 
    end
else
    dp=1;%Spectral peak is estimated on a data basis
    %dp=0;
    alpha=2;
    Nori=NW/2;
    Namp=NW-1;
    mDi=cell(1,1);mImO=cell(1,1);mImF=cell(1,1);MWi=single(zeros(1,1));mIn=cell(1,1);mInOr=cell(1,1);q=cell(1,1);
    for s=1:2
        lims(1,s)=lims(1,s)-10;    
        lims(2,s)=lims(2,s)+10; 
    end
    for s=1:1
        mDi{s}=lims(1,s)-NW/2+1:lims(2,s)+NW/2;
        mImO{s}=mDi{s}-Nori;
        mImF{s}=mImO{s}+Namp;
        MWi(s)=length(mDi{s});
        mIn{s}=lims(1,s):lims(2,s);
        mInOr{s}=mIn{s}-lims(1,s);
        q{s}=NW-(1:NW)+1;
    end

    or=[NW NW]/2+1;%Center of the k-space in pixels
    rk=NW.*sT(1)/kT(1);%Radius of the spectral peak
    pk=[rk 0];%k-space position of the spectral peak
    ik=bsxfun(@plus,or,pk);%k-space position of the spectral peak in pixels     
    ik=repmat(ik,[NT(3) 1]);
    rpeak=[((1-rRange)*rk(1)).^2 ((1+rRange)*rk(1)).^2];%Range for ridge searching in the radial direction
    
    W=single(gausswin(NW,alpha));   
    fL=single(zeros([NT 2]));
      
    NA=[NW NT(3)];
    %v=single(zeros([NW NW MWi 2]));
    Xpeak=peakROI(or,NA,pi/2*ones(1,NT(3)),rpeak,aRange,1);
    
    TF=single(zeros([NW NT(2) MWi]));      
    fLaux=single(zeros(NT(1:2)));   
    
    [kGrid,kkGrid,rGrid,rkGrid]=generateGrids([NT(1:2) 1]);
    TGT=single(zeros(1,1,1,1,1,6));   
    orTRot=orT-pi/2;
    ind_flip=(orTRot>=pi/2);    
    orTRot(orTRot>=pi/2)=orTRot(orTRot>=pi/2)-pi;
    for s=1:NT(3)
        TGT(4)=-orTRot(s);
        et=precomputationsSincRigidTransform(kGrid,[],rkGrid,TGT,1,0);
        if ~ind_flip(s)
            ImT(:,:,s)=ImT(end:-1:1,end:-1:1,s);
        end
        ImT(:,:,s)=sincRigidTransform(ImT(:,:,s),et,1,0);
    end 
   
    if gpu
        fLaux=gpuArray(fLaux);Xpeak=gpuArray(Xpeak);TF=gpuArray(TF);ImT=gpuArray(ImT);W=gpuArray(W);fL=gpuArray(fL);
    end        
    
    for p=1:NT(4)
        for s=1:NT(3)
            for m=1:MWi(1)
                TF(:,:,m)=ImT(mImO{1}(m):mImF{1}(m),:,s,p);
            end
            TFS=bsxfun(@times,TF,W);            
            for m=1:1
                TFS=fftGPU(TFS,m,gpu);
            end
            %Takes certain time, maybe possible to accelerate
            [peak rtmax]=ridgeDetector(TFS,Xpeak(:,s),dp,ik(s,:),rk,BW,1);            
            %size(peak)
            %peak
            %size(rtmax)
            %pause
            if ~dp            
                peak=repmat(peak{1},[size(TFS,2) size(TFS,3)]);                
            end                       
            grip=cell(1,1);
            for m=1:1
                peak=permute(peak,[3 1 2]);
                grip{m}=single(1:NW);
                if gpu
                    grip{m}=gpuArray(grip{m});
                end
            end
            grip{1}=grip{1}';
            for m=1:1
                grip{m}=bsxfun(@minus,grip{m},peak).^2;
            end
            %size(grip{1})
            %size(rtmax)
            %pause
            X=single(grip{1}<=rtmax);  
            X=ifftshift(X,1);
            TFS=TFS.*X;        
            for m=1:1
                TFS=ifftGPU(TFS,m,gpu);
            end
            TFS=bsxfun(@times,W,TFS);                 
            TFS=permute(TFS,[3 2 1]);   
            fLaux(:)=0;
            for m=1:NW
                fLaux(mIn{1},:)=fLaux(mIn{1},:)+TFS(mInOr{1}+m,:,q{1}(m));
            end
            fL(:,:,s,p,1)=fLaux;
        end                
    end
    if gpu
        fL=gather(fL);
    end
    for s=1:NT(3)
        TGT(4)=orTRot(s);        
        et=precomputationsSincRigidTransform(kGrid,[],rkGrid,TGT,1,0);        
        fL(:,:,s,:,1)=sincRigidTransform(fL(:,:,s,:,1),et,1,0);
        if ~ind_flip(s)
            fL(:,:,s,:,1)=fL(end:-1:1,end:-1:1,s,:,1);
        end
    end
    for s=1:2
        lims(1,s)=lims(1,s)+10;    
        lims(2,s)=lims(2,s)-10; 
    end
end

limsb=lims;

fL=angle(fL);
fL=fL(lims(1,1):lims(2,1),lims(1,2):lims(2,2),:,:,1);
Wr=fL(:,:,:,:,1)+pi;
Wr(Wr>pi)=Wr(Wr>pi)-2*pi;
fL(:,:,:,:,2)=Wr;

telapsed=toc(tstart);
fprintf('Local phase estimation time: %f\n',telapsed)
