function [ImT,kT,orT,Tred,raC,MTs,limsOu,ErrMGT,EccMGT,ErcMGT,lambdas]=generatePhantomHARP(NOR)

N=[192 192];
cent(1:2)=N/2-0.5;
[X(:,:,1),X(:,:,2)] =ndgrid(1:N(1),1:N(2));
X=bsxfun(@minus,X,permute(cent,[1 3 2]));
raC=bsxfun(@times,X,sum(X.^2,3).^(-0.5));%Orientations
Rad=sum(X.^2,3).^0.5;
Ri=128*7/32;
Ro=128*5/16;
MTd=single(Rad>=Ri & Rad<Ro);

S=NOR;
orT=(1:S)*pi/S;
S=length(orT);
for s=1:S
    ak(s,:)=[-sin(orT(s)-pi/4) cos(orT(s)-pi/4)];%Orientation of the spectral peak
    kT(s)=7.15; %7.15
end

useMask=0;
%multiple peak simulation
%proposed by Osman {3.64, 1.00, -0.89, 0.72, -0.52, 0.32, -0.15}
for s=1:S    
    if useMask
        ImTOr(:,:,s)=MTd.*(1+0.5*cos(2*pi*(X(:,:,1)*ak(s,1)+X(:,:,2)*ak(s,2))/kT(s)));
        %ImTOr(:,:,s)=MTd.*(3.64+1*cos(2*pi*(X(:,:,1)*ak(s,1)+X(:,:,2)*ak(s,2))/kT(s)) - ...
        %            0.89*cos(4*pi*(X(:,:,1)*ak(s,1)+X(:,:,2)*ak(s,2))/kT(s)) + ...
        %            0.72*cos(6*pi*(X(:,:,1)*ak(s,1)+X(:,:,2)*ak(s,2))/kT(s)));
    else
        ImTOr(:,:,s)=1+0.5*cos(2*pi*(X(:,:,1)*ak(s,1)+X(:,:,2)*ak(s,2))/kT(s));
        %ImTOr(:,:,s)=3.64+1*cos(2*pi*(X(:,:,1)*ak(s,1)+X(:,:,2)*ak(s,2))/kT(s)) - ...
        %    0.89*cos(4*pi*(X(:,:,1)*ak(s,1)+X(:,:,2)*ak(s,2))/kT(s)) + ...
        %    0.72*cos(6*pi*(X(:,:,1)*ak(s,1)+X(:,:,2)*ak(s,2))/kT(s));
       
       
    end
    %figure(s)
    %imshow(ImT(:,:,s),[])
end

lambdas=0:0.05:0.8;
L=length(lambdas);

ImT=cell(1,L);
for l=1:L
    ImT{l}=zeros([N S]);
    lambda=lambdas(l);
    ri=Ri*(1-lambda);
    radCov=Ri^2+(Rad.^2-ri^2);
    %rad=ri^2+(Rad.^2-Ri^2);
    radCov(radCov<0)=0;
    %rad(rad<0)=0;
    radCov=sqrt(radCov);
    %rad=sqrt(rad);
    x=bsxfun(@times,radCov,raC);
    %xx=bsxfun(@times,rad,raC);

    MTs{l}=interpn(X(:,:,1),X(:,:,2),MTd,x(:,:,1),x(:,:,2),'nearest',0);%single(rad>=Ri & rad<Ro);
    for s=1:S
        ImT{l}(:,:,s)=interpn(X(:,:,1),X(:,:,2),ImTOr(:,:,s),x(:,:,1),x(:,:,2),'linear',0);
    end
    Tred=single(zeros(2,S));

    FGT(:,:,1,1,1)=x([2:end 1],:,1)-x([end 1:end-1],:,1);
    FGT(:,:,1,2,1)=x([2:end 1],:,2)-x([end 1:end-1],:,2);
    FGT(:,:,1,1,2)=x(:,[2:end 1],1)-x(:,[end 1:end-1],1);
    FGT(:,:,1,2,2)=x(:,[2:end 1],2)-x(:,[end 1:end-1],2);
    FGT=FGT/2;

    widen=0;
    lims=limitROI(ones(size(MTs{l})),widen);
    %[ErrSGT,EccSGT,ErcSGT]=ejectionStrain(fGT(:,:,1,:,:),MTs,lims,raC,1);%Material case%Spatial would be with 1 in the last input
    [ErrMGT{l},EccMGT{l},ErcMGT{l}]=ejectionStrain(FGT(:,:,1,:,:),MTs{l},lims,raC,1);%Material case%Spatial would be with 1 in the last input
    widen=5;
    limsOu{l}=limitROI(MTs{l},widen);
end
