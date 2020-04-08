function [x,xT]=EncodingElastic(rec,x,T,forw)
%encoding MFDST forward (1), backward (0)

%ac activate fast resampling
%if nargin<5;ac=0;end

gpu=isa(x,'gpuArray');
if gpu==1;gpu=2;end

if (isfield(rec,'ParX')==0) || (rec.ParX.flagcoil==0)
    if forw
        x=applyTransform(x,T,rec.NDimen); %T
        if nargout>=2;xT=x;end
        x=bsxfun(@times,x,rec.S); %S
        if gpu==0
            x=resampling(x,rec.NY); %D downsampling
            for m=1:rec.NDimen;x=fftGPU(x,m,gpu);end    %F
        else
            for m=1:rec.NDimen;x=fftGPU(x,m,gpu,rec.F{m});end
        end
        x=sum(bsxfun(@times,x,rec.A),5); %M apply shot mask
    else
        x=bsxfun(@times,x,rec.A); %M apply shot mask
        if gpu==0
            for m=1:rec.NDimen;x=ifftGPU(x,m,gpu);end %F
            x=resampling(x,rec.NX); %D upsampling
        else
            for m=1:rec.NDimen;x=ifftGPU(x,m,gpu,rec.F{m}');end 
            %/(rec.NY(m)/(rec.NX(m)/rec.NY(m)));end
            %HAVEN'T CHECKED BUT THIS COULD PROVIDE QUICKER RECONSTRUCTIONS AS WELL, 
            %WHICH MAY BE COMPUTATIONALLY APPEALING IF MOTION ESTIMATION IS PERFORMED 
            %AT, LET'S SAY, HALF THE RESOLUTION AND THERE IS A FINAL RECONSTRUCTIONS 
            %STEP AT FULL RESOLUTION
        end
        x=sum(bsxfun(@times,x,conj(rec.S)),4); %Sconj
        x=sum(applyTransform(x,T,rec.NDimen),5); % Tinv
    end
    
elseif rec.ParX.flagcoil==1 %compute x coil by coil
    if forw
        x2=zeros([rec.NY rec.Nco],'like',x);
        %if nargout>=2; xT=applyTransform(x,T,rec.NDimen);end
        for co=1:rec.Nco
            for sh=1:size(T,5)
                x1=applyTransform(x,T(:,:,:,:,sh),rec.NDimen);
                if nargout>=2; if co==1, xT(:,:,:,:,sh)=x1;end;end
                x1=x1.*rec.S(:,:,:,co);
                if gpu==0 %ojo posible bug en estas 2 lineas
                    x1=resampling(x1,rec.NY); %D downsampling
                    for m=1:rec.NDimen;x1=fftGPU(x1,m,gpu);end    %F
                else
                    for m=1:rec.NDimen;x1=fftGPU(x1,m,gpu,rec.F{m});end
                end
                x2(:,:,:,co)=x2(:,:,:,co)+bsxfun(@times,x1,rec.A(:,:,:,1,sh)); %M apply shot mask
            end
        end
        x=x2;
    else
        x=bsxfun(@times,x,rec.A); %M apply shot mask
        x2=zeros(rec.NX ,'like',x);
        for co=1:rec.Nco
            for sh=1:size(T,5)
                x1=x(:,:,:,co,sh);
                if gpu==0
                    for m=1:rec.NDimen;x1=ifftGPU(x1,m,gpu);end %F
                    x1=resampling(x1,rec.NX); %D upsampling
                else
                    for m=1:rec.NDimen;x1=ifftGPU(x1,m,gpu,rec.F{m}');end
                end
                x2=x2+applyTransform(x1.*conj(rec.S(:,:,:,co)),T(:,:,:,:,sh),rec.NDimen); %Sconj Tinv
            end
        end
        x=x2;
    end
end

 
