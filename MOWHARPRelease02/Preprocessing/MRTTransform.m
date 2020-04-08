function To=MRTTransform(Ti,sT,T)

%MRTTRANSFORM   Transforms MR-T images
%   TO=MRTTRANSFORM(TI,ST,T) transforms a set of differently oriented MR-T
%   masks according to a precomputed set of transforms that align them to 
%   the MR-C reference
%   TI are the input masks
%   ST is the MR-T resolution
%   T the transform parameters (translations)
%   It returns:
%   TO, the transformed masks
%

NT=size(Ti);
NT(end+1:5)=1;

tstart=tic;
[xT(:,:,1),xT(:,:,2)]=ndgrid(1:sT(1):NT(1)*sT(1),1:sT(2):NT(2)*sT(2));
xTact=bsxfun(@plus,xT,permute(T,[3 4 1 2]));

To=single(zeros(size(Ti)));

for u=1:NT(5)
    for t=1:NT(4)  
        for s=1:NT(3)
            To(:,:,s,t,u)=interpn(xT(:,:,1),xT(:,:,2),Ti(:,:,s,t,u),xTact(:,:,1,s),xTact(:,:,2,s),'nearest');
        end
    end
end

To(isnan(To))=0; 
telapsed=toc(tstart);
fprintf('Time transforming MR-T: %f\n',telapsed)
