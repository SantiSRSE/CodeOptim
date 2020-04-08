function [indic,rtmax]=ridgeDetector(fL,Xpeak,dp,ik,rk,BW,st)

%RIDGEDETECTOR   Detects the spectral peak 
%   [INDIC,RTMIN]=RIDGEDETECTOR(FL,XPEAK,DP,IK,RK,BW) determines the
%   spectral peak for local phase estimation
%   FL is the spectrum/windowed spectrum
%   XPEAK is the ROI to search the peak
%   DP is a flag that indicates whether the spectral peak is taken as the
%   prescribed pattern frequency (0) or is estimated on a data basis (1)
%   IK is the k-space position of the spectral peak in pixels
%   RK is the radius of the spectral peak
%   BW is the normalized bandwidth
%   ST indicates whether the steered FT is used
%   It returns:
%   INDIC, the peak location
%   RTMAX, the square of the radius of the BP filter
%

if ~exist('st','var')
    st=0;
end

NT=size(fL);
NT(end+1:4)=1;
indic=cell(1,2);
rtmax=((BW*rk(1))^2)*ones(NT);
if ~dp 
    for m=1:2
        indic{m}=repmat(ik(:,m),[1 NT(4)]);
    end
else
    if ~st
        Mod=fftshift(fftshift(abs(fL),1),2);
    else
        Mod=fftshift(abs(fL),1);
    end
    ModX=bsxfun(@times,Mod,Xpeak);       
    if ~st
        ModX=reshape(ModX,[NT(1)*NT(2) NT(3:4)]);
        [maxim ind]=max(ModX,[],1);    
        ind=shiftdim(ind,1);
        [indic{1},indic{2}] = ind2sub(NT(1:2),ind);    
    else
        [maxim ind]=max(ModX,[],1);
        indic=shiftdim(ind,1);
    end    
end
