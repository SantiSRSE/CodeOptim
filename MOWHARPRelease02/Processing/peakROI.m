function Xpeak=peakROI(or,NT,orT,rPeak,aPeak,st)

%PEAKROI   Establishes the ROI where to search for the spectral peak 
%   XPEAK=PEAKROI(OR,NT,ORT,RPEAK,APEAK) 
%   OR is the center of the k-space in pixels
%   NT is the size of the MR-T images
%   ORT is the orientation of the tags (in radians)
%   RPEAK is the range for ridge searching in the radial direction
%   APEAK is the range for ridge searching in the angular direction
%   ST indicates whether the steered FT is used
%   It returns:
%   XPEAK, the ROI to search the peak
%

if ~exist('st','var')
    st=0;
end
NT(end+1:4)=1;

if ~st
    [vP(:,:,1),vP(:,:,2)]=ndgrid(1-or(1):NT(1)-or(1),1-or(2):NT(2)-or(2));
    vR=sum(vP.^2,3);
    Xrad=repmat(single(vR<=rPeak(2) & vR>=rPeak(1)),[1 1 NT(3)]);
    VT=single(zeros(NT(1:3)));
    for s=1:NT(3)
        %VT(:,:,s)=atan2(vP(:,:,2),-vP(:,:,1))+pi-orT(s);
        VT(:,:,s)=atan2(-vP(:,:,1),vP(:,:,2))-orT(s);
        %VT(:,:,s)=atan2(-vP(:,:,2),vP(:,:,1))-orT(s);
    end
    VT(VT>pi)=VT(VT>pi)-2*pi;
    VT(VT<-pi)=2*pi+VT(VT<-pi);
    Xpeak=single(abs(VT)<aPeak);
    Xpeak=Xpeak.*Xrad;    
    Xpeak=repmat(Xpeak,[1 1 1 NT(4)]);
else
    vP=ndgrid(1-or(1):NT(1)-or(1));
    vR=vP.^2;
    Xrad=repmat(single(vR<=rPeak(2) & vR>=rPeak(1)),[1 NT(3)]);
    Xpeak=Xrad.*single(vP>0);    
    Xpeak=repmat(Xpeak,[1 NT(2) 1 NT(4)]);
end

