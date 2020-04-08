function lims=limitROI(Adat,widen)

%LIMITROI   Obtains a myocardial mask
%   LIMITANTS=LIMITROI(ADAT,WIDEN) obtains a ROI around the myocardial mask
%   ADAT are the input masks
%   WIDEN is the security margin to extend the mask, defaults to 0
%   It returns:
% %   LIMS, the limitant region
%

if ~exist('widen','var')
    widen=0;
end

N=size(Adat);
lims(1,:)=ones(1,length(N));
lims(2,:)=N(1:length(N));

a(1)=1;
a(2)=-1;
for l=1:2    
    for n=1:length(N)
        Bdat=shiftdim(Adat,n-1);
        while isempty(find(Bdat(lims(l,n),:,:,:)~=0,1))
            lims(l,n)=lims(l,n)+a(l);
        end
    end
end
lims(1,1:2)=lims(1,1:2)-widen;
lims(2,1:2)=lims(2,1:2)+widen;
