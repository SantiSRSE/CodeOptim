function [T,Tred]=MRTRegistration(ImT,ImC,sT,sC,gpu)

%MRTREGISTRATION   Registers MR-T images for different orientations 
%   [T,TRED]=MRTREGISTRATION(IMT,IMC,ST,SC,GPU) registers MR-T images for
%   different orientations
%   IMT is the the MR-T mask at all cardiac phases (orientations in the 3rd
%   dimension)
%   IMC is the MR-C mask at ED
%   ST is the MR-T resolution
%   SC is the MR-C resolution
%   GPU is a flag that determines whether to use gpu (1) or cpu (0) 
%   computation
%   It returns:
%   T, the estimated translations (in pixels)
%   TRED, the rounded estimated translations (in pixels)
%

tstart=tic;
NT=size(ImT);
NC=size(ImC);
S=NT(3);%Number of orientations
[xC(:,:,1),xC(:,:,2)]=ndgrid(1:sC(1):NC*sC(1),1:sC(2):NC*sC(2));
[xT(:,:,1),xT(:,:,2)]=ndgrid(1:sT(1):NT*sT(1),1:sT(2):NT*sT(2));
xC=single(xC);
xT=single(xT);
T=single(zeros(2,S));%Just translations are estimated... this should be improved
IT=single(zeros([NC(1:2) S]));%Transformed images for each orientation
NitMax=100;
H=single(zeros(NitMax,S));%Energy for each iteration
dx=single(zeros([NC(1:2) 2 S]));%Gradient of the image for each location
dH=single(zeros(2,S));%Gradient of the energy for each parameter
iter=1;
w=1000/(NC(1)*NC(2));
stopIt=0;

if gpu
    xC=gpuArray(xC);T=gpuArray(T);ImT=gpuArray(ImT);xT=gpuArray(xT);IT=gpuArray(IT);dx=gpuArray(dx);ImC=gpuArray(ImC);
end

while ~stopIt
   %Coordinate update
   xTact=bsxfun(@plus,xC,permute(T,[3 4 1 2]));%Transforms for each orientation
   %Transformed images interpolation
   for s=1:S
       IT(:,:,s)=interpn(xT(:,:,1),xT(:,:,2),ImT(:,:,s),xTact(:,:,1,s),xTact(:,:,2,s),'linear');
   end     
   %Metric calculation
   IT(isnan(IT))=0;
   V=bsxfun(@minus,IT,ImC).^2;%Metric for each point
   if ~gpu
       H(iter,:)=permute(sum(sum(V,1),2),[1 3 2]);   
   else
       H(iter,:)=gather(permute(sum(sum(V,1),2),[1 3 2]));
   end
   %Gradient calculation
   for s=1:S
       [dx(:,:,2,s),dx(:,:,1,s)]=gradient(IT(:,:,s));
   end
   dV=permute(2*bsxfun(@minus,IT,ImC),[1 2 4 3]);
   dV=bsxfun(@times,dV,dx);%Gradient of the metric for each point
   if ~gpu
       dH=permute(sum(sum(dV,1),2),[3 4 1 2]); 
   else
       dH=gather(permute(sum(sum(dV,1),2),[3 4 1 2]));
   end
   %Transform update
   dH=-w*dH;
   if ~gpu
       T=T+dH;
   else
       T=T+gpuArray(dH);
   end
   %Output conditions
   if iter~=1
      %[iter max(sqrt(sum(D.*D,1))) min(H(iter,:)-H(iter-1,:)) -min(H(1,:))/100]
      if (max(sqrt(sum(dH.^2,1)))<0.05 && min(H(iter,:)-H(iter-1,:))>-min(H(1,:))/200) || iter>=NitMax
        stopIt=1;
      end
   end
   iter=iter+1;    
end
if gpu
    T=gather(T);
end
Tred=bsxfun(@times,sT',round(bsxfun(@times,T,(sT').^(-1))));

telapsed=toc(tstart);
fprintf('Time registering MR-T: %f\n',telapsed)
