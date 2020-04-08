function dfL=gradientWrapping(fL,sT)

%GRADIENTWRAPPING   Computes a wrapped gradient 
%   DFL=gradientWrapping(FL,ST) computes the wrapped gradient of the local
%   phase
%   FL is the estimated local phase
%   It returns:
%   DFL, the wrapped gradient of the local phase
%

tstart=tic;
M=size(fL);
dfL=single(zeros(M));
dfW=cell(1,2);dfWAux=cell(1,2);
for n=1:2
    for m=1:2
        %dfW{m}=diff(fL(:,:,:,:,m),1,3-n)/sT(n);
        dfW{m}=diff(fL(:,:,:,:,m),1,n)/sT(n);
        dfWAux{m}=abs(dfW{m});        
    end    
    dfWBis=dfW{1};            
    dfWBis(dfWAux{1}>dfWAux{2})=dfW{2}(dfWAux{1}>dfWAux{2});     
    dfWBisA=single(zeros(M(1:4)));
    dfWBisB=single(zeros(M(1:4)));
    if n==1
        dfWBisA(1:M(1)-1,:,:,:)=dfWBis;
        dfWBisB(2:M(1),:,:,:)=dfWBis;
    else
        dfWBisA(:,1:M(2)-1,:,:)=dfWBis;
        dfWBisB(:,2:M(2),:,:)=dfWBis;
    end
    dfL(:,:,:,:,n)=(dfWBisA+dfWBisB)/2;    
end

telapsed=toc(tstart);
fprintf('Gradient wrapping time: %f\n',telapsed)
