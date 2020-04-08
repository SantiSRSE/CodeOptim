%Author: Michael Jacob Mathew; An Implementation of Nalder Mead Simplex
%Algorithm for n independent variables
%Reference: Lagarias, Jeffrey C., et al. "Convergence properties of the Nelder--Mead simplex method in low dimensions." 
%SIAM Journal on Optimization 9.1 (1998): 112-147.

%The required function should entered in the function "f" below 
%value of n should be modified; n= no of independent variables +1
%%
function [Minima]=NelMead()
clc;

StdVal=10; %any value for convergence

n=5; %value of N+1
P=1; %reflection coefficient
Chi=2; %expansion coefficient
Gamma=0.5; %contraction coefficient
Sig=0.5; %shrink coefficient

tic;
i=1;
minmas=0;
%coefAux1 =[1 zeros(1,11) 14.3];
%coefAux2 =[zeros(1,3) 1 zeros(1,8) 14.3];
%coefAux3 =[zeros(1,6) 1 zeros(1,5) 14.3];
%coefAux4 =[zeros(1,9) 1 zeros(1,2) 14.3];
%coefAux1 =[1 0 0 0 0 0];
%coefAux2 =[0 0 0 1 0 0];
%coefAux3 =[0 1 0 0 0 0];
%coefAux4 =[0 0 0 0 1 0];
BW=0.35;
NW=32;
Gdef=16
sT=[1 1];
[ImT,kT,orT,Tred,raC,MTs,limsOu,ErrMGT,EccMGT,ErcMGT,lambdas]=generatePhantomHARP2(1, [0 0 inf inf]);
ImTGT1=generatePhantomHARP2(1, [0 0 inf inf]);
ImTGT2=generatePhantomHARP2(1, [-pi/4 pi/4 kT/2 kT/2]);
%ImTGT3=generatePhantomHARP2(1, coefAux3);
%ImTGT4=generatePhantomHARP2(1, coefAux4);
fLGT=zeros(72,72,4,1,2); dfL=zeros(72,72,4,1,2);
size(localPhase(ImTGT1{Gdef},'WFT',BW,sT,kT,orT-pi/4,limsOu{Gdef},Tred,NW,0))
[fLGT(:,:,1,:,:), lims]=localPhase(ImTGT1{Gdef},'WFT',BW,sT,kT,orT-pi/4,limsOu{Gdef},Tred,NW,0);
fLGT(:,:,2,:,:)=localPhase(ImTGT1{Gdef},'WFT',BW,sT,kT,orT+pi/4,limsOu{Gdef},Tred,NW,0);
fLGT(:,:,3,:,:)=localPhase(ImTGT2{Gdef},'WFT',BW,sT,kT,orT-pi/4,limsOu{Gdef},Tred,NW,0);
fLGT(:,:,4,:,:)=localPhase(ImTGT2{Gdef},'WFT',BW,sT,kT,orT+pi/4,limsOu{Gdef},Tred,NW,0);
dfL(:,:,1,:,:)=gradientWrapping(fLGT(:,:,1,:,:),sT); dfL(:,:,1,:,:)=MRTTransform(dfL(:,:,1,:,:),sT,Tred);
dfL(:,:,2,:,:)=gradientWrapping(fLGT(:,:,2,:,:),sT); dfL(:,:,2,:,:)=MRTTransform(dfL(:,:,2,:,:),sT,Tred);
dfL(:,:,3,:,:)=gradientWrapping(fLGT(:,:,3,:,:),sT); dfL(:,:,3,:,:)=MRTTransform(dfL(:,:,3,:,:),sT,Tred);
dfL(:,:,4,:,:)=gradientWrapping(fLGT(:,:,4,:,:),sT); dfL(:,:,4,:,:)=MRTTransform(dfL(:,:,4,:,:),sT,Tred);
kv(1,:)=2*pi*[-sin(orT-pi/4) cos(orT-pi/4)]/kT; kv(2,:)=2*pi*[-sin(orT+pi/4) cos(orT+pi/4)]/kT;
kv(3,:)=2*pi*[-sin(orT-pi/4) cos(orT-pi/4)]/kT; kv(4,:)=2*pi*[-sin(orT+pi/4) cos(orT+pi/4)]/kT;
fL2=IRWLSMaterial(kv,dfL(:,:,:,1,:),MTs{Gdef},lims,0.01);
    
SortVertices = CreateInitialSimplex(n,fL2);

while(StdVal >= 10^(0))
   
    SortVertices = BubbleSort(SortVertices,n);
    SortVertices.coord
    SortVertices.value
    Centroid = CalculateCentroid(SortVertices,n,fL2);
    StdVal = CalculateStd(SortVertices,n)
    minmas(i)=SortVertices(1).value; i=i+1;
    %Simplex_2D(SortVertices); drawnow; 
    
    Reflect.coord = (1+P).*Centroid.coord - P.*SortVertices(n).coord; %Reflect
    Reflect.value = f(Reflect,fL2);
     
    if(SortVertices(1).value <= Reflect.value && Reflect.value < SortVertices(n-1).value)
        SortVertices(n)=Reflect;
        continue; %if the above criterion is sattisfied, then terminate the iteration
    end
    
    if(Reflect.value < SortVertices(1).value) %Expand
        
        Expand.coord = (1-Chi).*Centroid.coord+Chi.*Reflect.coord;
        Expand.value = f(Expand,fL2);
        
        if(Expand.value < Reflect.value)
            
            SortVertices(n) = Expand;
            continue;
        else
            SortVertices(n) = Reflect;
            continue;
        end 
    end
    
    if(SortVertices(n-1).value <= Reflect.value) %Contract
                  
            ContractOut.coord = (1-Gamma).*Centroid.coord + Gamma.*Reflect.coord; %Contract Outside
            ContractOut.value = f(ContractOut,fL2);
            
            ContractIn.coord = (1-Gamma).*Centroid.coord + Gamma.*SortVertices(n).coord;  %Contract Inside
            ContractIn.value= f(ContractIn,fL2);
            
            if(SortVertices(n-1).value <= Reflect.value && Reflect.value < SortVertices(n).value && ContractOut.value <= Reflect.value)
                SortVertices(n) = ContractOut;
                continue;
            elseif(SortVertices(n).value <= Reflect.value && ContractIn.value < SortVertices(n).value) %Contract Inside
                SortVertices(n) = ContractIn;
                continue;
            else
                for i=2:n
                    SortVertices(i).coord = (1-Sig).*SortVertices(1).coord + Sig.*SortVertices(i).coord;
                    SortVertices(i).value = f(SortVertices(i),fL2);   
                end
            end
    end    
end
toc
Minima=SortVertices(1);
save MinimaCFinal16 Minima
ImTF=generatePhantomHARP2(1, Minima.coord');
imshow(log(abs(fftshift(fft2(ImTF{Gdef})))),[])
% plot(minmas);
% a=text(50,250, strcat('No: of iterations =  ', num2str(i)));
% set(a, 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 18, 'Color', 'red');
% xlabel('iterations');
% ylabel('error');
% title('Error Vs Iterations');
end

%%
function [value]=f(V,fL2) %Write your function in matrix form
%in case of any of the three trial functions un comment two lines  

    coef   =V.coord;
    Gdef=16;
        
    [ImT,kT,orT,Tred,raC,MTs,limsOu,ErrMGT,EccMGT,ErcMGT,lambdas]=generatePhantomHARP2(1, coef');
        
    BW=0.35;
    NW=32;
    sT=[1 1];
    
    fL1=zeros(72,72,4,1,2); dfL=zeros(72,72,4,1,2);
    fL1(:,:,1,:,:)=localPhase(ImT{Gdef},'WFT',BW,sT,kT,orT-pi/4,limsOu{Gdef},Tred,NW,0);
    fL1(:,:,2,:,:)=localPhase(ImT{Gdef},'WFT',BW,sT,kT,orT+pi/4,limsOu{Gdef},Tred,NW,0);
    fL1(:,:,3,:,:)=localPhase(ImT{Gdef},'WFT',BW,sT,coef(3),orT+coef(1),limsOu{Gdef},Tred,NW,0);
    [fL1(:,:,4,:,:), lims]=localPhase(ImT{Gdef},'WFT',BW,sT,coef(4),orT+coef(2),limsOu{Gdef},Tred,NW,0);
    dfL(:,:,1,:,:)=gradientWrapping(fL1(:,:,1,:,:),sT); dfL(:,:,1,:,:)=MRTTransform(dfL(:,:,1,:,:),sT,Tred);
    dfL(:,:,2,:,:)=gradientWrapping(fL1(:,:,2,:,:),sT); dfL(:,:,2,:,:)=MRTTransform(dfL(:,:,2,:,:),sT,Tred);
    dfL(:,:,3,:,:)=gradientWrapping(fL1(:,:,3,:,:),sT); dfL(:,:,3,:,:)=MRTTransform(dfL(:,:,3,:,:),sT,Tred);
    dfL(:,:,4,:,:)=gradientWrapping(fL1(:,:,4,:,:),sT); dfL(:,:,4,:,:)=MRTTransform(dfL(:,:,4,:,:),sT,Tred);
    kv(1,:)=2*pi*[-sin(orT) cos(orT)]/kT; kv(2,:)=2*pi*[-sin(orT+pi/2) cos(orT+pi/2)]/kT;
    kv(3,:)=2*pi*[-sin(orT+coef(1)) cos(orT+coef(1))]/coef(3); 
    kv(4,:)=2*pi*[-sin(orT+coef(2)) cos(orT+coef(2))]/coef(4);
    fL=IRWLSMaterial(kv,dfL(:,:,:,1,:),MTs{Gdef},lims,0.01);
    
    mediana=sqrt(sum(sum((fL-fL2).^2,4),5)); 
    median(mediana(MTs{Gdef}(lims(1,1):lims(2,1),lims(1,2):lims(2:2))>0.5))
    sum((fL(:)-fL2(:)).^2)
    
    value=sum((fL(:)-fL2(:)).^2);
    
    %x=V.coord(1); y=V.coord(2);  %Rosenbrock's Function
    %value=(1-x)^2+100*(y-x^2)^2;
    
%    x1=V.coord(1);x2=V.coord(2);x3=V.coord(3);x4=V.coord(4); %Powell's Quadratic Function
%    value=(x1+10*x2)^2+5*(x3-x4)^2+(x2-2*x3)^4+10*(x1-x4)^4;
      
%    x1=V.coord(1); x2=V.coord(2); x3=V.coord(3); %Fletcher and Powell's Helical Valley Function
%     theta=atan2(x1,x2);
%     value=100*(x3-10*theta)^2+(sqrt(x1^2+x2^2)-1)^2+x3^2;
    
%   value=([1,-1;0,1]*V.coord-[4;1])'*V.coord; %f(x,y)=x^2-4*x+y^2-y-x*y;
%   value=([1,0,-1;0,1,0;1,0,1]*V.coord-[3;10;-7])'*V.coord; %f(x,y,z)=x^2+y^2+z^2-3*x-10*y+7
    
end
%%
function [Vertices]=CreateInitialSimplex(n,fL2)

    %ExpectMin=rand(n-1,1).*2-1;
    ExpectMin=[0 pi/2 7.15 7.15]';
    Vertices(1).coord=ExpectMin; % expected minima
    Vertices(1).value=f(Vertices(1),fL2);
    ExpectMin=[0 pi/2 7.15/2 7.15/2]';
    Vertices(2).coord=ExpectMin; % expected minima
    Vertices(2).value=f(Vertices(2),fL2);
    ExpectMin=[-pi/4 pi/4 7.15 7.15]';
    Vertices(3).coord=ExpectMin; % expected minima
    Vertices(3).value=f(Vertices(3),fL2);
    ExpectMin=[-pi/4 pi/4 7.15/2 7.15/2]';
    Vertices(4).coord=ExpectMin; % expected minima
    Vertices(4).value=f(Vertices(4),fL2);
    ExpectMin=[-pi/4 pi/4 7.15/1.5 1.15/1.5]';
    Vertices(5).coord=ExpectMin; % expected minima
    Vertices(5).value=f(Vertices(5),fL2);
    %for i=2:n
    %    Vertices(i).coord=ExpectMin+(rand(n-1,1)-0.5).*[pi/2; pi/2; 2; 2 ]*2; %100 is the scale factor
    %    Vertices(i).value=f(Vertices(i),fL2);
    %end
    
end
%%
function [SortVertices]=BubbleSort(Vertices,n)
    
    while(n~=0)
        for i=1:n-1
            if(Vertices(i).value<=Vertices(i+1).value)
                continue;
            else
                temp=Vertices(i);
                Vertices(i)=Vertices(i+1);
                Vertices(i+1)=temp;
            end
        end
        n=n-1;          
    end
    SortVertices=Vertices;
    
end
%%
function [Centroid]=CalculateCentroid(Vertices,n,fL2)
    
    Sum=zeros((n-1),1);
    for i=1:n-1
        Sum=Sum+Vertices(i).coord;
    end
    Centroid.coord=Sum./(n-1);
    Centroid.value=f(Centroid,fL2);
    
end
%%
function[StdVal]=CalculateStd(Vertices,n) % this is the tolerance value, the standard deviation of the converging values
    
    for i=1:n
        ValueArray(i)=Vertices(i).value;
    end
    StdVal=std(ValueArray,1);
    
end
%%