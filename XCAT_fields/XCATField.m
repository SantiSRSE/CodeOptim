function [ TB ] = XCATField( Nshots, sm ,Dimen,gpu )

if Dimen==2
    C = LoadXCATField('../../data/KCL2D_vec_frame1_to_frame2.txt', 1,1,0);
    smV=[256 256 1 1 Nshots]; %tamano de campos de desplazamiento originales
elseif Dimen==3
    C = LoadXCATField('../../data/3D512/KCL3D256_vec_frame1_to_frame2.txt', [],1,0);
    smV=[256 256 256 1 Nshots]; %tamano de campos de desplazamiento originales
end
%Mult=(rand(1,Nshots)-0.5)*sm(1)/128; %multiplicador para campo de desplazamiento de cada shot aleatorio
Mult=linspace(-1,1,Nshots)*sm(1)/256;
TB=zeros([1 1 1 1 Nshots sm(1:3) 1 3],'single'); % campo de desplazamiento de cada shot
if gpu>0, TB=gpuArray(TB);end
for i=1:Nshots
    [~,TBplus] = ResampleField( C, smV,Mult(i),sm,Dimen,gpu); %leo el campo de desplazamiento de xcat
    TB(1,1,1,1,i,:,:,:,1,:)=single(dynInd(TBplus,2,5));  %asigno el campo de desplazamiento correspondiente a cada shot
    clear TBplus
end
clear C
TB=permute(TB,[6 7 8 10 5 1 2 3 4 9]);


