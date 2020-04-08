function M=buildShotMask(NY,Nsh,mode,gpu)
%Even and odd shots
%mode='sequential';
%mode='checker';
%mode='Randomchecker';

NY(end+1:4)=1;
M=single(zeros([NY Nsh]));
if any(NY(1:3)==1), Dimen=2; else, Dimen=3; end

if strcmp(mode,'sequential')
    for ns=1:Nsh
        M(ns:Nsh:end,:,1,1,ns)=1;
    end
elseif strcmp(mode,'checker')
    s=1;
    for m=1:sqrt(Nsh)
        BlockM=sqrt(Nsh);
        for n=1:sqrt(Nsh)
            BlockN=sqrt(Nsh);
            M(m:BlockM:end,n:BlockN:end,1,1,s)=1;
            s=s+1;
        end
    end
elseif strcmp(mode,'randomchecker')%THIS IS NOT READY PROBABLY, REQUIRES CHANGES AS THE DIMENSIONALITY IS NOT MATCHED
    NBlocksM=NY(1)/sqrt(Nsh);
    NBlocksN=NY(2)/sqrt(Nsh);
    for o=1:NBlocksM
        for p=1:NBlocksN
            indRand=randperm(Nsh);
            [IR,JR]=ind2sub([sqrt(Nsh) sqrt(Nsh)],indRand);
            s=1;
            for m=1:sqrt(Nsh)
                BlockM=sqrt(Nsh);
                for n=1:sqrt(Nsh)
                    BlockN=sqrt(Nsh);
                    M(IR(s)+(o-1)*BlockM,JR(s)+(p-1)*BlockN,1,1,s)=1;
                    s=s+1;
                end
            end
        end
    end 
else
    error('Mode %s not identified',mode);
end

% If dimension 2 do not exist, move to the third
if Dimen==2, if NY(2)==1, M=permute(M,[1 3 2 4 5]); end, end
% if 3D assume no subsampling in third dimension
if Dimen==3
    M=M(:,:,1,1,:); 
    %M=repmat(M(:,:,1,1,:),[1 1 NY(3) 1 1]); 
end

if gpu;M=gpuArray(M);end
