function [Err,Ecc Erc]=ejectionStrain(f,MTs,lims,raC,sp)

%EJECTIONSTRAIN   Computes the ejection Green strain tensor 
%   [ERR,ECC ERC]=ejectionStrain(F,LIMS,MTS,RAC,SP) computes the ejection
%   Green strain tensor form the deformation strain tensor
%   F is the estimated deformation gradient tensor
%   MTS is the MR-C mask in the MR-T space
%   LIMS is the limitant region for the myocardium
%   RAC, the radial direction coordinate vector
%   SP, flag to indicate whether the spatial or the material deformation
%   gradient tensor is introduced as input
%   It returns:
%   ERR, the RR component of the ejection Green strain tensor
%   ECC, the CC component of the ejection Green strain tensor
%   ERC, the RC component of the ejection Green strain tensor
%

f=shiftdim(f,3);
ciC(:,:,1)=-raC(:,:,2);ciC(:,:,2)=raC(:,:,1);
raC=shiftdim(raC,2);ciC=shiftdim(ciC,2);
M=size(f);
N=size(raC);
Err=single(zeros(N(2:3)));Ecc=single(zeros(N(2:3)));Erc=single(zeros(N(2:3)));
I=eye(2);
for m=1:2
    indi{m}=lims(1,m):lims(2,m);
end
for m=1:M(3)
    for n=1:M(4)
        if MTs(indi{1}(m),indi{2}(n))>0.5
            faux=f(:,:,m,n);
            if sp
                F=inv(faux);
            else
                F=faux;
            end
            E=0.5*(F'*F-I);
            raux=raC(:,indi{1}(m),indi{2}(n));
            caux=ciC(:,indi{1}(m),indi{2}(n));       
            Err(indi{1}(m),indi{2}(n))=raux'*E*raux;
            Ecc(indi{1}(m),indi{2}(n))=caux'*E*caux;
            Erc(indi{1}(m),indi{2}(n))=raux'*E*caux;            
        end
    end
end

Err=Err(lims(1,1):lims(2,1),lims(1,2):lims(2,2));
Ecc=Ecc(lims(1,1):lims(2,1),lims(1,2):lims(2,2));
Erc=Erc(lims(1,1):lims(2,1),lims(1,2):lims(2,2));
