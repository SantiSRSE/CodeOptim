function [r,MTs]=gridCalculation(MCs,MT,sC,sT)

%GRIDCALCULATION   Computes the orientation coordinates and the MR-C
%segmentation in the MR-T space
%   [R,MTS]=GRIDCALCULATION(MCS,MT,SC,ST) computes the orientation
%   coordinates and the MR-C segmentation in the MR-T space
%   MCS is the MR-C mask at ES
%   MT is the the MR-T mask at ED (orientations in the 3rd dimension)
%   ST is the MR-T resolution
%   SC is the MR-C resolution
%   It returns:
%   R, the radial direction coordinate vector
%   MTS, the MR-C mask in the MR-T space
%

NC=size(MCs);
NT=size(MT);

%Coordinate systems
[xC(:,:,1),xC(:,:,2)] =ndgrid(1:sC(1):NC*sC(1),1:sC(2):NC*sC(2));
[xT(:,:,1),xT(:,:,2)] =ndgrid(1:sT(1):NT*sT(1),1:sT(2):NT*sT(2));

%Transforms the cine segmentation to the tagging coordinate space
MTs=interpn(xC(:,:,1),xC(:,:,2),MCs,xT(:,:,1),xT(:,:,2),'linear');

sAux=regionprops(MTs>0.5, 'centroid');
s=sAux.Centroid(2:-1:1);

[rAux(:,:,1),rAux(:,:,2)]=ndgrid([1-s(1):NT(1)-s(1)],[1-s(2):NT(2)-s(2)]);
%Orientation vector
r=bsxfun(@times,rAux,sum(rAux.^2,3).^(-0.5));
