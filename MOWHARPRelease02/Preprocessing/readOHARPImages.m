function [ImC,ImT,eT,kT,orT]=readOHARPImages(folder,TagIm)

%READOHARPIMAGES   Reads OHARP images both for MR-C and MR-T 
%   [IMC,IMT,ET,KT,ORT]=READOHARPIMAGES(FOLDER,TAGIM)
%   reads images both for MR-C at ED and for MR-T at all cardiac phases 
%   (for the different orientations acquired)
%   FOLDER are a pair of paths to the MR-C and the MR-T data
%   TAGIM is a vector that contains the indexes of useful orientations
%   It returns:
%   IMC, the MR-C image at ED
%   IMT, the MR-T image at all cardiac phases (orientations in the 3rd
%   dimension, cardiac phases in the 4th dimension)
%   ET, DICOM info about MR-T
%   KT, spacing of the tags
%   ORT, orientation of the tags (in radians)
%

tstart=tic;
%Read MR-C image at ED
ImC=single(dicomread(sprintf('%s/im1',folder{1})));
S=length(TagIm);
eT=cell(1,S);
kT=single(zeros(1,S));
orT=single(zeros(1,S));
for s=1:S
    %Read MR-T images    
    ImT(:,:,s,:)=single(dicomread(sprintf('%s/im%d',folder{2},TagIm(s))));
    eT{s}=dicominfo(sprintf('%s/im%d',folder{2},TagIm(s)));    
    kT(s)=eT{s}.TagSpacingFirstDimension;
    orT(s)=eT{s}.TagAngleFirstAxis;
end
orT=pi*orT/180;
telapsed=toc(tstart);
fprintf('Time reading images: %f\n',telapsed)
