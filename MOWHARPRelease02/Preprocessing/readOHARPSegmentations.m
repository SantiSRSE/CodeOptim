function [MCd,MCs,MT]=readOHARPSegmentations(folder,tsis,TagIm)

%READOHARPSEGMENTATIONS   Reads OHARP segmentations both for MR-C and MR-T 
%   [MCD,MCS,MT]=READOHARPSEGM(FOLDER,TSIS,TAGIM)
%   reads segmentations both for MR-C at ED and ES and for MR-T at ED (for 
%   the different orientations acquired)
%   FOLDER are a pair of paths to the MR-C and the MR-T data
%   TSIS indexes the systolic phase
%   TAGIM is a vector that contains the indexes of useful orientations
%   It returns a set of masks for the segmentations:
%   MCD, the MR-C mask at ED
%   MCS, the MR-C mask at ES
%   MT, a 3D array of MR-T masks at ED for the different orientations (3rd
%   dimension)
%

tstart=tic;
%Read MR-C segmentations
MCd=single(dicomread(sprintf('%s/mioc1',folder{1}))~=0);%Diastole
MCs=single(dicomread(sprintf('%s/mioc%d',folder{1},tsis))~=0);%Systole

%Read MR-T segmentations
for s=1:length(TagIm)
    MT(:,:,s)=single(dicomread(sprintf('%s/mioc%d',folder{2},TagIm(s)))~=0);    
end
telapsed=toc(tstart);
fprintf('Time reading segmentations: %f\n',telapsed)
