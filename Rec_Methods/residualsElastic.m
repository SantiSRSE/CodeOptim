function [r,xT,tR]=residualsElastic(rec)

%get residuals from the new transformation
%M conj(FDSTx -y)

[rec.T,Taux]=transformBSpline(rec.T0,rec.tS,rec.P); % obtain new mesh from control points

if any(rec.Field.lambda~=0);tR=regularizeBSpline(Taux,rec.tS,rec.Field);else, tR=[];end

[r,xT]=EncodingElastic(rec,rec.x,rec.T,1); %forward application
r=conj(r-rec.y);%Note we conjugate the residuals for convenience

r=bsxfun(@times,r,rec.A); %this mask could be only applied to y
