function [ c ] = ConvField( T1,T2,M )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

edT=sqrt(sum(abs(T1-T2).^2,4)); %rmse


c=gather(multDimMax(bsxfun(@times,edT,M),1:3));
% NN=size(edT);NN(end+1:5)=1;
% if strcmp(rec.typEnc,'sequential') 
%     edT=reshape(edT,[NN(1) NN(2)*NN(5)]);
% else
%     Ns=sqrt(NN(5));
%     edT=reshape(edT,[NN(1:2) Ns Ns]);
%     edT=permute(edT,[1 3 2 4]);
%     edT=reshape(edT,Ns*NN(1:2));
% end
% rang=[min(edT(:)) max(edT(:))];
% c=10.^(gather(linspace(log10(rang(1)),log10(rang(2)),6)));c=round(c*10000)/10000;
% figure
% imshow(log10(edT),log10(rang))
% colormap(parula);  %Color palate named "bone"
% caxis(gather(log10(rang)));
% colorbar
% colorbar('FontSize',8,'YTick',log10(c),'YTickLabel',c);
% title('Euclidean error transformations (log10 scale)')
% set(gcf,'Position',get(0,'ScreenSize'))
% hold on
% v=[0.1 0.2 0.5 1 2 5 10];
% [C,h]=contour(edT,v,'Color','k');
% clabel(C,h,v,'Color','k','FontWeight','bold')

end

