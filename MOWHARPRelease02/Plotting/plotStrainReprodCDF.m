function plotStrainReprodCDF(NormFrob,MTsBis,FontSize,expNo)

FontSize=12;
S=size(NormFrob,5);
T=size(NormFrob,3);
markers = {'o','s','d','^','>','v','<','p','h','x'};
ColorSet=get(groot,'defaultAxesColorOrder');
ColorSet(end+1:S,:)=ColorSet(1:S-size(ColorSet,1),3:-1:1);
lineas = {'-.','-'};

cont=1;
for s=1:S-1
    for t=1:T
        Datos=NormFrob(:,:,t,9,s);%Case 0.35 0.5 es 8
        Datos=abs(Datos(MTsBis));                
        XIA=linspace(0,0.5,1000);
        
        FA=ksdensity(Datos,XIA,'support','positive','function','cdf');
%         if t==1
%             au{cont}=plot(XIA,FA,':','Color',ColorSet(s,:),'LineWidth',2,'Marker',markers{s});         
%         elseif t==2
%             au{cont}=plot(XIA,FA,'-.','Color',ColorSet(s,:),'LineWidth',2,'Marker',markers{s});
%         else
%             au{cont}=plot(XIA,FA,'-','Color',ColorSet(s,:),'LineWidth',2,'Marker',markers{s});
%         end
        if t==1
            au{cont}=plot(XIA,FA,'-','Color',[0 191/255 1],'LineWidth',2,'Marker',markers{s});
        elseif t==2
            au{cont}=plot(XIA,FA,'-','Color',[1 69/255 0],'LineWidth',2,'Marker',markers{s});
        elseif t==3
            %au{cont}=plot(XIA,FA,'-','Color','b','LineWidth',2,'Marker',markers{s});
            au{cont}=plot(XIA,FA,'-','Color','r','LineWidth',2,'Marker',markers{s});
        elseif t==4
            %au{cont}=plot(XIA,FA,'-','Color','r','LineWidth',2,'Marker',markers{s});
            au{cont}=plot(XIA,FA,'-','Color','b','LineWidth',2,'Marker',markers{s});
        end
        cont=cont+1;
        hold on         
    end
end
numMark=20;
for n=1:length(au)
    nummarkers(au{n},numMark)
end
% if expNo==1
%     AX=legend('HARP I=2','WHARP I=2','AWHARP I=2','HARP I=3','WHARP I=3','AWHARP I=3',...
%         'HARP I=6','WHARP I=6','AWHARP I=6','HARP I=9','WHARP I=9','AWHARP I=9',...
%         'HARP I=18','WHARP I=18','AWHARP I=18','Location','SouthEast');
% else
%     AX=legend('HARP rep. mea.','WHARP rep. mea.','AWHARP rep. mea.',...
%         'HARP ext. ori.','WHARP ext. ori.','AWHARP ext. ori.','Location','SouthEast');
% end
if expNo==1
%     AX=legend('HARP I=2','WHARP I=2','MP-HARP I=4', 'MP-AWHARP I=4','HARP I=4','WHARP I=4','MP-HARP I=8','MP-AWHARP I=8',...
%         'HARP I=8','WHARP I=8','MP-HARP I=16','MP-AWHARP I=16','HARP I=9','WHARP I=9','MP-HARP I=18','MP-AWHARP I=18',...
%         'HARP I=18','WHARP I=18','MP-HARP I=36','MP-AWHARP I=36','Location','SouthEast');
      AX=legend('HARP I=2','WHARP I=2','MP-AWHARP I=4', 'MP-HARP I=4','HARP I=4','WHARP I=4',...
          'MP-AWHARP I=8','MP-HARP I=8', 'HARP I=8','WHARP I=8','MP-AWHARP I=16','MP-HARP I=16',...
          'HARP I=9','WHARP I=9','MP-AWHARP I=18','MP-HARP I=18','Location','SouthEast');
else
    AX=legend('HARP rep. mea.','WHARP rep. mea.','MP-HARP rep. mea.','MP-AWHARP rep. mea.',...
        'HARP ext. ori.','WHARP ext. ori.','MP-HARP ext. ori.','MP-AWHARP ext. ori.','Location','SouthEast');
end 
axis([0 0.5 0 1])
grid on
xlabel('FND','FontSize',FontSize)
ylabel('CDF','FontSize',FontSize)
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',FontSize)
legend('boxoff')
set(gca,'fontsize',FontSize,'Color',[1 1 1])
%set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);      
set(gcf,'Color',[1 1 1]);
%export_fig('/home/lcg13/Articulos/02EnProceso/MRM14/figs/Results/cdfFND.png');
export_fig(sprintf('./Data/cdfFND%d.png',expNo));
