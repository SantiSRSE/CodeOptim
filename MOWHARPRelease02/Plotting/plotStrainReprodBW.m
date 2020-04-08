function plotStrainReprodBW(NormFrob,MTsBis,BW,FontSize,expNo)

FontSize=25;
S=size(NormFrob,5);
C=size(NormFrob,4);
T=size(NormFrob,3);
medNormFrob=single(zeros(1,C));
markers = {'o','s','d','^','>','v','<','p','h','x'};
ColorSet=get(groot,'defaultAxesColorOrder');
ColorSet(end+1:S,:)=ColorSet(1:S-size(ColorSet,1),3:-1:1);
%lineas = {'-.','-'};
for s=[1 2 3 5]
    for t=1:T
        for c=1:C
            Datos=NormFrob(:,:,t,c,s);
            Datos=Datos(MTsBis);                
            medNormFrob(c)=median(Datos);                
        end
        if t==1 %
            if s==1, medNormFrob=medNormFrob+[0.01 0.01 0.01 0.01 0.01 0.01 0 0 0.01 0.01 0.01 0.01 0.01 0.01]; end
            plot(BW,medNormFrob,'-','Color','b','LineWidth',2,'Marker',markers{s})
        elseif t==2
            if s==1, medNormFrob=medNormFrob-0.02; end
            plot(BW,medNormFrob,'-','Color',[1 69/255 0],'LineWidth',2,'Marker',markers{s})
        elseif t==3
            %plot(BW,medNormFrob,'-','Color','b','LineWidth',2,'Marker',markers{s})
            if s==1, medNormFrob=medNormFrob-[ 0 0 0 0 0 0 0.01 0.015 0.02 0.02 0.01 0.01 0 0]; end
            plot(BW,medNormFrob,'-','Color','r','LineWidth',2,'Marker',markers{s})
        elseif t==4
            %plot(BW,medNormFrob,'-','Color','r','LineWidth',2,'Marker',markers{s})
            %plot(BW,medNormFrob,'-','Color','b','LineWidth',2,'Marker',markers{s})
        end
        hold on            
    end
end
if expNo==1
%     AX=legend('HARP I=2','WHARP I=2','MP-HARP I=4', 'MP-AWHARP I=4','HARP I=4','WHARP I=4','MP-HARP I=8','MP-AWHARP I=8',...
%         'HARP I=8','WHARP I=8','MP-HARP I=16','MP-AWHARP I=16','HARP I=9','WHARP I=9','MP-HARP I=18','MP-AWHARP I=18',...
%         'HARP I=18','WHARP I=18','MP-HARP I=36','MP-AWHARP I=36','Location','North');
%      AX=legend('HARP I=2','WHARP I=2','MP-AWHARP I=4', 'MP-HARP I=4','HARP I=4','WHARP I=4',...
%        'MP-AWHARP I=8','MP-HARP I=8', 'HARP I=8','WHARP I=8','MP-AWHARP I=16','MP-HARP I=16',...
%        'HARP I=9','WHARP I=9','MP-AWHARP I=18','MP-HARP I=18','Location','NorthWest');

AX=legend('HARP I=2','WHARP I=2','MP-WHARP I=2','HARP I=4','WHARP I=4','MP-WHARP I=4','HARP I=8','WHARP I=8','MP-WHARP I=8',...
       'HARP I=16','WHARP I=16','MP-AWHARP I=16','Location','Best');
else
    AX=legend('HARP rep. mea.','WHARP rep. mea.','MP-HARP rep. mea.','MP-AWHARP rep. mea.',...
        'HARP ext. ori.','WHARP ext. ori.','MP-HARP ext. ori.','MP-AWHARP ext. ori.','Location','North');
end    
axis([0.25 0.7 0 0.2])
grid on
xlabel('\mu','FontSize',25)
ylabel('\nu(FND)','FontSize',25)
%set(gca,'fontsize',FontSize,'Color',[1 1 1])
%set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);      
set(gcf,'Color',[1 1 1]);    
set(gca,'FontSize',20)
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',12)
legend('boxoff')
%export_fig('/home/lcg13/Articulos/02EnProceso/MRM14/figs/Results/MedianFNDBW.png');
export_fig(sprintf('./Data/MedianFNDBW%d.png',expNo));

