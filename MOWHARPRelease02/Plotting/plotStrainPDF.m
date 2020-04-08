function plotStrainPDF( E,MTsBis,S,T,FontSize,range,file,label,expNo,loca)

FontSize=12;
cont=1;
markers = {'o','s','d','^','>','v','<','p','h','x'};
ColorSet=get(groot,'defaultAxesColorOrder');
ColorSet(end+1:S,:)=ColorSet(1:S-size(ColorSet,1),3:-1:1);
lineas = {'-.','-'};

for ti=1:2
    for s=1:S
        for t=1:T                            
            Datos=E{ti}{s}{6}{t};        %10
            Datos=Datos(MTsBis);                
            XIA=linspace(range(1),range(2),1000);           
            FA=ksdensity(Datos,XIA,'function','pdf');
%             if t==1
%                 au{cont}=plot(XIA,FA,':','Color',ColorSet(s,:),'LineWidth',2,'Marker',markers{s});
%             elseif t==2
%                 au{cont}=plot(XIA,FA,'-.','Color',ColorSet(s,:),'LineWidth',2,'Marker',markers{s});
%             else
%                 au{cont}=plot(XIA,FA,'-','Color',ColorSet(s,:),'LineWidth',2,'Marker',markers{s});
%             end
            if t==1
                au{cont}=plot(XIA,FA,lineas{2},'Color',[0 191/255 1],'LineWidth',2,'Marker',markers{s});
            elseif t==2
                au{cont}=plot(XIA,FA,lineas{2},'Color',[1 69/255 0],'LineWidth',2,'Marker',markers{s});
            elseif t==3
                au{cont}=plot(XIA,FA,lineas{2},'Color','b','LineWidth',2,'Marker',markers{s});
            else
                au{cont}=plot(XIA,FA,lineas{2},'Color','r','LineWidth',2,'Marker',markers{s});
            end
            if ti==1
                le(cont)=au{cont};
                cont=cont+1;
                hold on 
            else
                cont=cont+1;
            end
        end
    end      
end
numMark=10;
for n=1:length(au)
    nummarkers(au{n},numMark)
end
% if expNo==1
%     AX=legend('HARP I=2','WHARP I=2','AWHARP I=2','HARP I=3','WHARP I=3','AWHARP I=3', ...
%         'HARP I=6','WHARP I=6','AWHARP I=6','HARP I=9','WHARP I=9','AWHARP I=9',...
%         'HARP I=18','WHARP I=18','AWHARP I=18','Location',loca);
% else
%     AX=legend('HARP rep. mea.','WHARP rep. mea.','AWHARP rep. mea.',...
%         'HARP ext. ori.','WHARP ext. ori.','AWHARP ext. ori.','Location',loca);
% end
if expNo==1
    AX=legend('HARP I=2','WHARP I=2','MP-HARP I=6', 'MP-AWHARP I=6','HARP I=3','WHARP I=3','MP-HARP I=9','MP-AWHARP I=9',...
        'HARP I=6','WHARP I=6','MP-HARP I=18','MP-AWHARP I=18','HARP I=9','WHARP I=9','MP-HARP I=27','MP-AWHARP I=27',...
        'HARP I=18','WHARP I=18','MP-HARP I=54','MP-AWHARP I=54','Location','North');
else
    AX=legend('HARP rep. mea.','WHARP rep. mea.','MP-HARP rep. mea.','MP-AWHARP rep. mea.',...
        'HARP ext. ori.','WHARP ext. ori.','MP-HARP ext. ori.','MP-AWHARP ext. ori.','Location','North');
end
grid on
xlabel(sprintf('%s',label),'FontSize',FontSize,'Interpreter','latex')
ylabel('PDF','FontSize',FontSize)
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',FontSize)
legend('boxoff')
set(gca,'fontsize',FontSize,'Color',[1 1 1])
%set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1])
set(gcf,'Color',[1 1 1])
xlim(range);
%export_fig(sprintf('/home/lcg13/Articulos/02EnProceso/MRM14/figs/Results/%sPDF.png',file));
export_fig(sprintf('./Data/%sPDF%d.png',file,expNo));

