function plotStrainError(ErrE,lambdas,tit,useLog)

C=size(ErrE,1);%Deformations
S=size(ErrE,2);%Stripes
T=size(ErrE,3);%HARP/WHARP
I=[2 3 6 9 18];

ColorSet=get(groot,'defaultAxesColorOrder');
ColorSet(end+1:S,:)=ColorSet(1:S-size(ColorSet,1),3:-1:1);

FontSize=12;
markers = {'o','s','d','^','>','v','<','p','h','x'};
lineas = {'-.','-',':','--'}; 
if ~useLog
    sinic=2;
else
    sinic=1;
end
for s=sinic:S
    for t=1:T
        erro=ErrE(:,s,t);
        if t==1
            if useLog %[0 191/255 1]
                plot(lambdas,log(erro),lineas{s},'Color','b','LineWidth',2,'Marker',markers{s})
            else
                plot(lambdas,erro,'Color','b','LineWidth',2,'Marker',markers{s})
            end                
        elseif t==2
            if useLog
                plot(lambdas,log(erro),lineas{s},'Color',[1 69/255 0],'LineWidth',2,'Marker',markers{s})
            else
                plot(lambdas,erro,lineas{s},'Color',[1 69/255 0],'LineWidth',2,'Marker',markers{s})
            end
        elseif t==3
            if useLog
                %plot(lambdas,log(erro),lineas{s},'Color','b','LineWidth',2,'Marker',markers{s})
            else
                %plot(lambdas,erro,lineas{s},'Color','b','LineWidth',2,'Marker',markers{s})
            end
        elseif t==4
            if useLog
                plot(lambdas,log(erro),lineas{s},'Color','r','LineWidth',2,'Marker',markers{s})
            else
                plot(lambdas,erro,lineas{s},'Color','r','LineWidth',2,'Marker',markers{s})
            end
        elseif t==5
            if useLog
                plot(lambdas,log(erro),lineas{s},'Color','k','LineWidth',2,'Marker',markers{s})
            else
                plot(lambdas,erro,lineas{s},'Color','k','LineWidth',2,'Marker',markers{s})
            end
        end
        hold on            
    end
end
if ~useLog
    %AX=legend('HARP','WHARP','MP-HARP','MP-AWHARP','MOHARP','MOWHARP','MOP-HARP','MOP-AWHARP','Location','NorthWest');    
else   
    %AX=legend('HARP','WHARP','MP-HARP','MP-AWHARP','MOHARP','MOWHARP','MOP-HARP','MOP-AWHARP','Location','NorthWest');
end
AX=legend('HARP','WHARP','MP-WHARP','MOHARP','MOWHARP','MOP-WHARP','Location','NorthWest');
%AX=legend('MOHARP','MOWHARP','MOP-WHARP','Location','SouthWest');
grid on
xlabel('\gamma','FontSize',FontSize)
if ~useLog
    ylabel(sprintf('B(%s)',tit),'FontSize',FontSize);
else
    ylabel(sprintf('log(MSE(%s))',tit),'FontSize',FontSize)    
end
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',FontSize)
legend('boxoff')
set(gca,'fontsize',FontSize,'Color',[1 1 1])
%set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);      
set(gcf,'Color',[1 1 1]);      
%export_fig('/home/lcg13/Articulos/02EnProceso/MRM14/figs/Results/MedianFNDBW.png');
%export_fig(sprintf('/home/lcg13/Articulos/02EnProceso/MedIA15/UltimateFigs/MedianFNDBW%d.pdf',expNo));

