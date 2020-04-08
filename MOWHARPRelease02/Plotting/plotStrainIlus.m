function plotStrainIlus( ErrM,S,C,FontSize,ErrLims,x,I,file,num,realData,NoOr )

for s=1:S+1   
    au=0;
    for c=2:2:C
        au=au+1;
        subtightplot(floor(C/2),S+1,(au-1)*(S+1)+s,[0 0])
        if s~=S+1        
            if realData
                imshow(ErrM{1}{s}{c}{num}(10:end-10,1+11+10:end-11-10),ErrLims)
            else
                %imshow(abs(ErrM{1}{1}{c}{s}-ErrM{1}{1}{c}{S}),ErrLims)
                imshow(ErrM{1}{1}{c}{s},ErrLims)
            end
            if c==2
                if realData
                    title(sprintf('I=%d',I(s)),'FontSize',FontSize)                    
                else
                    title(sprintf('%s',I{s}),'FontSize',FontSize)    
                end
            end
        %text(0.1*size(Err{s}{n},2),0.96*size(Err{s}{n},1),sprintf('H=%1.3f',entropy(Err{s}{n}(Err{s}{n}~=0))))
        %text(0.5*size(Err{s}{n},2),0.96*size(Err{s}{n},1),sprintf('eK=%1.3f',kurtosis(Err{s}{n}(Err{s}{n}~=0))-3))
        else
            imshow([],ErrLims)                    
            if realData
                text(0.65,0.5,strcat('\mu',sprintf('=%0.2f',x(c))),'FontSize',FontSize)
            else
                text(0.85,0.5,strcat('\mu',sprintf('=%0.2f',x(c))),'FontSize',FontSize)
            end
            colorbar('West')
            set(gca,'fontsize',FontSize,'Color',[1 1 1])
        end   
    end
end
%colormap('default')    
colormap('jet')   
%%export_fig(sprintf('/home/lcg13/Articulos/02EnProceso/MRM14/figs/Results/%s.png',file));
if realData
    set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);
    export_fig(sprintf('./Data/%s.png',file));
else
    %set(gcf,'Color',[1 1 1]);
    %export_fig(sprintf('/home/lcg13/Articulos/02EnProceso/MedIA15/UltimateFigs/Synth%s.pdf',file));
    tit=sprintf('%s-YeAM-I=%02d',file,NoOr);
    text(0.5,1,tit,'FontSize',FontSize+8)
    set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);
    export_fig(sprintf('./Synthetic/%s.png',tit));
end
