function  ErrARR  = PeakEvol( MTsBis, lambdas )

FontSize=12;
markers = {'o','s','d','^','>','v','<','p','h','x'};
lineas = {'-.','-',':','--'}; 

kk=[14:2:28 110 1102];
for k=1:8
    load(['MOWHARPStrainDataC0' num2str(kk(k))],'ErrM');
    for c=1:17
    noPoints=sum(MTsBis{c}(:));
            for t=1:4%Different methods                    
                    ErrRR{c}(:,:,t)=ErrM{1}{1}{c}{t}-ErrM{1}{1}{c}{end};
            end        
            ErrARR(c,:,k)=permute(sqrt(sum(sum(ErrRR{c}.^2,1),2)/noPoints),[5 3 4 1 2]);%Error
            BiasARR(c,:,k)=permute(sum(sum(ErrRR{c},1),2)/noPoints,[5 3 4 1 2]);%Error
            StdARR(c,:,k)=sqrt(ErrARR(c,:,k).^2-BiasARR(c,:,k).^2);
          
    end

end

A=squeeze(ErrARR(:,4,1:8));
B=squeeze(BiasARR(:,4,1:8));
C=squeeze(StdARR(:,4,1:8));
%A(:,9:10)=squeeze(ErrARR(:,2,9:10));

for i=1:17, leyenda{i}=['\gamma=' num2str(lambdas(i))]; end
leyenda{18}='First peak';
leyenda{19}='Second peak';

ColorSet=get(groot,'defaultAxesColorOrder');

pp=[1     6   8    11    14    16];
%pp=[1     4     7    10    13    16];
for i=1:6
    if i==1, resta=8;
    elseif i==2, resta=5;
    elseif i==3, resta=[4 3.5 4 4 4.33 4 4 4];
    elseif i==4, resta=[3 3 3 3 3.33 3.66 3 3];
    elseif i==5, resta=2;
    elseif i==6, resta=1; end
    plot(kk(1:8)/100-0.14,log(A(pp(i),1:8))-resta,'Color',ColorSet(i,:),'LineWidth',2,'Marker',markers{1} )
    hold on
end
legend(leyenda{1:3:17},'Location','Best')
grid on
xlabel('Distance between peaks (1/mm)','FontSize',12)
ylabel('log(MSE(ERR))','FontSize',12)
set(gca,'FontSize',25)
legend('boxoff')

figure(2)
for i=1:6
    plot(kk(1:8)/100-0.14,B(pp(i),1:8),'Color',ColorSet(i,:),'LineWidth',2,'Marker',markers{1} )
    hold on
end
legend(leyenda{1:3:17},'Location','SouthWest')
grid on
xlabel('Distance between peaks (1/mm)','FontSize',12)
ylabel('B(ERR)','FontSize',12)
set(gca,'FontSize',25)
legend('boxoff')

figure(3)
for i=1:6
    plot(kk(1:8)/100-0.14,log(C(pp(i),1:8)),'Color',ColorSet(i,:),'LineWidth',2,'Marker',markers{1} )
    hold on
end
legend(leyenda{1:3:17},'Location','SouthWest')
grid on
xlabel('Distance between peaks (1/mm)','FontSize',12)
ylabel('log(Var(ERR))','FontSize',12)
set(gca,'FontSize',25)
legend('boxoff')