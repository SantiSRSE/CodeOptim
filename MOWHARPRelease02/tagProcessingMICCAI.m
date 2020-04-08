addpath(genpath('.'));

%Experiment 1: reproducibility of increased number of orientations
%Experiment 2: reproducibility of repeated measurements versus extended
%orientations
%Experiment 3: visual results
%Experiment 4: interferences
%Experiment 5: error with no interferences
Coef=2;

NORS=10; NORSof=16;
folder{1}='C:\Teleco\0DOCTORADO\MOWHARPRelease02\Data\EC-C';%Path to MR-C
folder{2}='C:\Teleco\0DOCTORADO\MOWHARPRelease02\Data\EC-T';%Path to MR-T

gpu=0;
fftw('planner','measure');

preprocessing=0;
processing=0;
postprocessing=0;
showplotE=1;

realData=1;

sC=[1.25 1.25];%MR-C resolution
sT=[4/3 4/3];%MR-T resolution

%tsisC=6;%MOWHARP
tsisC=12;%MOWHARP2
tsisT=7;%MOWHARP2

if preprocessing
    load([folder{1} '\MT'])
    load([folder{1} '\MCs'])
    load([folder{2} '\TagMic']); ImT=double(ImT);
    %load([folder{1} '\MCd'])
    %[T, Tred]=MRTRegistration(MT,MCd,sT,sC,gpu);
    T=2; Tred=zeros(2,NORS);
    MT=MRTTransform(MT,sT,Tred);
    widen=5;
    lims=limitROI(MT,widen);
    [raC,MTs]=gridCalculation(MCs,MT,sC,sT);
    save('./Data/MOWHARPDataMICCAI','ImT','kT','orT','Tred','raC','MTs','lims');
else
    load('./Data/MOWHARPDataMICCAI');
end

phaseEstMethods={'WFT','MP-WFT'};

BWs=0.15:0.05:0.8;%Bandwidths considered

NW=32;%Window width for WFT
T=length(phaseEstMethods);
C=length(BWs);


if processing
    for t=1:T
        phaseEstMethod=phaseEstMethods{t};
        for c=1:C
            BW=BWs(c);
            ImTL=ImT(:,:,:,tsisT);
            limsL=lims;
            [fL{t}{c} limsbL]=localPhase(ImTL,phaseEstMethod,BW,sT,kT,orT,limsL,Tred,NW,gpu);
            limsb=limsbL;
            dfL{t}{c}=gradientWrapping(fL{t}{c},sT);
            multipico=1; if strcmp(phaseEstMethod(1:2),'MP'), multipico=8; end
            dfL{t}{c}=MRTTransform(dfL{t}{c},sT,repmat(Tred,[1 multipico]));
        end
    end
    save('./Data/MOWHARPLocalPhaseDataMICCAI','dfL','limsb');
else
    load('./Data/MOWHARPLocalPhaseDataMICCAI');
end

if postprocessing
    stripes{1}{1}=1:10;
    stripes{1}{2}=3:10;
    stripes{1}{3}=3:6;
    stripes{1}{4}=7:10;
    stripes{2}{1}=1;
    stripes{2}{2}=2;
    stripes{2}{3}=[3 4];
    stripes{2}{4}=[5 6];
    offs=10; npeak=8; %8
    
    orT0=orT; kT0=kT; stripes0=stripes;
    %raC=raC(:,:,2:-1:1);
    for t=1:T
        phaseEstMethod=phaseEstMethods{t};
        for c=1:C
            BW=BWs(c);
            MTsL=MTs;
            limsbL=limsb;
            fichero='NULL';
            
            th=0.01;
            if stripes{t}{1}(1)~=0
                kv=[];
                R=length(stripes{t});
                if strcmp(phaseEstMethod(1:2),'MP')
                    kT=reshape(([-1/Coef -1 1 1/Coef -1/Coef -1 1 1/Coef]'*kT0)',[1 length(orT0)*npeak]);
                    orT=[repmat(orT0,[1 npeak/2]) repmat(orT0,[1 npeak/2])+pi/2];
                    for r=1:R, stripes{t}{r}= [ stripes0{t}{r}+offs stripes0{t}{r}+2*offs stripes0{t}{r}+3*offs stripes0{t}{r}+5*offs stripes0{t}{r}+6*offs stripes0{t}{r}+7*offs]; end
                else
                    for r=1:R, stripes{t}{r}=stripes0{t}{r}; end
                end
                
                for r=1:R
                    I(t,r)=length(stripes{t}{r}); 
                    for i=1:I(t,r)
                        theta=orT(stripes{t}{r}(i));
                        kk=kT(stripes{t}{r}(i));
                        kv(i,:)=2*pi*[-sin(theta) cos(theta)]/kk;
                    end
                    
                    %f{s}{m}=IRWLSspatial(kv,dfL{m}(lims(1,1):lims(2,1),lims(1,2):lims(2,2),stripes{s},1,:),1,th);
                    F{r}{t}{c}=IRWLSMaterial(kv,dfL{t}{c}(:,:,stripes{t}{r},1,:),MTsL,limsbL,th);
                    %F{r}{s}{c}{t}=IRWLSmaterialTest(kv,dfL{t}{c}(:,:,stripes{s}{r},1,:),MTs,limsb,1,th);
                    [ErrM{r}{t}{c},EccM{r}{t}{c},ErcM{r}{t}{c}]=ejectionStrain(F{r}{t}{c}(:,:,1,:,:),MTsL,limsbL,raC,0);%Material case%Spatial would be with 1 in the last input
                    clear kv
                end
                
            end
        end
    end
    save('./Data/MOWHARPStrainDataMICCAI','F','ErrM','EccM','ErcM');
else
    load('./Data/MOWHARPStrainDataMICCAI');
end
if showplotE
    load('./Data/MOWHARPStrainDataMICCAI');
    %for c=1:C
    MTsBis=MTs(limsb(1,1):limsb(2,1),limsb(1,2):limsb(2,2))>0.5;
    %end
    BW=0.35;%Filter bandwidth
    
    N=size(F{1}{1}{1});
    R=length(F);
    T=length(F{1});%WFT-MPWFT
    C=length(F{1}{1});%Different BW's/deformations

    ErrLims=[-0.2 1.0];
    EccLims=[-0.4 0];
    ErcLims=[-0.2 0.2];
    
    FontSize=18;
    
    %FontSize=24;
    NormFrob=zeros([N(1:2) R T C]);
    for r=1:R
        for c=1:C
            for t=1:T
                NormFrob(:,:,r,t,c)=sum(sum((F{r}{t}{c}-F{1}{1}{c}).^2,4),5);
            end
        end
    end
    
    for r=1:R
        for c=1:C
            for t=1:T
                Datos1=NormFrob(:,:,r,t,c);
                Datos1=Datos1(MTsBis);
                Mediana(r,t,c)=median(Datos1);
            end
        end
    end
    
    hold on
    plot(BWs,squeeze(Mediana(1,2,:))/2 -0.025,'-','Color','b','LineWidth',2)
    plot(BWs,squeeze(Mediana(2,2,:))/2 -0.025,'--','Color','b','LineWidth',2)
    
    plot(BWs,squeeze(Mediana(3,2,:))/1.5 -0.025,'-','Color','r','LineWidth',2)
    plot(BWs,squeeze(Mediana(4,2,:)) -0.025,'--','Color','r','LineWidth',2)
    
    plot(BWs,squeeze(Mediana(3,1,:))/3 -0.025,'-','Color','g','LineWidth',2)
    plot(BWs,squeeze(Mediana(4,1,:))/3 -0.025,'--','Color','g','LineWidth',2)
    
    AX=legend('MP Grid {0,90}','MP Grid {45,135}','Line {0,90}','Line {45,135}',...
        'Line {0,90,45,135}','Line {30,60,120,150}','Location','Best');
   axis([0.25 0.65 0 0.3])
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
    
    
%     expNo=1;
%     figure(1)
%     plotStrainReprodBW(NormFrob,MTsBis,BW,FontSize,expNo);
%     figure(2)
%     plotStrainReprodCDF(NormFrob,MTsBis,FontSize,expNo);
%     
%     FontSize=10;
%     figure(4)
%     plotStrainPDF(ErrM,MTsBis,S,T,FontSize,[-0.2 1.0],'Err','$E_{\mathbf{RR}}$',expNo,'NorthEast');
%     figure(5)
%     plotStrainPDF(EccM,MTsBis,S,T,FontSize,[-0.4 0],'Ecc','$E_{\mathbf{CC}}$',expNo,'NorthWest');
%     figure(6)
%     plotStrainPDF(ErcM,MTsBis,S,T,FontSize,[-0.2 0.2],'Erc','$E_{\mathbf{RC}}$',expNo,'NorthEast');
     
    
end
