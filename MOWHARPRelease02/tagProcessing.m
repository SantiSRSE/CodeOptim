addpath(genpath('.'));

%Experiment 1: reproducibility of increased number of orientations
%Experiment 2: reproducibility of repeated measurements versus extended
%orientations
%Experiment 3: visual results
%Experiment 4: interferences
%Experiment 5: error with no interferences
expNo=1;
%Coef=2;

if expNo==4
    NORS=[2 3 4 5 6 7 8 9 10 12 14 15 16 18 20 21 24 25 27 28 30 32 35 36];
else
    NORS=36;
end
for nNOR=1:length(NORS)
    NOR=NORS(nNOR);
    folder{1}='.\Data\EC-C';%Path to MR-C
    folder{2}='.\Data\EC-T';%Path to MR-T

    gpu=0;
    fftw('planner','measure');

    preprocessing=1;
    processing=1;
    postprocessing=1;
    showplotE=1;

    if expNo>3
        realData=0;
    else
        realData=1;
    end    
    if realData
        sC=[1.25 1.25];%MR-C resolution
        sT=[4/3 4/3];%MR-T resolution
    else
        sC=[1 1];
        sT=[1 1];
    end

    %tsisC=6;%MOWHARP
    tsisC=9;%MOWHARP2
    TagIm=[1:15 45 17:24 46 26:30 47 32:33 48 35:44];%Case MOWHARP2

    if preprocessing
        if realData
            [MCd,MCs,MT]=readOHARPSegmentations(folder,tsisC,TagIm);
            [~,ImT,eT,kT,orT]=readOHARPImages(folder,TagIm);    
            [T Tred]=MRTRegistration(MT,MCd,sT,sC,gpu);
            MT=MRTTransform(MT,sT,Tred);
            widen=5;
            lims=limitROI(MT,widen);
            [raC,MTs]=gridCalculation(MCs,MT,sC,sT);
            save('./Data/MOWHARPData','ImT','kT','orT','Tred','raC','MTs','lims');
        else
            [ImT,kT,orT,Tred,raC,MTs,lims,ErrMGT,EccMGT,ErcMGT,lambdas]=generatePhantomHARP(NOR);
            save('./Data/MOWHARPSyntC228','ImT','kT','orT','Tred','raC','MTs','lims','ErrMGT','EccMGT','ErcMGT','lambdas');
        end
    else
        if realData
            load('./Data/MOWHARPData');
        else
            load('./Data/MOWHARPSyntC110');
        end
    end

    if expNo==4
        phaseEstMethods={'FT','SFT','WFT','WSFT','MP-FT','MP-WFT'}; %,'MPe'
        %{'WFT','WSFT'};%{'FT','SFT'};%{'FT','SFT','WFT','WSFT'};%Phase estimation methods considered
        %phaseEstMethods={'WFT','AWFT'};
    else
        if realData, phaseEstMethods={'FT','WFT','MP-FT','MP-WFT'};
        else phaseEstMethods={'FT','WFT','MP-FT','MP-WFT'}; end
    end
    if expNo==5
        BWs=0.35; %0.5 para AWFT
    else
        BWs=0.15:0.05:0.8;%Bandwidths considered
    end
    NW=32;%Window width for WFT
    T=length(phaseEstMethods);
    if expNo==5
        C=length(lambdas);
    else
        C=length(BWs);
    end
        
    if processing
        if realData
            %for s=1:size(ImT,4)
            %    figure(s)
            %    imshow(ImT(:,:,1,s),[])
            %end
            %pause
            %This parameter has been fixed by the lines that appear commented just
            %before
            tsisT=6;
        elseif expNo==5
            tsisT=1;
        elseif expNo==4
            tsisT=10;
        end
        fL=cell(T,C);
        dfL=cell(T,C);        
        for t=1:T
            phaseEstMethod=phaseEstMethods{t};             
            for c=1:C
                if expNo==5
                    BW=BWs;
                    ImTL=ImT{c}(:,:,:,tsisT);
                    limsL=lims{c};
                elseif expNo==4
                    BW=BWs(c);
                    ImTL=ImT{tsisT};
                    limsL=lims{tsisT};
                else
                    BW=BWs(c);
                    ImTL=ImT(:,:,:,tsisT);
                    limsL=lims;
                end      
                [fL{t}{c} limsbL]=localPhase(ImTL,phaseEstMethod,BW,sT,kT,orT,limsL,Tred,NW,gpu);
                if expNo==5;
                    limsb{c}=limsbL;
                else
                    limsb=limsbL;
                end
                dfL{t}{c}=gradientWrapping(fL{t}{c},sT);
                multipico=1; if strcmp(phaseEstMethod(1:2),'MP'), multipico=8; end
                dfL{t}{c}=MRTTransform(dfL{t}{c},sT,repmat(Tred,[1 multipico]));                    
            end
        end
        if realData, save('./Data/MOWHARPLocalPhaseData3','dfL','limsb');
        else save('./Data/MOWHARPLocalPhaseDataSyntC228','dfL','limsb'); end
    else
        if realData, load('./Data/MOWHARPLocalPhaseData');
        else load('./Data/MOWHARPLocalPhaseDataSyntC110'); end
    end

    %figure(3);for s=1:44;subtightplot(4,11,s,[0 0]);imshow(dfL{1}{6}(:,:,s,1,1),[]);end
    %figure(2);for s=1:44;subtightplot(4,11,s,[0 0]);imshow(dfL{2}{6}(:,:,s,1,1),[]);end
    %return

if postprocessing
    switch expNo
        %-(orT*180/pi-180) 
        case 1,
%             stripes{1}{1}=[1 19];%0/90
%             stripes{1}{2}=[28 10];%-45/45
%             stripes{2}{1}=[25 1 13];%-60/0/60
%             stripes{2}{2}=[31 7 19];%-30/30/90
%             stripes{3}{1}=[25 31 1 7 13 19];%-60/-30/0/30/60/90
%             stripes{3}{2}=[22 28 34 4 10 16];%-75/-45/-15/15/45/75
%             stripes{4}{1}=[23 27 31 35 3 7 11 15 19];%-70/-50/-30/-10/10/30/50/70/90
%             stripes{4}{2}=[21 25 29 33 1 5 9 13 17];%-80/-60/-40/-20/0/20/40/60/80
%             stripes{5}{1}=[21 23 25 27 29 31 33 35 1 3 5 7 9 11 13 15 17 19];%-80/-70/-60/-50/-40/-30/-20/-10/0/10/20/30/40/50/60/70/80/90
%             stripes{5}{2}=[20 22 24 26 28 30 32 34 36 2 4 6 8 10 12 14 16 18];%%-85/-75/-65/-55/-45/-35/-25/-15/-5/5/15/25/35/45/55/65/75/85
            stripes{1}{1}=[1 19];%0/90
            stripes{1}{2}=[28 10];%-45/45
            stripes{2}{1}=[1 10 19 28];% 0/45/90-45 
            stripes{2}{2}=[4 13 34 25];%15/60/-15/-60
            stripes{3}{1}=[23 27 31 1 7 11 15 19];%-70/-50/-30/0/30/50/70/90
            stripes{3}{2}=[21 25 29 34 4 9 13 17];%-80/-60/-40/-15/15/40/60/80
            stripes{4}{1}=[23 27 31 35 3 7 11 15 19];%-70/-50/-30/-10/10/30/50/70/90
            stripes{4}{2}=[21 25 29 33 1 5 9 13 17];%-80/-60/-40/-20/0/20/40/60/80
            stripes{5}{1}=[21 23 25 27 29 31 33 35 1 3 5 7 9 11 13 15 17 19];%-80/-70/-60/-50/-40/-30/-20/-10/0/10/20/30/40/50/60/70/80/90
            stripes{5}{2}=[20 22 24 26 28 30 32 34 36 2 4 6 8 10 12 14 16 18];%%-85/-75/-65/-55/-45/-35/-25/-15/-5/5/15/25/35/45/55/65/75/85
            offs=44; npeak=4;
        case 2,
            stripes{1}{1}=[1 37 38 19 39 40];%0/0/0/90/90/90
            stripes{1}{2}=[28 41 42 10 43 44];%-45/-45/-45/45/45/45
            stripes{2}{1}=[25 31 1 7 13 19];%-60/-30/0/30/60/90
            stripes{2}{2}=[22 28 34 4 10 16];%-75/-45/-15/15/45/75
            offs=44; npeak=4;
        case 3,
            stripes{1}{1}=[1 19];%0/90
            stripes{2}{1}=[25 1 13];%-60/0/60
            stripes{3}{1}=[28 1 10 19];%-45/0/45/90
            stripes{4}{1}=[25 31 1 7 13 19];%-60/-30/0/30/60/90
            stripes{5}{1}=[21 25 29 33 1 5 9 13 17];%-80/-60/-40/-20/0/20/40/60/80
            stripes{6}{1}=[22 25 28 31 34 1 4 7 10 13 16 19];%-75/-60/-45/-30/-15/0/15/30/45/60/75/90
            stripes{7}{1}=[21 23 25 27 29 31 33 35 1 3 5 7 9 11 13 15 17 19];%-80/-70/-60/-50/-40/-30/-20/-10/0/10/20/30/40/50/60/70/80/90
            stripes{8}{1}=[20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19];%-80/-75/-70/-65/-60/-55/-50/-45/-40/-35/-30/-25/-20/-15/-10/-5/0/5/10/15/20/25/30/35/40/45/50/55/60/65/70/75/80/85/90
            offs=44; npeak=4;
        case 4,            
            %stripes{1}{1}=1:size(dfL{1}{1},3);%0/90 
            stripes{1}{1}=[1 19];% 1 19 y 28 10
            offs=36; npeak=4;
        case 5,            
            stripes{1}{1}=[18 36];
            %stripes{2}{1}=[1 19 28 10];
            %stripes{3}{1}=[21 25 29 34 4 9 13 17];
            stripes{2}{1}=1:36;
            offs=36; npeak=4; %8
    end
    orT0=orT; kT0=kT; stripes0=stripes;
    %raC=raC(:,:,2:-1:1);
    for t=1:T
        phaseEstMethod=phaseEstMethods{t};
        for c=1:C
            if expNo==5
                BW=BWs;
                MTsL=MTs{c};
                limsbL=limsb{c};
            elseif expNo==4
                BW=BWs(c); 
                tsisT=10;
                MTsL=MTs{tsisT};       
                limsbL=limsb;
            else
                BW=BWs(c);
                MTsL=MTs;       
                limsbL=limsb;
            end               
            fichero='NULL';               
                  
            th=0.01;         
            S=length(stripes);
            for s=1:S
                if stripes{s}{1}(1)~=0
                    kv=[];
                    R=length(stripes{s});
                    if strcmp(phaseEstMethod(1:2),'MP')
                    if strcmp(phaseEstMethod(1:3),'MP-')
                        if npeak==6
                            kT=reshape(([-1/3 -1/2 -1 1 1/2 1/3]'*kT0)',[1 length(orT0)*npeak]); orT=repmat(orT0,[1 npeak]);
                        elseif npeak==8
                            kT=reshape(([-1/Coef -1 1 1/Coef -1/Coef -1 1 1/Coef]'*kT0)',[1 length(orT0)*npeak]); 
                            orT=[repmat(orT0,[1 npeak/2]) repmat(orT0,[1 npeak/2])+pi/2];
                        elseif npeak==4
                            kT=reshape(([-1/Coef -1 1 1/Coef]'*kT0)',[1 length(orT0)*npeak]); orT=repmat(orT0,[1 npeak]);
                        end                        
                        %for r=1:R, stripes{s}{r}= [stripes0{s}{r} stripes0{s}{r}+offs stripes0{s}{r}+2*offs stripes0{s}{r}+3*offs stripes0{s}{r}+4*offs stripes0{s}{r}+5*offs]; end
                        
                        %for r=1:R, stripes{s}{r}= [ stripes0{s}{r}+2*offs stripes0{s}{r}+3*offs stripes0{s}{r}+6*offs stripes0{s}{r}+7*offs]; end
                        for r=1:R, stripes{s}{r}= [ stripes0{s}{r}+offs stripes0{s}{r}+2*offs stripes0{s}{r}+3*offs]; end
                    elseif strcmp(phaseEstMethod(1:3),'MPe')
                        %mp=[1 1 1 1]; mo=[0 pi/4 pi/2 3*pi/4];
                        %mp=[1 2.13 1 2.13]; mo=[0 pi/4 pi/2 3*pi/4];
                        %mp=[1 2 1 2]; mo=[0 0 pi/2 pi/2]; 
                        %kT=reshape(((1./mp)'*kT0)',[1 length(orT0)*npeak]); 
                        %orT=[orT0+mp(1) orT0+mp(2) orT0+mp(3) orT0+mp(4)];                                                
                        %for r=1:R, stripes{s}{r}= [stripes0{s}{r} stripes0{s}{r}+offs stripes0{s}{r}+2*offs stripes0{s}{r}+3*offs]; end
                    end
                    end
                    I(s)=length(stripes{s}{1});
                    for r=1:R
                        for i=1:I(s)
                            theta=orT(stripes{s}{r}(i));
                            kk=kT(stripes{s}{r}(i)); 
                            kv(i,:)=2*pi*[-sin(theta) cos(theta)]/kk;
                        end                           
                        %f{s}{m}=IRWLSspatial(kv,dfL{m}(lims(1,1):lims(2,1),lims(1,2):lims(2,2),stripes{s},1,:),1,th);
                        F{r}{s}{c}{t}=IRWLSMaterial(kv,dfL{t}{c}(:,:,stripes{s}{r},1,:),MTsL,limsbL,th);                      
                        %F{r}{s}{c}{t}=IRWLSmaterialTest(kv,dfL{t}{c}(:,:,stripes{s}{r},1,:),MTs,limsb,1,th);     
                        [ErrM{r}{s}{c}{t},EccM{r}{s}{c}{t},ErcM{r}{s}{c}{t}]=ejectionStrain(F{r}{s}{c}{t}(:,:,1,:,:),MTsL,limsbL,raC,0);%Material case%Spatial would be with 1 in the last input                        
                    end
                %else
                %    R=length(stripes{s});                                    
                %    for r=1:R   
                %        F{r}{s}{c}{t}=0;                      
                %        ErrM{r}{s}{c}{t}=ErrMGT(limsb(1,1):limsb(2,1),limsb(1,2):limsb(2,2));
                %        EccM{r}{s}{c}{t}=EccMGT(limsb(1,1):limsb(2,1),limsb(1,2):limsb(2,2));
                %        ErcM{r}{s}{c}{t}=ErcMGT(limsb(1,1):limsb(2,1),limsb(1,2):limsb(2,2));
                %    end
                end                       
            end           
        end   
    end    
    if ~realData                        
        t=length(phaseEstMethods)+1;
        S=length(stripes);
        for s=1:S
            R=length(stripes{s});                                    
            for c=1:C
                for r=1:R   
                    F{r}{1}{c}{t}=0;                      
                    if expNo==4
                        ErrM{r}{s}{c}{t}=ErrMGT{tsisT}(limsb(1,1):limsb(2,1),limsb(1,2):limsb(2,2));
                        EccM{r}{s}{c}{t}=EccMGT{tsisT}(limsb(1,1):limsb(2,1),limsb(1,2):limsb(2,2));
                        ErcM{r}{s}{c}{t}=ErcMGT{tsisT}(limsb(1,1):limsb(2,1),limsb(1,2):limsb(2,2));
                    else
                        ErrM{r}{s}{c}{t}=ErrMGT{c}(limsb{c}(1,1):limsb{c}(2,1),limsb{c}(1,2):limsb{c}(2,2));
                        EccM{r}{s}{c}{t}=EccMGT{c}(limsb{c}(1,1):limsb{c}(2,1),limsb{c}(1,2):limsb{c}(2,2));
                        ErcM{r}{s}{c}{t}=ErcMGT{c}(limsb{c}(1,1):limsb{c}(2,1),limsb{c}(1,2):limsb{c}(2,2));
                    end
                end
            end
        end
    end 
    %save('./Data/MOWHARPStrainDataC228','F','ErrM','EccM','ErcM');
end
if showplotE
    %load('./Data/MOWHARPStrainDataC228');
    if expNo==4
        tsisT=10;
        MTsBis{tsisT}=MTs{tsisT}(limsb(1,1):limsb(2,1),limsb(1,2):limsb(2,2))>0.5;
        BW=0.5;%Filter bandwidth
    elseif expNo~=5
        MTsBis=MTs(limsb(1,1):limsb(2,1),limsb(1,2):limsb(2,2))>0.5;
        BW=0.15:0.05:0.8;%Filter bandwidth
    else
        for c=1:C
            MTsBis{c}=MTs{c}(limsb{c}(1,1):limsb{c}(2,1),limsb{c}(1,2):limsb{c}(2,2))>0.5;
        end
        BW=0.35;%Filter bandwidth
    end
    N=size(F{1}{1}{1}{1});
    S=length(F{1});%Different I's    
    if expNo~=4
        I=[2 3 4 6 9 12 18 36];   
    elseif expNo==4
        %I={'FT','SFT','WFT','SWFT','GT'};
        I={'FT','SFT','WFT','WSFT','MP-FT','MP-AWFT','GT'};
        %I={'HARP','WHARP','GT'};  
    end
    C=length(F{1}{1});%Different BW's/deformations
    T=length(F{1}{1}{1});%FT-WFT
    ErrLims=[-0.2 1.0];
    EccLims=[-0.4 0];
    ErcLims=[-0.2 0.2];
    
    FontSize=18; 
    if expNo==3
        FontSize=24;
        figure(1)
        plotStrainIlus(ErrM,S,C,FontSize,ErrLims,BW,I,'ERR_HARP',1,realData);    
        figure(2)
        plotStrainIlus(ErrM,S,C,FontSize,ErrLims,BW,I,'ERR_WHARP',2,realData);
        figure(3)
        plotStrainIlus(ErrM,S,C,FontSize,ErrLims,BW,I,'ERR_MPHARP',3,realData);
        figure(4)
        plotStrainIlus(ErrM,S,C,FontSize,ErrLims,BW,I,'ERR_MPAWHARP',4,realData);
        
        figure(5)
        plotStrainIlus(EccM,S,C,FontSize,EccLims,BW,I,'ECC_HARP',1,realData);    
        figure(6)
        plotStrainIlus(EccM,S,C,FontSize,EccLims,BW,I,'ECC_WHARP',2,realData);
        figure(7)
        plotStrainIlus(EccM,S,C,FontSize,EccLims,BW,I,'ECC_MPHARP',3,realData);
        figure(8)
        plotStrainIlus(EccM,S,C,FontSize,EccLims,BW,I,'ECC_MPAWHARP',4,realData);
        
        figure(9)
        plotStrainIlus(ErcM,S,C,FontSize,ErcLims,BW,I,'ERC_HARP',1,realData);    
        figure(10)
        plotStrainIlus(ErcM,S,C,FontSize,ErcLims,BW,I,'ERC_WHARP',2,realData);  
        figure(11)
        plotStrainIlus(ErcM,S,C,FontSize,ErcLims,BW,I,'ERC_MPHARP',3,realData);
        figure(12)
        plotStrainIlus(ErcM,S,C,FontSize,ErcLims,BW,I,'ERC_MPAWHARP',4,realData);
    elseif expNo==4
        FontSize=12;
        figure(4)
        plotStrainIlus(ErrM,T,C,FontSize,ErrLims,BWs,I,'ERR',1,realData,length(stripes{1}{1})); 
        figure(5)
        plotStrainIlus(EccM,T,C,FontSize,EccLims,BWs,I,'ECC',1,realData,length(stripes{1}{1}));
        figure(6)
        plotStrainIlus(ErcM,T,C,FontSize,ErcLims,BWs,I,'ERC',1,realData,length(stripes{1}{1})); break
    elseif expNo==5
        for c=1:C%Different deformations
            noPoints=sum(MTsBis{c}(:));
            for s=1:S%Different stripe numbers
                for t=1:T-1%Different methods                    
                    ErrRR{c}(:,:,s,t)=ErrM{1}{s}{c}{t}-ErrM{1}{s}{c}{T};
                    ErrCC{c}(:,:,s,t)=EccM{1}{s}{c}{t}-EccM{1}{s}{c}{T};
                    ErrRC{c}(:,:,s,t)=ErcM{1}{s}{c}{t}-ErcM{1}{s}{c}{T};
                end
            end        
            ErrARR(c,:,:)=permute(sqrt(sum(sum(ErrRR{c}.^2,1),2)/noPoints),[5 3 4 1 2]);%Error
            ErrACC(c,:,:)=permute(sqrt(sum(sum(ErrCC{c}.^2,1),2)/noPoints),[5 3 4 1 2]);
            ErrARC(c,:,:)=permute(sqrt(sum(sum(ErrRC{c}.^2,1),2)/noPoints),[5 3 4 1 2]);
            BiasARR(c,:,:)=permute(sum(sum(ErrRR{c},1),2)/noPoints,[5 3 4 1 2]);%Error
            BiasACC(c,:,:)=permute(sum(sum(ErrCC{c},1),2)/noPoints,[5 3 4 1 2]);
            BiasARC(c,:,:)=permute(sum(sum(ErrRC{c},1),2)/noPoints,[5 3 4 1 2]);
            StdARR(c,:,:)=sqrt(ErrARR(c,:,:).^2-BiasARR(c,:,:).^2);
            StdACC(c,:,:)=sqrt(ErrACC(c,:,:).^2-BiasACC(c,:,:).^2);
            StdARC(c,:,:)=sqrt(ErrARC(c,:,:).^2-BiasARC(c,:,:).^2);
        end
        figure(1)
        plotStrainError(StdARR,lambdas,'ERR',1);
        figure(2)
        plotStrainError(StdACC,lambdas,'ECC',1);
        figure(3)
        plotStrainError(StdARC,lambdas,'ERC',1);                  
        figure(4)
        plotStrainError(BiasARR,lambdas,'ERR',0);
        figure(5)
        plotStrainError(BiasACC,lambdas,'ECC',0);
        figure(6)
        plotStrainError(BiasARC,lambdas,'ERC',0);
        figure(7)
        plotStrainError(ErrARR,lambdas,'ERR',1);
        figure(8)
        plotStrainError(ErrACC,lambdas,'ECC',1);
        figure(9)
        plotStrainError(ErrARC,lambdas,'ERC',1);
        
    
        

        
    else
        %FontSize=24;
        NormFrob=zeros([N(1:2) T C S]);
        for s=1:S
            for c=1:C
                for t=1:T
                    NormFrob(:,:,t,c,s)=sum(sum((F{1}{s}{c}{t}-F{2}{s}{c}{t}).^2,4),5);
                end
            end
        end
        for s=1:S
            for r=1:R
                for c=1:C
                    for t=1:T
                        Datos1=NormFrob(:,:,t,c,s);
                        Datos2=NormFrob(:,:,t,c,r);
                        Datos1=Datos1(MTsBis);
                        Datos2=Datos2(MTsBis);
                        Mediana(s,r,c,t)=100*(median(Datos2)-median(Datos1))/median(Datos2);
                        Difer(s,r,c,t)=ranksum(Datos1,Datos2)+j*(median(Datos2)-median(Datos1));
                    end
                end
            end
        end     
                        
        figure(1)
        plotStrainReprodBW(NormFrob,MTsBis,BW,FontSize,expNo);    
        figure(2)
        plotStrainReprodCDF(NormFrob,MTsBis,FontSize,expNo);   
    
        FontSize=10;
        figure(4)
        plotStrainPDF(ErrM,MTsBis,S,T,FontSize,[-0.2 1.0],'Err','$E_{\mathbf{RR}}$',expNo,'NorthEast');       
        figure(5)
        plotStrainPDF(EccM,MTsBis,S,T,FontSize,[-0.4 0],'Ecc','$E_{\mathbf{CC}}$',expNo,'NorthWest');
        figure(6)
        plotStrainPDF(ErcM,MTsBis,S,T,FontSize,[-0.2 0.2],'Erc','$E_{\mathbf{RC}}$',expNo,'NorthEast');        

        Mediana(:,:,2:2:end,2)

        real(Difer(:,:,2:2:end,2))    
        imag(Difer(:,:,2:2:end,2))
        real(permute(Difer(2,1,2:2:end,2),[1 3 2]))
    end
%end
    if expNo==4
        close all
        clearvars -EXCEPT NORS nNOR expNo
    end
end
end