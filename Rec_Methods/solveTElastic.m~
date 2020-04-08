function [rec,EH]=solveTElastic(rec)
%SOLVET   Estimates the rigid transforms each shot has been subject to

NT=size(rec.P);NT(end+1:5)=1;
Nn=prod(rec.NX(1:3))*rec.Nco;

permP=[6 7 8 9 5 1 2 3 4]; %take care of dimensions
rec.Pprev=rec.P; 

EH=[]; 
%Iterations
for n=1:rec.parT.nT    
    rec.parT.w(:)=1;%=rec.parT.w/multA;%We divide to try to jump quicker   
    flagw=zeros([1 1 1 1 rec.Nsh]);%Nothing has converged
    
    %RESIDUALS
    [r,xT,tR]=residualsElastic(rec); 
    
    %ENERGY
    Eprev=multDimSum(abs(r).^2,1:4)/Nn;
    if n==1;EH=Eprev;else EH=cat(1,EH,Eprev);end
    
    dt=zeros(size(rec.P),'like',rec.P);
    for indice=1:rec.tS.N %memory issues shot-by-shot
        if rec.parT.conv(:,:,:,:,indice)
        %JACOBIAN flag jacobian >=1
        [~,J,~]=gradientBSpline(rec,xT,tR,r,indice); 
        %GRADIENT DIRECTION
        dt(:,:,:,:,indice)=2*permute(multDimSum(real(bsxfun(@times,J,r(:,:,:,:,indice))),1:4),permP)/Nn;
        end
    end
    
    %BOTH JACOBIAN & GRADIENT DIRECTION flag jacobian =0
    %dt=permute(gradientBSpline(rec,xT,tR,r,0),permP)*2/Nn;  
    %sum(abs(dt(:)-dt2(:))), pause
    clear xT tR J
        
    if n==1;dtm1=abs(max(dt(:)))/rec.parT.maxD;end
    %May be attributed a physical meaning: maximum displacement at first iteration, useful to normalize
    dt=dt/dtm1;
    
    %SEARCH DIRECTION
    if strcmp(rec.parT.solver,'GD')
        st=dt; %Gradient descent
    elseif strcmp(rec.parT.solver,'CG')
        if n==1
            st=dt;
        else
            if rec.parT.CGstyle==1 %conjugate gradient separated
                for s=1:NT(5)
                    dts=dynInd(dt,s,5);dtps=dynInd(dtp,s,5);stps=dynInd(stp,s,5);
                    if strcmp(rec.parT.CGrule,'FR');be=(dts(:)'*dts(:))/((dtps(:)'*dtps(:))+eps);
                    elseif strcmp(rec.parT.CGrule,'PR');be=(dts(:)'*(dts(:)-dtps(:)))/((dtps(:)'*dtps(:))+eps);
                    elseif strcmp(rec.parT.CGrule,'HS');be=(dts(:)'*(dts(:)-dtps(:)))/((stps(:)'*(dts(:)-dtps(:)))+eps);
                    elseif strcmp(rec.parT.CGrule,'DY');be=(dts(:)'*dts(:))/((stps(:)'*(dts(:)-dtps(:)))+eps);
                    else error('CG rule %s not contemplated',rec.parT.CGrule);
                    end 
                    be=max(be,0);  
                    be(flagwd(s))=0;
                    st=dynInd(st,s,5,dts+be*stps);
                end
            else %conjugate gradient joint
                if strcmp(rec.parT.CGrule,'FR');be=(dt(:)'*dt(:))/((dtp(:)'*dtp(:))+eps);
                elseif strcmp(rec.parT.CGrule,'PR');be=(dt(:)'*(dt(:)-dtp(:)))/((dtp(:)'*dtp(:))+eps);
                elseif strcmp(rec.parT.CGrule,'HS');be=(dt(:)'*(dt(:)-dtp(:)))/((stp(:)'*(dt(:)-dtp(:)))+eps);
                elseif strcmp(rec.parT.CGrule,'DY');be=(dt(:)'*dt(:))/((stp(:)'*(dt(:)-dtp(:)))+eps);
                else error('CG rule %s not contemplated',rec.parT.CGrule);
                end
                be=max(be,0);
                st=dt+be*stp;
            end
        end
        dtp=dt;stp=st;
    else
        error('Undefined %s solver',rec.parT.solver);
    end
    
        %LINE SEARCH
    P=rec.P;A=rec.A;
    E=Eprev; E0=E;
    while any(flagw(:)<2)  
        rec.P=P-bsxfun(@times,st,rec.parT.w); %salto en w
        rec.P=dynInd(rec.P,flagw(:)~=2,5); % cojo solo los shots que no conveergen   
        rec.A=dynInd(A,flagw(:)~=2,5);    
        r=residualsElastic(rec);       %residuo para backtracking
        E=dynInd(E,flagw(:)~=2,5,multDimSum(abs(r).^2,1:4)/Nn); %energia
        flagw(E>Eprev & flagw==1)=2;
        flagw(E<=Eprev & flagw==0)=1;
        Eprev(flagw==1)=E(flagw==1);        
        rec.parT.w(flagw~=2)=rec.parT.w(flagw~=2)/rec.parT.mult;%We divide the jump 
        if all(rec.parT.w(flagw~=2)<rec.parT.tolW);break;end             
    end
    %[squeeze(E0) squeeze(E) squeeze(flagw)]
    %pause
    rec.P=P;rec.A=A;     
    flagup=ismember(flagw(:),1:2); % update only the good ones
    rec.parT.w(flagup)=rec.parT.w(flagup)*rec.parT.mult;
    %flagup=(flagw(:)==2);  
    if any(flagup); rec.P=dynInd(rec.P,flagup,5,dynInd(rec.P,flagup,5)-bsxfun(@times,dynInd(st,flagup,5),dynInd(rec.parT.w,flagup,5)));end
    
    %IMPLEMENT CONVERGENCE CRITERIA INTERNAL
    %SANTI, PODEMOS ASUMIR QUE ESTOS CASOS HAN CONVERGIDO, HA HABIDO UN CG
    %SIN ENCONTRAR UN DECRECIMIENTO DE ENERGIA, SE HA INTENTADO UN GD (PARA EL CASO POR SEPARADO) Y NI
    %CON ESAS, SE QUITAN ESTOS SHOTS Y SE COMPUTAN LOS PARAMETROS PARA EL
    %RESTO DE SHOTS
    if n>1 
    %    rec.parT.w<rec.parT.tolW & flagwd
        T_old=rec.T;
        rec.T=transformBSpline(rec.T0,rec.tS,rec.P); %Compare deformation fields between iterations
        conv=ConvField(T_old,rec.T,rec.M);
        if max(conv)<rec.parT.toler
            rec.parT.out=1;
            break,
        else
            rec.parT.out=0;    
        end
        clear T_old
    else
        rec.T=transformBSpline(rec.T0,rec.tS,rec.P);
    end
    flagwd=(rec.parT.w<rec.parT.tolW);%Cases on which a decrease was not found
    %IMPLEMENT MEAN SUBTRACTION TO REMOVE DRIFTING
    if rec.parT.meanT
        %caso Elastico remove mean
        PuntosMed=multDimMea(rec.P,5);
        rec.P=bsxfun(@minus,rec.P,PuntosMed); 
        TBmed=transformBSpline(rec.T0,rec.tS,PuntosMed); %Transformo la imagen a -Media
        rec.x=applyTransform(rec.x,TBmed,rec.NDimen); clear TBmed PuntosMed
        rec.x=rec.x.*rec.M;
    end
end


