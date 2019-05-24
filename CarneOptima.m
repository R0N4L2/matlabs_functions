clear;
clc;
poolobj = parpool('local');
pA=xlsread('C:\Users\ronal\OneDrive\Documentos\carnes\optimizacion_carnes2.xls','piernaA');
pB=xlsread('C:\Users\ronal\OneDrive\Documentos\carnes\optimizacion_carnes2.xls','piernaB');
bA=xlsread('C:\Users\ronal\OneDrive\Documentos\carnes\optimizacion_carnes2.xls','brazoA');
[pes,~,peso]=xlsread('C:\Users\ronal\OneDrive\Documentos\carnes\optimizacion_carnes2.xls','pesaje');
[dem,~,demstr]=xlsread('C:\Users\ronal\OneDrive\Documentos\carnes\optimizacion_carnes2.xls','Demanda');
[fo,~,fostr]=xlsread('C:\Users\ronal\OneDrive\Documentos\carnes\optimizacion_carnes2.xls','NivelServicio');
fo=fo(:,[3:4,6]);
dem=dem(:,3:4);
dem=sumcond(dem(:,[2,1]));
fr=promcond(fo(:,[3,1]));
dem=dem(:,[2,1]);
pbA=pes(1);
ppA=pes(2);
ppB=pes(3);
[ti1,MpA,rel1,C1,E1]=matriz_generador(pA,'pA',ppA);
[ti2,MpB,rel2,C2,E2]=matriz_generador(pB,'pB',ppB);
[ti3,MbA,rel3,C3,E3]=matriz_generador(bA,'bA',pbA);
[s0,d0]=size(MpA);
[s1,d1]=size(MpB);
[s2,d2]=size(MbA);
Meq=[[MpA,zeros(s0,d1+d2)];[zeros(s1,d0),MpB,zeros(s1,d2)];[zeros(s2,d0+d1),MbA]];
EE=[[E1,zeros(s0-1,d1+d2)];[zeros(s1-1,d0),E2,zeros(s1-1,d2)];[zeros(s2-1,d0+d1),E3]];
beq=zeros(s0+s1+s2,1);
tpA=find(ti1<0);
tpB=find(ti2<0)+d0;
tbA=find(ti3<0)+d0+d1;
tipo=[ti1,ti2,ti3];
f=zeros(length(tipo),1);
fp=promponcond(fo(:,[3,2]),fo(:,[3,1]));
[la,idx]=ismember(tipo',fp(:,1));
f(la)=fp(idx(la),2)+1;
f(find(ti1<0))=-mean(sum(C1(2:end,:)>0,2))-10;
f(find(ti2<0)+d0)=-mean(sum(C2(2:end,:)>0,2))-5;
f(find(ti3<0)+d0+d1)=-mean(sum(C3(2:end,:)>0,2))-1;
if any(dem(:,1)>0)
    demp=dem(dem(:,1)>0,:);
    [la,idx]=ismember(demp(:,2),fr(:,1));
    Mi=sparse(size(demp,1),length(tipo));
    fi=fr(idx(la),2);
    bi=demp(:,1);
    for j=1:size(demp,1)
        Mi(j,ismember(tipo,demp(j,2)))=sum(unique(Meq(:,ismember(tipo,demp(j,2))),'rows'),1);
    end
else
    demp=[];
    Mi=[];
    bi=[];
    fi=[];
end
if any(dem(:,1)==0)
    dem0=dem(dem(:,1)==0,:);
    Meq0=sparse(size(dem0,1),length(tipo));
    beq0=dem0(:,1);
    for j=1:size(dem0,1)
        Meq0(j,ismember(tipo,dem0(j,2)))=sum(unique(Meq(:,ismember(tipo,dem0(j,2))),'rows'),1);
    end
else
    dem0=[];
    Meq0=[];
end 
lb=zeros(d0+d1+d2,1);
%limpieza de memoria haciendo matrices sparse:
Mi=sparse(Mi);
Meq=sparse(Meq);
Meq0=sparse(Meq0);
delete(poolobj);
ub=1e3*ones(d0+d1+d2,1);
[la,idx]=ismember(tipo',dem(:,2));
ub(la)=1e3*(dem(idx(la),1)>0);
%ub(ub==0)=.08*sum(bi)/sum(ub==0);
%.9 de cuplimiento (oferta minima) y 1.08 de oferta maxima
% Setup Options
opts = optiset('solver','clp'); %CLP is a row constraint solver
if ~isempty(demp)&&~isempty(dem0)
    A=[Mi;sum([Mi;Meq0],1);Meq;Meq0;EE];
    ru=[1.08*bi;1.08*sum(bi);zeros(size([Meq;Meq0;EE],1),1)];
    rl=[bi.*fi;.9*sum(bi);zeros(size([Meq;Meq0;EE],1),1)];
elseif isempty(dem0)
    A=[Mi;sum(Mi,1);Meq;EE];
    ru=[inf(size(bi,1),1);1.08*sum(bi);zeros(size([Meq;EE],1),1)];
    rl=[bi.*fi;.9*sum(bi);zeros(size([Meq;EE],1),1)];
else
    A=[Meq;Meq0;EE];
    ru=zeros(size([Meq;Meq0;EE],1),1);
    rl=ru;
end
% Build OPTI Object
Opt = opti('f',-f,'lin',A,rl,ru,'bounds',lb,ub,'options',opts);
% Solve Problem
[x0,fval] = solve(Opt);
%variables innecesarias
varsi=find(x0<1e3&(tipo'>0|x0>0));
tipo=tipo(:,varsi);
Mi=Mi(:,varsi);
w1=find(any(Mi>0,2));
bi=bi(w1,:);
fi=fi(w1,:);
Mi=Mi(w1,:);
Meq0=Meq0(:,varsi);
Meq=Meq(:,varsi);
w=find(any(Meq>0,2));
Meq=Meq(w,:);
EE=EE(:,varsi);
EE=EE(any(EE>0,2),:);
f=f(varsi,:);
lb=lb(varsi,:);
ub=ub(varsi,:);
ti1=ti1(intersect(1:size(ti1,2),varsi));
rel1=rel1(intersect(1:size(rel1,2),varsi));
ti2=ti2(intersect((1:size(ti2,2))+d0,varsi)-d0);
rel2=rel2(intersect((1:size(rel2,2))+d0,varsi)-d0);
ti3=ti3(intersect((1:size(ti3,2))+d0+d1,varsi)-d0-d1);
rel3=rel3(intersect((1:size(rel3,2))+d0+d1,varsi)-d0-d1);
%recalcular problema eliminando las variables innecesarias
% Setup Options
opts = optiset('solver','clp'); %CLP is a row constraint solver
if ~isempty(Mi)&&~isempty(Meq0)
    A=[Mi;sum([Mi;Meq0],1);Meq;Meq0;EE];
    ru=[1.08*bi;1.08*sum(bi);zeros(size([Meq;Meq0;EE],1),1)];
    rl=[bi.*fi;.9*sum(bi);zeros(size([Meq;Meq0;EE],1),1)];
elseif isempty(Meq0)
    A=[Mi;sum(Mi,1);Meq;EE];
    ru=[inf(size(bi,1),1);1.08*sum(bi);zeros(size([Meq;EE],1),1)];
    rl=[bi.*fi;.9*sum(bi);zeros(size([Meq;EE],1),1)];
else
    A=[Meq;Meq0;EE];
    ru=zeros(size([Meq;Meq0;EE],1),1);
    rl=ru;
end
% Build OPTI Object
Opt = opti('f',-f,'lin',A,rl,ru,'bounds',lb,ub,'options',opts);
% Solve Problem
[x0,fval] = solve(Opt);
intsol=-fval;
%configuracion 0
r1=size(C1,1)-1;
Conf0pA=C1([1;intersect(1:r1,w)+1],:);
Conf0pA=[Conf0pA,[0;x0(find(ti1<0))]];
Conf0pA=Conf0pA([1;find(Conf0pA(2:end,end)>0)+1],:);
Conf0pA=full(Conf0pA(:,any(Conf0pA(2:end,:)>0,1)));
a=unique(Conf0pA(1,:),'sorted');
X=[];
for j=a
    X=[X,[j;sum(Conf0pA(2:end,Conf0pA(1,:)==j),2)]];
end
Conf0pA=X';
Conf0pA=num2cell(Conf0pA);
Conf0pA{1,1}='Cantidad de pierna A';
r2=size(C2,1)-1;
Conf0pB=C2([1;intersect((1:r2)+r1,w)-r1+1],:);
d0=size(ti1,2);
Conf0pB=[Conf0pB,[0;x0(find(ti2<0)+d0)]];
Conf0pB=Conf0pB([1;find(Conf0pB(2:end,end)>0)+1],:);
Conf0pB=full(Conf0pB(:,any(Conf0pB(2:end,:)>0,1)));
a=unique(Conf0pB(1,:),'sorted');
X=[];
for j=a
    X=[X,[j;sum(Conf0pB(2:end,Conf0pB(1,:)==j),2)]];
end
Conf0pB=X';
Conf0pB=num2cell(Conf0pB);
Conf0pB{1,1}='Cantidad de pierna B';
r3=size(C3,1)-1;
Conf0bA=C3([1;intersect((1:r3)+r1+r2,w)-r1-r2+1],:);
d1=size(ti2,2);
Conf0bA=[Conf0bA,[0;x0(find(ti3<0)+d0+d1)]];
Conf0bA=Conf0bA([1;find(Conf0bA(2:end,end)>0)+1],:);
Conf0bA=full(Conf0bA(:,any(Conf0bA(2:end,:)>0,1)));
a=unique(Conf0bA(1,:),'sorted');
X=[];
for j=a
    X=[X,[j;sum(Conf0bA(2:end,Conf0bA(1,:)==j),2)]];
end
Conf0bA=X';
Conf0bA=num2cell(Conf0bA);
Conf0bA{1,1}='Cantidad de brazo A';
[a0,s0]=size(Conf0pA);
[a1,s1]=size(Conf0pB);
[a2,s2]=size(Conf0bA);
Conf0=vertcat(horzcat('Configuracion de corte para pierna tipo A',repmat({' '},[1,max([s0,s1,s2])-1])),...
    horzcat(Conf0pA,repmat({' '},[a0,max([s0,s1,s2])-s0])),repmat({' '},[1,max([s0,s1,s2])]),...
    horzcat('Configuracion de corte para pierna tipo B',repmat({' '},[1,max([s0,s1,s2])-1])),...
    horzcat(Conf0pB,repmat({' '},[a1,max([s0,s1,s2])-s1])),repmat({' '},[1,max([s0,s1,s2])]),...
    horzcat('Configuracion de corte para brazo tipo A',repmat({' '},[1,max([s0,s1,s2])-1])),...
    horzcat(Conf0bA,repmat({' '},[a2,max([s0,s1,s2])-s2])));
xlswrite('resultadoscarne.xls',Conf0,'configuracion 0');
xlswrite('resultadoscarne.xls',[tipo',x0],'resultado de O.L.');
if ~isempty(demp)&&~isempty(dem0)
    p1=setdiff(1:size(demp,1),w1);
    p2=setdiff(1:size(dem0,1),w2);
    ut=[[demp(p1,2),zeros(size(p1,2),1)];[demp(w1,2),Mi*x0];...
        [dem0(w2,2),Meq0*x0];[dem0(p2,2),zeros(size(p2,2),1)]];
elseif isempty(demp)
    p2=setdiff(1:size(dem0,1),w2);
    ut=[[dem0(w2,2),Meq0*x0];[dem0(p2,2),zeros(size(p2,2),1)]];
else
    p1=setdiff(1:size(demp,1),w1);
    ut=[[demp(p1,2),zeros(size(p1,2),1)];[demp(w1,2),Mi*x0]];
end    
[la,idx]=ismember(cell2mat(demstr(2:end,4)),ut(:,1));
demstr=horzcat(demstr,vertcat('solucion lineal 0 (cumplimiento de oferta de 90%, sobredemanda de 8% total y oferta minima x corte)',num2cell(ut(idx(la),2))));
conf=[1:sum(tipo<0);full(tipo(tipo<0))]';
d2=size(ti3,2);
rel=find([rel1,rel2,rel3]);
intcon=repmat('I',[1,(d0+d1+d2)]);
intcon(rel)='C';
% Integer Constraints
xtype=intcon;
intcon=find(intcon=='I');
if ~isempty(Mi)&&~isempty(Meq0)
    [x2,fval2,exitflag2,OUTPUT]=intlinprog(-f,intcon,[-Mi;[-1;1]*sum([Mi;Meq0],1)],[-bi.*fi;[-.9;1.08]*sum(bi)],[Meq;Meq0;EE],zeros(size([Meq;Meq0;EE],1),1),lb,ub);
    A=[Mi;sum([Mi;Meq0],1);Meq;Meq0;EE];
    ru=[inf(size(bi,1),1);1.08*sum(bi);zeros(size([Meq;Meq0;EE],1),1)];
    rl=[bi.*fi;.9*sum(bi);zeros(size([Meq;Meq0;EE],1),1)];
elseif isempty(Meq0)
    [x2,fval2,exitflag2,OUTPUT]=intlinprog(-f,intcon,[-Mi;[-1;1]*sum(Mi,1)],[-bi.*fi;[-.9;1.08]*sum(bi)],[Meq;EE],zeros(size([Meq;EE],1),1),lb,ub);
    A=[Mi;sum(Mi,1);Meq;EE];
    ru=[inf(size(bi,1),1);1.08*sum(bi);zeros(size([Meq;EE],1),1)];
    rl=[bi.*fi;.9*sum(bi);zeros(size([Meq;EE],1),1)];
else
    [x2,fval2,exitflag2,OUTPUT]=intlinprog(-f,intcon,[],[],[Meq;Meq0;EE],zeros(size([Meq;Meq0;EE],1),1),lb,ub);
    A=[Meq;Meq0;EE];
    ru=zeros(size([Meq;Meq0;EE],1),1);
    rl=ru;
end
% Build OPTI Object
Opt=opti('f',-f,'lin',A,rl,ru,'bounds',lb,ub,'xtype',xtype);
% Solve the MILP problem
[x1,fval1,exitflag,info]=solve(Opt);
if ~isempty(x2)&&(norm(round(x1)-x1)>norm(round(x2)-x2))
    x2=x1;
    fval1=fval2;
end
intsol=-fval1;
if ~isempty(x1)
    intsol=[intsol;-fval1];
    if ~isempty(demp)&&~isempty(dem0)
        p1=setdiff(1:size(demp,1),w1);
        p2=setdiff(1:size(dem0,1),w2);
        ut=[[demp(p1,2),zeros(size(p1,2),1)];[demp(w1,2),Mi*x1];...
            [dem0(w2,2),Meq0*x1];[dem0(p2,2),zeros(size(p2,2),1)]];
    elseif isempty(demp)
        p2=setdiff(1:size(dem0,1),w2);
        ut=[[dem0(w2,2),Meq0*x1];[dem0(p2,2),zeros(size(p2,2),1)]];
    else
        p1=setdiff(1:size(demp,1),w1);
        ut=[[demp(p1,2),zeros(size(p1,2),1)];[demp(w1,2),Mi*x1]];
    end   
    [la,idx]=ismember(cell2mat(demstr(2:end,4)),ut(:,1));
    demstr=horzcat(demstr,vertcat('solucion lineal entera 1 (cumplimiento de oferta de 90%, sobredemanda de 8% total y oferta minima x corte)',num2cell(ut(idx(la),2))));
    %configuracion 1    
    r1=size(C1,1)-1;
    Conf0pA=C1([1;intersect(1:r1,w)+1],:);
    Conf0pA=[Conf0pA,[0;x1(find(ti1<0))]];
    Conf0pA=Conf0pA([1;find(Conf0pA(2:end,end)>0)+1],:);
    Conf0pA=full(Conf0pA(:,any(Conf0pA(2:end,:)>0,1)));
    a=unique(Conf0pA(1,:),'sorted');
    X=[];
    for j=a
        X=[X,[j;sum(Conf0pA(2:end,Conf0pA(1,:)==j),2)]];
    end
    Conf0pA=X';
    Conf0pA=num2cell(Conf0pA);
    Conf0pA{1,1}='Cantidad de pierna A';
    r2=size(C2,1)-1;
    Conf0pB=C2([1;intersect((1:r2)+r1,w)-r1+1],:);
    d0=size(ti1,2);
    Conf0pB=[Conf0pB,[0;x1(find(ti2<0)+d0)]];
    Conf0pB=Conf0pB([1;find(Conf0pB(2:end,end)>0)+1],:);
    Conf0pB=full(Conf0pB(:,any(Conf0pB(2:end,:)>0,1)));
    a=unique(Conf0pB(1,:),'sorted');
    X=[];
    for j=a
        X=[X,[j;sum(Conf0pB(2:end,Conf0pB(1,:)==j),2)]];
    end
    Conf0pB=X';
    Conf0pB=num2cell(Conf0pB);
    Conf0pB{1,1}='Cantidad de pierna B';
    r3=size(C3,1)-1;
    Conf0bA=C3([1;intersect((1:r3)+r1+r2,w)-r1-r2+1],:);
    d1=size(ti2,2);
    Conf0bA=[Conf0bA,[0;x1(find(ti3<0)+d0+d1)]];
    Conf0bA=Conf0bA([1;find(Conf0bA(2:end,end)>0)+1],:);
    Conf0bA=full(Conf0bA(:,any(Conf0bA(2:end,:)>0,1)));
    a=unique(Conf0bA(1,:),'sorted');
    X=[];
    for j=a
        X=[X,[j;sum(Conf0bA(2:end,Conf0bA(1,:)==j),2)]];
    end
    Conf0bA=X';
    Conf0bA=num2cell(Conf0bA);
    Conf0bA{1,1}='Cantidad de brazo A';
    [a0,s0]=size(Conf0pA);
    [a1,s1]=size(Conf0pB);
    [a2,s2]=size(Conf0bA);
    Conf0=vertcat(horzcat('Configuracion de corte para pierna tipo A',repmat({' '},[1,max([s0,s1,s2])-1])),...
        horzcat(Conf0pA,repmat({' '},[a0,max([s0,s1,s2])-s0])),repmat({' '},[1,max([s0,s1,s2])]),...
        horzcat('Configuracion de corte para pierna tipo B',repmat({' '},[1,max([s0,s1,s2])-1])),...
        horzcat(Conf0pB,repmat({' '},[a1,max([s0,s1,s2])-s1])),repmat({' '},[1,max([s0,s1,s2])]),...
        horzcat('Configuracion de corte para brazo tipo A',repmat({' '},[1,max([s0,s1,s2])-1])),...
        horzcat(Conf0bA,repmat({' '},[a2,max([s0,s1,s2])-s2])));    
    xlswrite('resultadoscarne.xls',Conf0,'configuracion 1');
end    
% Meqr=-Meq;
% Meqr(:,tipo<0)=.95*Meqr(:,tipo<0);
xlswrite('resultadoscarne.xls',demstr,'demanda');
xlswrite('resultadoscarne.xls',[tipo',x1],'resultado de O.L.E.');