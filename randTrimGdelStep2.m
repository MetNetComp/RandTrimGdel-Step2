function [gvalueList] = randTrimGdelStep2(model,targetMet,givenGvalue,max_it)
% randTrimGdelStep2 randomly and iteatively trimms gene deletions to 
% obtain the minimal gene deletion strategies when gDel_minRN (Step 1 of RandTrimGdel)
% strategy is given.
%
%INPUTS
%model:     COBRA model structure containing the following required fields to perform gDel_minRN.
%   rxns                    Rxns in the model
%   mets                    Metabolites in the model
%   genes               Genes in the model
%   grRules            Gene-protein-reaction relations in the model
%   S                       Stoichiometric matrix (sparse)
%   b                       RHS of Sv = b (usually zeros)
%   c                       Objective coefficients
%   lb                      Lower bounds for fluxes
%   ub                      Upper bounds for fluxes
%   rev                     Reversibility of fluxes
%
% targetMet:   target metabolites
% givenGvalue: gene deletion strategy obtained by gDel_minRN                          
%              An example is available in "gDel-minRN-strategy.mat"
%
% max_it: the number of iterations
%
%OUTPUTS
% gvalueList: a list of gene deletion strategies
%
% See example.m for an example.
%  
%
%  Mar. 1, 2023   Takeyuki TAMURA
%


tic

options=cplexoptimset('cplex');
options.mip.tolerances.integrality=10^(-12);
%sss=sprintf('example1.mat');

[model,targetRID,extype] = modelSetting(model,targetMet)

m=size(model.mets,1);
n=size(model.rxns,1);
g=size(model.genes,1);
gid=find(model.c);
pid=targetRID;

replaceList=findReplacableGenes(model);

gvalue=givenGvalue;
for i=1:g
    if gvalue{i,2}==0 && replaceList(i,1)~=0 && gvalue{replaceList(i,1),2}==1
        gvalue{i,2}=1;
        gvalue{replaceList(i,1),2}=0;
    end
end


model2=model;

[grRules0] = calculateGR(model,gvalue);
lb2=model.lb;
ub2=model.ub;

for i=1:n
    if grRules0{i,4}==0
        lb2(i)=0;
        ub2(i)=0;
    end
end
[opt0.x, opt0.f, opt0.stat, opt0.output] = ...
    cplexlp(-model.c, [],[], model.S, zeros(m,1),lb2, ub2);
[opt0.x(gid) opt0.x(pid)]

GR0=-opt0.f;
lb2(gid)=GR0;
ub2(gid)=GR0;
model2.c(gid)=0;
model2.c(pid)=1;
[opt1.x, opt1.f, opt1.stat, opt1.output] = ...
    cplexlp(-model2.c, [],[], model.S, zeros(m,1),lb2, ub2);
[opt1.x(gid) opt1.x(pid)]
GRLB=opt1.x(gid);
PRLB=opt1.x(pid);
if GRLB<0.001 || PRLB<0.001
    display('error in gene replacement.\n');return;
end
%GRLB=0.001+ it*(GRLB-0.001)/max_it;
%PRLB=0.001+ it*(PRLB-0.001)/max_it;

[term,ng,nt,nr,nko,reactionKO,reactionKO2term] = readGeneRules(model);
[f,intcon,A,b,Aeq,beq,lb,ub,xname] = geneReactionMILP(model,term,ng,nt,nr,nko,reactionKO);

lp.Aeq=Aeq;
lp.beq=[zeros(size(lp.Aeq,1),1)];
j=1;
for i=1:size(model.grRules,1)
    if isempty(model.grRules{i,:})==0
        ind(1,j)=i;
        j=j+1;
    end
end
z1=-diag(model.ub);
z2=diag(model.lb);
z3=eye(n);
lp.A=A;
lp.b=b;
lp.lb=lb;
lp.ub=ub;

for i=1:ng
    if gvalue{i,2}==1
        lp.lb(i,1)=1;
        lp.ub(i,1)=1;
    end
end

[grRules0] = calculateGR(model,gvalue);

j=1;
for i=1:n
    if isempty(model.grRules{i,1})==0
        lp.lb(ng+nt+j,1)=grRules0{i,4};
        lp.ub(ng+nt+j,1)=grRules0{i,4};
        j=j+1;
    end
end

lp.f=[-ones(ng,1); zeros(nt,1); zeros(nko,1)];
for i=1:n
    s2=repelem('B',ng+nt+nko);
    lp.ctype=sprintf('%s%s',s2);
end

[opt.x, opt.f, opt.stat, opt.output] = ...
    cplexmilp(lp.f, lp.A, lp.b, lp.Aeq, lp.beq,[],[],[],lp.lb, lp.ub,lp.ctype,[],options);


if opt.stat>0
    for i=1:ng
        vg(i,1)=opt.x(i);
        gvalue{i,2}=opt.x(i);
    end
    for i=1:nt
        vt(i,1)=opt.x(ng+i);
    end
    for i=1:nko
        vko(i,1)=opt.x(ng+nt+i);
    end
end




model2=model;
grprlist(1,1)=opt1.x(gid);
grprlist(1,2)=opt1.x(pid);
grprlist(1,3)=opt1.x(gid);
grprlist(1,4)=opt1.x(pid);

for i=1:g
    if gvalue{i,2}==0 && replaceList(i,1)~=0 && gvalue{replaceList(i,1),2}==1
        gvalue{i,2}=1;
        gvalue{replaceList(i,1),2}=0;
    end
end

gvalue0=gvalue;

for it=1:max_it
    it
    gvalue=gvalue0;
    
    k=2;
    loop=0;
    c=g;
    while loop<c
        loop=loop+1;
        if loop==1
            zeroList=find(cell2mat(gvalue(:,2))==0);
            if isempty(zeroList)==1
                gvalueList{1,1}=gvalue;
                I2=1;I3=1;I4=1;
                return;
            end
            nonReplace=find(replaceList==0);
            candidates=intersect(zeroList,nonReplace);
            c=size(candidates,1); 
            randCandidates=candidates(randperm(c));
        end
      
        i=randCandidates(loop,1);
        if gvalue{i,2}==1
            save('a.mat');display('error');return;
        end
        gvalue{i,2}=1;
        
        [grRules2] = calculateGR(model,gvalue);
        lb2=model.lb;
        ub2=model.ub;
        
        for j=1:n
            if grRules2{j,4}==0
                lb2(j)=0;
                ub2(j)=0;
            end
        end
        [opt2.x, opt2.f, opt2.stat, opt2.output] = ...
            cplexlp(-model.c, [],[], model.S, zeros(m,1),lb2, ub2);
        grprlist(k,1)=opt2.x(gid);
        grprlist(k,2)=opt2.x(pid);
        GR2=-opt2.f;
        lb2(gid)=GR2;
        ub2(gid)=GR2;
        model2.c(gid)=0;
        model2.c(pid)=1;
        [opt3.x, opt3.f, opt3.stat, opt3.output] = ...
            cplexlp(model2.c, [],[], model.S, zeros(m,1),lb2, ub2);
        grprlist(k,3)=opt3.x(gid);
        grprlist(k,4)=opt3.x(pid);

        if  ((opt3.x(gid) < 0.001) || (opt3.x(pid) < 0.001))
            gvalue{i,2}=0;
            grprlist(k,:)=grprlist(k-1,:);
        else
            loop=0;
        end
        k=k+1;
    end

    
    
    for i=1:g
        if gvalue{i,2}==0 && replaceList(i,1)~=0 && gvalue{replaceList(i,1),2}==1
            gvalue{i,2}=1;
            gvalue{replaceList(i,1),2}=0;
        end
    end
    

    %gvalueList=horzcat(givenGvalue(:,2),gvalue0(:,2),gvalue(:,2));
    %size1=size(find(cell2mat(givenGvalue(:,2))==0),1);
    %size2=size(find(cell2mat(gvalue0(:,2))==0),1);
    size3=size(find(cell2mat(gvalue(:,2))==0),1);
    
    if it==1
        gvalueList{it,1}=gvalue;
        gvalueList{it,2}=size(find(cell2mat(gvalue(:,2))==0),1);
        gvalueList{it,3}=grprlist(k-1,3);
        gvalueList{it,4}=grprlist(k-1,4);
    else
        flag=0;
        s=size(gvalueList,1);
        for i=1:s
            if isequal(cell2mat(gvalueList{i,1}(:,2)),cell2mat(gvalue(:,2)))==1
                flag=1;
            end
        end
        if flag==0
            gvalueList{s+1,1}=gvalue;
            gvalueList{s+1,2}=size(find(cell2mat(gvalue(:,2))==0),1);
            gvalueList{s+1,3}=grprlist(k-1,3);
            gvalueList{s+1,4}=grprlist(k-1,4);
            %save(sss,'gvalueList');
        end
    end
    
end


gs=size(gvalueList,1);

for ggg=1:gs
    
    gvalue=gvalueList{ggg,1};
    
    zeroList=find(cell2mat(gvalue(:,2))==0);
    rep=find(replaceList==0);
    candidates=intersect(zeroList,rep);
    c=size(candidates,1);
    
    if c~=0
        
        model2=model;
        
        [grRules0] = calculateGR(model,gvalue);
        lb2=model.lb;
        ub2=model.ub;
        
        for i=1:n
            if grRules0{i,4}==0
                lb2(i)=0;
                ub2(i)=0;
            end
        end
        [opt0.x, opt0.f, opt0.stat, opt0.output] = ...
            cplexlp(-model.c, [],[], model.S, zeros(m,1),lb2, ub2);
        [opt0.x(gid) opt0.x(pid)]
        
        GR0=-opt0.f;
        lb2(gid)=GR0;
        ub2(gid)=GR0;
        model2.c(gid)=0;
        model2.c(pid)=1;
        [opt1.x, opt1.f, opt1.stat, opt1.output] = ...
            cplexlp(-model2.c, [],[], model.S, zeros(m,1),lb2, ub2);
        [opt1.x(gid) opt1.x(pid)]
        GRLB=opt1.x(gid);
        PRLB=opt1.x(pid);
        if GRLB<0.001 || PRLB<0.001
            display('error in gene replacement.\n');return;
        end
        %GRLB=0.001+ it*(GRLB-0.001)/max_it;
        %PRLB=0.001+ it*(PRLB-0.001)/max_it;
        
        
        
        gvalue0=gvalue;
        
        
        gvalue=gvalue0;
        
        k=2;
        loop=0;
        c=g;
        while loop<c
            loop=loop+1;
            if loop==1
                zeroList=find(cell2mat(gvalue(:,2))==0);
                if isempty(zeroList)==1
                    gvalueList{1,1}=gvalue;
                    I2=1;I3=1;I4=1;
                    return;
                end
                rep=find(replaceList==0);
                candidates=intersect(zeroList,rep);
                c=size(candidates,1);
                
                
                randCandidates=candidates(randperm(c));
            end
            
            i=randCandidates(loop,1);
            if gvalue{i,2}==1
                save('a.mat');display('error');return;
            end
            gvalue{i,2}=1;
            
            [grRules2] = calculateGR(model,gvalue);
            lb2=model.lb;
            ub2=model.ub;
            
            for j=1:n
                if grRules2{j,4}==0
                    lb2(j)=0;
                    ub2(j)=0;
                end
            end
            [opt2.x, opt2.f, opt2.stat, opt2.output] = ...
                cplexlp(-model.c, [],[], model.S, zeros(m,1),lb2, ub2);
            grprlist(k,1)=opt2.x(gid);
            grprlist(k,2)=opt2.x(pid);
            GR2=-opt2.f;
            lb2(gid)=GR2;
            ub2(gid)=GR2;
            model2.c(gid)=0;
            model2.c(pid)=1;
            [opt3.x, opt3.f, opt3.stat, opt3.output] = ...
                cplexlp(model2.c, [],[], model.S, zeros(m,1),lb2, ub2);
            grprlist(k,3)=opt3.x(gid);
            grprlist(k,4)=opt3.x(pid);
            
            if  ((opt3.x(gid) < 0.001) || (opt3.x(pid) < 0.001))
                gvalue{i,2}=0;
                grprlist(k,:)=grprlist(k-1,:);
            else
                loop=0;
            end
            k=k+1;
        end
        
        
        for i=1:g
            if gvalue{i,2}==0 && replaceList(i,1)~=0 && gvalue{replaceList(i,1),2}==1
                gvalue{i,2}=1;
                gvalue{replaceList(i,1),2}=0;
            end
        end
        
        flag=0;
        for i=1:ggg-1
            if isequal(cell2mat(gvalueList{i,1}(:,2)),cell2mat(gvalue(:,2)))==1
                flag=1;
            end
        end
        if flag==0
            gvalueList{ggg,1}=gvalue;
            gvalueList{ggg,2}=size(find(cell2mat(gvalue(:,2))==0),1);
            [GR ,PR] = GRPRchecker(model,targetMet,gvalue);
            gvalueList{ggg,3}=GR;
            gvalueList{ggg,4}=PR;
        else
            gvalueList{ggg,1}=[];
            gvalueList{ggg,2}=99999;
            gvalueList{ggg,3}=0;
            gvalueList{ggg,4}=0;
        end
    end
end


tmp=9999*cell2mat(gvalueList(:,2))-cell2mat(gvalueList(:,4))-0.0001*cell2mat(gvalueList(:,3));
[B,I2]=sort(tmp);
gvalueList=gvalueList(I2,:);

system('rm -f clone*.log');
time=toc;
%save(sss,'gvalueList');
return;
end

