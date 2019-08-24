% function [Pure_Growth,Comm_Growth]=common_com_sim(community,n_spc,vld,uv)
function [Comm_Growth,rgt]=common_com_sim(community,n_spc)

global lxini spc 
global SxZ_tot kmax KM ke alpha beta subs_cn biom_inx met_udf met_udf_co n_carbon
global Y_BM Y_SUB maxEnzyme kl n_s
n_s=n_spc;
kl={};
 for s=1:n_s
     spc{s}=community{s};
     if ismember('expdata',fields(spc{s}))
         Expdata{s}=spc{s}.expdata;
     else 
         disp('no available experimental data for "''i"th species!');
         x_ini{s}=spc{s}.x_ini;
         Time{s}=spc{s}.Time;
     end
     name{s}=spc{s}.id;
     SxZ{s}=spc{s}.S;
     kmax{s}=spc{s}.kmax;
     K_MM{s}=spc{s}.KM;
     ke{s}=spc{s}.ke;
     alpha{s}=spc{s}.alpha;
     beta{s}=spc{s}.beta;
     ke{s}=spc{s}.ke;
     met_udf{s}=spc{s}.met_udf;
     e_rel{s}=spc{s}.e_rel;
     subs_cn{s}=spc{s}.subs_cn;
     biom_coef{s}=spc{s}.biom_coef;
     n_carbon{s}=spc{s}.n_carbon;  
     [cr,cl]=size(n_carbon{s}); 
     if (cr>1)&&(cl>1)
        disp('this species use more than one cl species');
     end
     biom_inx{s}=strmatch('BIOM',met_udf{s});
 
    if ismember('expdata',fields(spc{s}))
        Time{s}=Expdata{s}(:,1);    %% this should be fixed format of species_info
        Data{s}=Expdata{s}(:,2:end);
        Data{s}(find(isnan(Time{s})),:)=[];
        Time{s}(find(isnan(Time{s})))=[];
        
        [rexp cexp]=size(Data{s});
        Biomass{s}=Data{s}(:,biom_inx{s}).*biom_coef{s};    
        Data{s}(:,biom_inx{s})=Data{s}(:,biom_inx{s}).*biom_coef{s};
        Substrate{s}=Data{s}(:,subs_cn{s}); 
        temp_ini=Data{s}(1,:);
        x_ini{s}=temp_ini';
        Pure_Growth{s}=num2cell([Time{s},Biomass{s}]);
        lx{s}=length(x_ini{s});
    else
    end
    Tspan{s}=Time{s};
    [sr,sl]=size(SxZ{s});
    n_ezm=sl;
    kl{s}=length(kmax{s});
    KM{s}=[];
    for i=1:length(subs_cn{s})
        KM{s}(:,i)=K_MM{s}(:,i).*ones(kl{s},1);  
    end
    alpha{s}=alpha{s}.*ones(kl{s},1);   
    beta{s}=beta{s}.*ones(kl{s},1);
    ke{s}=ke{s}.*ones(kl{s},1); 
    e_rel{s}=e_rel{s}.*ones(kl{s},1); 

    Y_BM{s}=SxZ{s}(1,:)';
    Y_SUB{s}=SxZ{s}(subs_cn{s},:)';     
    maxmue{s}=kmax{s}.*Y_BM{s};
    maxEnzyme{s}=(ke{s}+alpha{s})./(beta{s}+maxmue{s}); 
    e0{s}=e_rel{s}.*maxEnzyme{s};      % absolute enzyme level
    para{s}=kmax{s};
    options_ode=[];
 end

met_udf_co={};
bm_co={};
for s=1:n_s
        met_udf_co=union(met_udf_co,met_udf{s},'stable');
end
met_udf_co(1)=[];
for s=1:n_s
    met_udf{s}{1}=[met_udf{s}{1} '_' name{s}];
end
for s=1:n_s
        bm_co1=union(bm_co,met_udf{s}(1),'stable');
        bm_co=[bm_co;met_udf{s}(1)];
end
clear_index=[];
SxZ_tot={};
temp_xini={};
if length(bm_co1)==1   
    clear_index=[];
    for s=1:n_s
        SxZ_tot{s}(3:lx{s}+1,:)=SxZ{s}(2:end,:);
        SxZ_tot{s}(s,:)=SxZ{s}(1,:);
        x_ini{s}(3:lx{s}+1)=x_ini{s}(2:end);
        x_ini{s}(s)=x_ini{s}(1);
    end
    SxZ_tot{1}(2,:)=0;
    SxZ_tot{2}(1,:)=0;
    x_ini{1}(2)=0;
    x_ini{2}(1)=0;
    met_udf_co=[bm_co;met_udf_co]; 

else    
for i=1:length(met_udf_co)  
    if ~isempty(strmatch('BIOMASS',met_udf_co{i}))
        clear_index=[clear_index;i];
    end
end
met_udf_co(clear_index)=[]; 
met_udf_co=union(bm_co,met_udf_co,'stable');

[cm_met,SxZ_tot,x_ini_mdf]=calc_SxZ(met_udf_co,met_udf,SxZ,n_s,x_ini);
x_ini=x_ini_mdf;

end

x_ini_co=zeros(length(x_ini{1}),1);
inx_hex=find(strcmpi('Hexose',met_udf_co))
inx_ac=find(strcmpi('Acetate',met_udf_co))

for s=1:n_s
    x_ini_co=x_ini_co+x_ini{s}; 
end
if x_ini_co(inx_hex)>80
    x_ini_co(inx_hex)=80;
end
if x_ini_co(inx_ac)>60
    x_ini_co(inx_ac)=60;
end
lxini=length(x_ini_co);
e_ini=[];
for s=1:n_s
    e_ini=[e_ini;e0{s}];
end
y0=[x_ini_co; e_ini];
time=[];
for s=1:n_s
    time(s)=Time{s}(end);
end
% tmax=max(time);
% tspan=0:0.05:tmax;
tspan=[0 max(time)];
TSpan=union(Time{1},Time{2});
TSpan=0:0.02:24;

time_same=intersect(Time{1},Time{2});
[y,inx1]=ismember(time_same,Time{1});
[y,inx2]=ismember(time_same,Time{2});
Pure_Growth=num2cell([Biomass{1}(inx1),Biomass{2}(inx2)]);
id_pair=name;
Pure_Growth=[id_pair;Pure_Growth];
options_ode=[];
[T Y]=ode15s(@qssm,TSpan,y0,options_ode,kmax,n_s);

Comm_Growth=[Y(:,1),Y(:,2)];

ut={};  %%% cell of u vector of all species
vt={};  %%% cell of u vector of all species
rgt=[];
for s=1:n_s   
    u=[];
    v=[];
    rg=[];
    for i=1:length(TSpan)
       [ui,vi,rgi]=plot_uv(s,Y(i,:));
       u=[u;ui'];
       v=[v;vi'];
       rg=[rg;rgi'];
     end
    ut{s}=u;
    vt{s}=v;
    rgt(:,s)=rg;
end

function dy=qssm(t,y,para,ns)
global lxini spc 
global SxZ_tot kmax KM ke alpha beta subs_cn biom_inx  met_udf met_udf_co n_carbon
global Y_BM Y_SUB maxEnzyme kl
% i_strain=mysys.i_strain;
% kmax1=para(1:kl1);
% kmax2=para(kl1+1:kl1+kl2);
% kmax3=para(kl1+kl2+1:end);

for i=1:length(y)
    if y(i)<0,y(i)=eps;end
end

%%%%%%%% to identify substrates for each species %%%%%%%%%%%%%%%%
for i=1:ns
    for j=1:length(subs_cn{i})
       sub_inx(j)=strmatch(met_udf{i}(subs_cn{i}(j)),met_udf_co);%% strmatch cannot identify two strings at one time
    end
       substrate{i}=y(sub_inx);
end

lx=length(met_udf_co);
for i=1:ns
    e{i}=y(lx+1:lx+length(kmax{i}));
    lx=lx+length(kmax{i});
end

%%% if one species use more than one substrate, make sure KM should also more than one column
for i=1:ns
    if length(subs_cn{i}>1)
        [rkm,ckm]=size(KM{i});
        if ckm==1
            for j=1:length(subs_cn{i})
                KM{i}(:,j)=KM{i}(:,1);
            end
        end
    end
end
u_co={};
v_co={};
rg_co={};
dxdt={};
re_co={};
for s=1:ns
    e_rel=e{s}./maxEnzyme{s};
    rkin=kmax{s};
    re_s=ke{s};
    sub=substrate{s};
    KS=KM{s};
    BIOM=y(s);  
    SxZ=SxZ_tot{s};
    klz=length(kmax{s});
S=[];
for i=1:length(sub)   %%% number of substrate
    S(:,i)=sub(i).*ones(klz,1); 
end
    for ss=1:length(sub)  
        rkin=rkin.* S(:,ss) ./(KS(:,ss) + S(:,ss)); %%%%  
        re_s=re_s.*S(:,ss) ./(KS(:,ss) + S(:,ss));
    end 
    r=rkin.*e_rel;
    rc=n_carbon{s};
    rc=-rc.*Y_SUB{s};
    rc=sum(rc,2);  %% if more than one column, add each row
    r_c=r.*rc;
    roi=e_rel.*rkin;
            roi=rc.*roi;

    roi=max(roi,zeros(kl{s},1));
    pu=max(roi,zeros(kl{s},1)); pv=max(roi,zeros(kl{s},1));
    sumpu=sum(pu); maxpv=max(pv);
    if sumpu>0, u=pu/sumpu; else u=zeros(kl{s},1); end
    if maxpv>0, v=pv/maxpv; else v=zeros(kl{s},1); end
    u_co{s}=u;
    v_co{s}=v;
    rM=v.*e_rel.*rkin;
    rg_co{s}=sum(Y_BM{s}.*v.*r);
    re_co{s}=re_s;
    diagV=diag(v);
    dxdt{s}=SxZ*diagV*r*BIOM;
end
    dy=zeros(length(y),1);
    lx=length(met_udf_co);
    for s=1:ns        
        dy(1:lx)=dy(1:lx)+dxdt{s}(1:lx);
    end

    for s=1:ns        
        dy(lx+1:lx+length(kmax{s}))=alpha{s}+re_co{s}.*u_co{s}-(beta{s}+rg_co{s}).*e{s};
        lx=lx+length(kmax{s});
    end
 
%         lx=lx+s.length(kmax{s});
%     
%             dy=zeros(length(y),1);
%         dy(1)=dxdt1(1);  %% biomass of bth
%         dy(2)=dxdt2(1);  %% biomass of bif
%         dy(3)=dxdt3(1);  %% biomass of eha
%         dy(4:10)=dxdt1(2:8)+dxdt2(2:8)+dxdt3(2:8);
%         dy(11:11+kl1-1)=alpha{s}+re1.*u1-(beta_bth+rg1).*e1;
%         dy(11+kl1:11+kl1+kl2-1)=alpha_bif+re2.*u2-(beta_bif+rg2).*e2;
%         dy(11+kl1+kl2:end)=alpha_eha+re3.*u3-(beta_eha+rg3).*e3;
%     

