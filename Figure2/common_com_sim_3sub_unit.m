function [Tmodel Ymodel]=common_com_sim_nouv(community,n_spc,SPC_EM,x_ini0)
global lxini spc 
global SxZ_tot kmax KM ke alpha beta subs_cn biom_inx met_udf met_udf_co n_carbon
global Y_BM Y_SUB maxEnzyme kl
n_s=n_spc;
kl={};
kmax={};
KM={};
ke={};
alpha={};
beta={};
subs_cn={};
biom_coef={};
n_carbon={};
e_rel={};
met_udf={};
biom_inx={};
Y_BM={};
Y_SUB={};
maxEnzyme={};
spc={};
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
     n_carbon{s}=spc{s}.n_carbon;  %%
     [cr,cl]=size(n_carbon{s}); %%
     if (cr>1)&&(cl>1)
        disp('this species use more than one cl species');%%
     end
     biom_inx{s}=strmatch('BIOM',met_udf{s});
 
    if ismember('expdata',fields(spc{s}))
        Time{s}=Expdata{s}(:,1);    %%
        Data{s}=Expdata{s}(:,2:end);
        [rexp cexp]=size(Data{s});
        Biomass{s}=Data{s}(:,biom_inx{s});    
        Data{s}(:,biom_inx{s})=Data{s}(:,biom_inx{s}).*biom_coef{s};
        Substrate{s}=Data{s}(:,subs_cn{s}); %% 
        temp_ini=Data{s}(1,:);
        x_ini{s}=temp_ini';
    else
    end
    Tspan{s}=Time{s};
    [sr,sl]=size(SxZ{s});
    n_ezm=sl;


for i=1:size(SPC_EM,2)
    if i==1
        stat_inx=1;
        end_inx=SPC_EM(s,i);
    else
        stat_inx=end_inx+1;
        end_inx=end_inx+SPC_EM(s,i-1);
    end
    kl{s}(i)=SPC_EM(s,i);
end

    Y_BM{s}=SxZ{s}(1,:)';
    Y_SUB{s}=SxZ{s}(subs_cn{s},:)';     %%%%
    maxmue{s}=kmax{s}.*Y_BM{s};
    maxEnzyme{s}=(ke{s}+alpha{s})./(beta{s}+maxmue{s}); 
    e0{s}=e_rel{s}.*maxEnzyme{s};      %
    para{s}=kmax{s};
    options_ode=[];
 end
  KM=K_MM;
%%%
for s=1:n_s
    met_udf{s}{1}=[met_udf{s}{1} '_' name{s}];
end

met_udf_co={};
bm_co={};
for s=1:n_s
        met_udf_co=union(met_udf_co,met_udf{s},'stable');
        bm_co=union(bm_co,met_udf{s}(1),'stable');
%         bm_co(i)=met_udf{s}(1);
end
clear_index=[];
for i=1:length(met_udf_co)   
    if ~isempty(strmatch('BIOMASS',met_udf_co{i}))
        clear_index=[clear_index;i];
    end
end
met_udf_co(clear_index)=[]; %
met_udf_co=union(bm_co,met_udf_co,'stable'); %% 

[cm_met,SxZ_tot,x_ini_mdf]=calc_SxZ(met_udf_co,met_udf,SxZ,n_s,x_ini)
x_ini=x_ini_mdf;

x_ini_co=zeros(length(x_ini{1}),1);
for s=1:n_s
    x_ini_co=x_ini_co+x_ini{s}; 
end
x_ini_co=x_ini0;
lxini=length(x_ini_co);
e_ini=[];
for s=1:n_s
    e_ini=[e_ini;e0{s}];
end

y0=[x_ini_co; e_ini];

%%% deal with simulation Time %%%
time=[];
for s=1:n_s
    time(s)=Time{s}(end);
end

tspan=[0 30];
options_ode=[];
[T Y]=ode15s(@qssm,tspan,y0,options_ode,kmax,n_s);
ut={};  %%%
vt={};  %%%

Tmodel=T;
Ymodel=Y(:,1:length(x_ini_co));

for ij=1:2
    Ymodel(:,ij)=Ymodel(:,ij)./biom_coef{ij};
end
%%

[ys,yl]=size(Ymodel);
if yl==cexp
    bmflag=1;  %%%%% plot biom
else
    bmflag=0;
end

for p=1:ys
    for q=1:yl
        if Ymodel(p,q)<0
            Ymodel(p,q)=0;
        end
    end
end


 
    
end
 
function dy=qssm(t,y,para,ns)
global lxini spc 
global SxZ_tot kmax KM ke alpha beta subs_cn biom_inx  met_udf met_udf_co n_carbon
global Y_BM Y_SUB maxEnzyme kl


for i=1:length(y)
    if y(i)<0,y(i)=eps;end
end

%%%%%%%% to identify substrates for each species %%%%%%%%%%%%%%%%
Sub=[];
for i=1:ns
    for j=1:length(subs_cn{i})
       sub_inx(j)=strmatch(met_udf{i}(subs_cn{i}(j)),met_udf_co);%%
    end
       Sub_tt=y(sub_inx);
       for j=1:length(Sub_tt)   %%% number of substrate
          S{j}=Sub_tt(j).*ones(kl{i}(j),1);  %%%
          Sub=[Sub;S{j}];
       end
       substrate{i}=Sub;
       sub_inx=[];
       Sub=[];
end

lx=length(met_udf_co);
for i=1:ns
    e{i}=y(lx+1:lx+length(kmax{i}));
    lx=lx+length(kmax{i});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    BIOM=y(s);  %%
    SxZ=SxZ_tot{s};
    klz=length(kmax{s});
S=[];

    S=sub;

        rkin=rkin.* S ./(KS + S); %%%%  
        re_s=re_s.*S ./(KS + S);
    r=rkin.*e_rel;
    rc=n_carbon{s};
     for i=1:length(kl{s})
        rc(:,i)=-Y_SUB{s}(:,i).*n_carbon{s};
    end
    rc=sum(rc,2);  %% if more than one column, add each row
    r_c=r.*rc;
    roi=e_rel.*rkin;
            roi=rc.*roi;

    roi=max(roi,zeros(klz,1));
    pu=max(roi,zeros(klz,1)); pv=max(roi,zeros(klz,1));
    sumpu=sum(pu); maxpv=max(pv);
    if sumpu>0, u=pu/sumpu; else u=zeros(klz,1); end
    if maxpv>0, v=pv/maxpv; else v=zeros(klz,1); end
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
 
end
%     
