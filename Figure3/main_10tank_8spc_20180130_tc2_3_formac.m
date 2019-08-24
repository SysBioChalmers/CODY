function run_10tank
k=[1.2    23    12    0.01     10     6      10     6];
K=[45    15    10     120     12    120     35    120 ];
main_10tank(k,K,1);
function main_10tank(k_hyd_5spc,K_XZ_5spc,round)
addpath('/Users/gejun/Desktop/Dynamic_Gut_Microbiome/Research/BMGF project/Gut_infant/Dynami_Gut_Model/Species_Model')
addpath('/Users/gejun/Desktop/Dynamic_Gut_Microbiome/Research/BMGF project/Gut_infant/Dynami_Gut_Model')
SpeciesPath='/Users/gejun/Desktop/Dynamic_Gut_Microbiome/Research/BMGF project/Gut_infant/Dynami_Gut_Model/Species_Model/VLD_comm';

CurrentPath=pwd;
cd(SpeciesPath)
spc_id={'B_theta_20170903_6_1_VLD','B_fragilis_prop_20170907_5_VLD','Bifido_longum_20170901_VLD','Bifido_breve_20170907_3_VLD', 'Bifido_ado_20170907_3_VLD' 'Eubacterium_hallii_20170918' 'F_prausnitzii_1212_VLD' 'Roseburia_intes_VLD'};

n_species=length(spc_id);

load 8SPC_0405

e_ini=[];
e_ini_X=[];
 for s=1:n_species
     spc{s}=SPC{s};
     SxZ{s}=spc{s}.S;
     kmax{s}=spc{s}.kmax;
     e_rel{s}=spc{s}.e_rel;
     ke{s}=spc{s}.ke;
     alpha{s}=spc{s}.alpha;
     beta{s}=spc{s}.beta;
     Y_BM{s}=SxZ{s}(1,:)';
    maxmue{s}=kmax{s}.*Y_BM{s};
    maxEnzyme{s}=(ke{s}+alpha{s})./(beta{s}+maxmue{s}); 
    e0{s}=e_rel{s}.*maxEnzyme{s}; 
     met_udf{s}=spc{s}.met_udf;
     e_ini=[e_ini;e0{s}];
     e_ini_X=[e_ini_X;e0{s}];
 end
 
met_co={};
for s=1:n_species
        met_co=union(met_co,met_udf{s},'stable');
end
met_co(1)=[];   

Metabolites_list{1}={'HMO'  'hexose' 'succinate' 'acetate' 'lactate' 'formate' 'butyrate' 'ethanol' 'H2' 'propionate'};
Metabolites_list{2}={1000  180.16   118.09       59.04      90.08      46.03    87.098     46.07    2    73.07};  %% ?????????????????????mol_mass????????????????????????????1000?????????????????????????????????????????hmo??????????????1000??????????????????????????????????????????????????1
Pn_list={power(10,-6); 2*power(10,-6);power(10,-6);5*power(10,-6);power(10,-6);5*power(10,-6);5*power(10,-6);power(10,-6);0; 5*power(10,-6)};    %% permeability for each metabolite is different,[glucose succinate acetate lactate formate ethanol buytrate H2], unit: cm/sec,so transfer to dm/h,*3600/10;glycerol: 1.6*10^-7;ethanol:2.1*E^-3;water:1.2*E^-3

Km_list={0.3;         0.2;          0.1;        0.1;           0.1;         0.1;           0.2;          0.1;        0;   0.1      };
Pn_list={power(10,-6); 2*power(10,-6);power(10,-6); 10*power(10,-5);power(10,-6);5*power(10,-6);10*power(10,-5);power(10,-6);0; 10*power(10,-5)};

met_co_hmo=['HMO';met_co];
met_length=length(met_co);

for i=1:length(met_co_hmo)
    met_index(i)=find(strcmpi(met_co_hmo(i),Metabolites_list{1}));
end
Mass=cell2mat(Metabolites_list{2}(met_index));
Pn=cell2mat(Pn_list(met_index));
Km_lumen=cell2mat(Km_list(met_index));
ratio=[2.42  6.88   17.2   6.49    2.05    0.003  0.007  0.008]';   
total_ratio=sum(ratio);
ratio=ratio./total_ratio;
x50=zeros(n_species,1);          
x10=5.*ones(n_species,1).*ratio;      
x20=5.*ones(n_species,1).*ratio;
x30=5.*ones(n_species,1).*ratio;
x40=5.*ones(n_species,1).*ratio;

x50=0.00001.*ones(n_species,1).*ratio;
XB10=15.*ones(n_species,1).*ratio;     
XB20=10.*ones(n_species,1).*ratio;
XB30=4.*ones(n_species,1).*ratio;
XB40=2.*ones(n_species,1).*ratio;n_rct=10;
V_colon_t=4.5;  
V_colon=[0.15 0.15 0.15 0.15 0.03 0.03 0.03 0.03 0.15 1]';  

GLU_in=150;    
F0=[0.1  0.1  0  0.5]';  
f=[0  0.01  0.01  0.01]';   
kd=[0.075 0.03 0.015 0.008
    16 16 16 16]  
sim_time=600;             
ini_conc=1.*ones(met_length,1);
ini_conc_MCS=0.5.*ones(met_length,1);   
ini_conc_V5=zeros(met_length+1,1);      
ini_conc_BLD=0.05.*ones(met_length,1);       
M_biomass=113;  
k_hyd=1.2*power(10,3)/2400;
K_XZ=10;

Y_SZ=5.*ones(n_species,1);
HMO=80;  
hyd_spc=ones(n_species,1); 
k_hyd=k_hyd_5spc;   
K_XZ=K_XZ_5spc;
hmo_para=[HMO;hyd_spc; k_hyd; K_XZ; Y_SZ];
HMO_ini=0.5;
ini_conc=[HMO_ini*10;ini_conc];
ini_conc_MCS=[HMO_ini*5;ini_conc_MCS];
ini_conc_BLD=[HMO_ini;ini_conc_BLD];
ini_stat_Lumen1=[x10;ini_conc;e_ini];
ini_stat_Lumen2=[x20;ini_conc;e_ini];
ini_stat_Lumen3=[x30;ini_conc;e_ini];
ini_stat_Lumen4=[x40;ini_conc;e_ini];
ini_stat_Mucosa1=[XB10;ini_conc_MCS;e_ini_X];
ini_stat_Mucosa2=[XB20;ini_conc_MCS;e_ini_X];
ini_stat_Mucosa3=[XB30;ini_conc_MCS;e_ini_X];
ini_stat_Mucosa4=[XB40;ini_conc_MCS;e_ini_X];
lg_standard=length(ini_stat_Lumen1);   
ini_stat_V5=[x50;ini_conc_V5;e_ini;V_colon(9);x50.*V_colon(9);ini_conc_V5.*V_colon(9)];  
ini_stat_Blood=[ini_conc_BLD];
ini_stat=[ini_stat_Lumen1;ini_stat_Lumen2;ini_stat_Lumen3;ini_stat_Lumen4;ini_stat_Mucosa1;ini_stat_Mucosa2;ini_stat_Mucosa3;ini_stat_Mucosa4;ini_stat_V5;ini_stat_Blood;];

Km_water={1/2 1/6 1/40 1/70};

Pn=Pn.*3600/100; 
Pn=Pn.*10;
av=1/0.05; 
Pn_mucosa=Pn.*av;

density=0.7*1000;

[T,Y,Y_model,Y_mets,met_udf_co]=dynamic_simulation_4to12_11_formac_mdf(SPC,n_species,n_rct,GLU_in,F0,f,kd,V_colon,ini_stat,sim_time,lg_standard,Km_lumen,Pn_mucosa,Km_water,density,Mass,hmo_para);

ReferenceStates=StoreResults_10tank(SPC,n_species,n_rct,GLU_in,F0,f,kd,V_colon,ini_stat,sim_time,lg_standard,Km_lumen,Pn_mucosa,Km_water,density,T,Y,Y_model,Y_mets,met_udf_co,lg_standard);


