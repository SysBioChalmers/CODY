function [Tmodel,Ymodel,Time,Data,met_udf]=common_pure_sim_mulsub(species_info,spc_EM,vldsim)

global spc SPC_EM
global KM Y_BM Y_SUB x_ini 
global ke alpha beta 
global maxEnzyme e_rel0
global SxZ subs_cn n_ezm kl n_carbon biom_inx
global spc SPC_EM

 
spc=species_info;
SPC_EM=spc_EM;

 if ismember('expdata',fields(spc))
      Expdata=spc.expdata;
 else 
     disp('no available experimental data for this species!');
     x_ini=spc.x_ini;
     Time=spc.Time;
 end
 
 SxZ=spc.S;
 kmax=spc.kmax;
 K_MM=spc.KM;

 ke=spc.ke;
 alpha=spc.alpha;
 beta=spc.beta;
 met_udf=spc.met_udf;
 e_rel0=spc.e_rel;
 subs_cn=spc.subs_cn;
 biom_coef=spc.biom_coef;
 
 n_carbon=spc.n_carbon;  
 [cr,cl]=size(n_carbon); 
 if (cr>1)&&(cl>1)
 end
 
 biom_inx=strmatch('BIOM',met_udf);
 
 if ismember('expdata',fields(spc))
     Time=Expdata(:,1);    
     Data=Expdata(:,2:end);
     Data(any(isnan(Data)'),:) = [];  
     Time(isnan(Time))=[];
     Data(isnan(Data))=[];

     [mYexp nYexp]=size(Data);
     Biomass=Data(:,biom_inx);    
     Data(:,biom_inx)=Data(:,biom_inx).*biom_coef;
     Substrate=Data(:,subs_cn); 
     x_ini=Data(1,:);
 else
 end

Tspan=Time;
[sr,sl]=size(SxZ);
n_ezm=sl;
kl=length(kmax);
KM=[];
KM=K_MM;
for i=1:length(SPC_EM)
    if i==1
        stat_inx=1;
        end_inx=SPC_EM(i);
    else
        stat_inx=end_inx+1;
        end_inx=end_inx+SPC_EM(i-1);
    end

    kl(i)=SPC_EM(i);

end


Y_BM=SxZ(1,:)';
Y_SUB=SxZ(subs_cn,:)';     
maxmue=kmax.*Y_BM;
maxEnzyme=(ke+alpha)./(beta+maxmue); 
e0=e_rel0.*maxEnzyme;      %
para=kmax;
options_ode=[];
Y0=[x_ini e0'];
Tspan=[0 Time(end)];   
[T,Y]=ode15s(@diff,Tspan,Y0,options_ode,para);
Tmodel=T;
Ymodel=Y(:,1:length(x_ini));
Ymodel(:,biom_inx)=Ymodel(:,biom_inx)./biom_coef;
Data(:,biom_inx)=Data(:,biom_inx)./biom_coef;
for ij=2:4
    Ymodel(:,ij)=Ymodel(:,ij)./Ymodel(1,ij).*2;   %%
    Data(:,ij)=Data(:,ij)./Data(1,ij).*2;
end
[ys,yl]=size(Ymodel);

if yl==nYexp
    bmflag=1;  
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



function dy=diff(t,y,para) 
global KM Y_BM Y_SUB x_ini 
global ke alpha beta 
global maxEnzyme 
global SxZ subs_cn n_ezm kl n_carbon biom_inx spc 

kmax=para;
kmax_leg=length(kmax);
for i=1:length(y)
    if y(i)<0
        y(i)=eps;
    end
end

BIOM=y(biom_inx);
sub=y(subs_cn);
Sub_tt=[];
for i=1:length(sub)   %%% number of substrate
    S{i}=sub(i).*ones(kl(i),1);  
    Sub_tt=[Sub_tt;S{i}];
end

e=y(length(x_ini)+1:end); 
e_rel=e./maxEnzyme;
es=size(e);

rkin=kmax;
re=ke;

    rkin=rkin.* Sub_tt ./(KM+ Sub_tt); %%%%  
    re=re.*Sub_tt ./(KM + Sub_tt);    
    
    r=rkin.*e_rel;
    rc=n_carbon;
    for i=1:length(sub)
        rc(:,i)=-Y_SUB(:,i).*n_carbon;
    end
    rc=sum(rc,2);
    r_c=r.*rc;
    
    roi=e_rel.*rkin;
            roi=rc.*roi;

    roi=max(roi,zeros(kmax_leg,1));
    pu=max(roi,zeros(kmax_leg,1)); pv=max(roi,zeros(kmax_leg,1));
    sumpu=sum(pu); maxpv=max(pv);
    if sumpu>0, u=pu/sumpu; else u=zeros(kmax_leg,1); end
    if maxpv>0, v=pv/maxpv; else v=zeros(kmax_leg,1); end
    rM=v.*e_rel.*rkin;
    rg=sum(Y_BM.*v.*r);
    diagV=diag(v);
    dxdt=SxZ*diagV*r*BIOM;
    dy=zeros(length(y),1);
    dy(1:length(x_ini))=dxdt(1:length(x_ini));
    dy(length(x_ini)+1:end)=alpha+re.*u-(beta+rg).*e;
    
    

function K=set_K(nz,x_ode,x_kin_all,K_all,nx_kin_all)


nx=length(x_ode);
K=zeros(nz,nx);

nx_kin_sum=zeros(nz,1);
sum=0; 
for i=2:nz
    sum=sum+nx_kin_all(i-1);
    nx_kin_sum(i)=sum;
end

for i=1:nz
    
	x_kin=x_kin_all(nx_kin_sum(i)+1:nx_kin_sum(i)+nx_kin_all(i));

    nx_kin=length(x_kin);
    for j=1:nx_kin
        ix_kin(j)=strmatch(x_kin(j),x_ode,'exact');
        K(i,ix_kin(j))=K_all(nx_kin_sum(i)+j);
    end
end

