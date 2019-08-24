function Model=StoreResults(SPC,n_species,n_rct,GLU_in,F0,f,kd,V_colon,ini_stat,sim_time,length_each_tank,km,Pn,Km_water,density,T,Y,Y_Model,Y_met,met_co,lg_standard)
Model.SpeciesList = SPC;
Model.SpeciesNumber = n_species;
Model.RectorNumber=n_rct;
Model.FeedSubstrate=GLU_in;
Model.Feeding=F0;
Model.backflux=f;
Model.detachrate=kd;
Model.Volume=V_colon;
Model.InitCondition=ini_stat;
% Model.BMcoefficient=bm_co;
Model.TotalTime=sim_time;
Model.ylg=length_each_tank;
Model.km=km;
Model.Pn=Pn;   %% already include av
Model.Km_water=Km_water;
Model.rho=density;
Model.SimulationTime=T;
Model.SimulationResult=Y;
Model.Result_cell=Y_Model;
Model.SimulationMetabolites=Y_met;
Model.met_udf_co=met_co;
Model.length_tank=lg_standard;
