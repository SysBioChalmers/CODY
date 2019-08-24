%% this code is to simulate community culture of user-defined species in batch condition %%
function main_comm_batch

clc
clear
close all

currentDepth = 1; % get the supper path of the current path
currPath = fileparts(mfilename('fullpath'));% get current path
cd(currPath);


main_Co_Bth_Pure;
pause(10);
main_Co_Bad_Pure;
pause(10);
main_Co_Bth_Bad_Coculture;
