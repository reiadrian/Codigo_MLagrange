% *************************************************************************
% MULTISCALE PROGRAM - MAIN FILE
% *************************************************************************
clear all; close all
main_example_path = pwd ;%'C:\Users\Hp 15-K118\Desktop\Projects\Multiscale_Application';%/home/javiermro/Projects/Examples'; 
%% ********* TESTS **********
first_mode=1;   nmode = 1;   Macro0 = cell(nmode,1);  
TEST_DATA(1).path_file= [main_example_path '/Examples/Articulo2/RVE64_4E/']; TEST_DATA(1).file = 'Macro.mfl' ;
isMICRO.MICRO =0; % For MS analysis 
TEST_DATA(1).nLab = 0; % Parallel pools
Snaps=[];

for iTEST=first_mode:nmode
    clc
    disp(['*** TEST NUMBER: ' num2str(iTEST) ' ***']);
    isMICRO.epsilon_Macro0 = Macro0{iTEST};
    analysis(TEST_DATA(iTEST).path_file,TEST_DATA(iTEST).file,TEST_DATA(iTEST).nLab,isMICRO,Snaps);
end

