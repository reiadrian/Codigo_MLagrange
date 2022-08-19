function [ct,BiotM,beta,PermK,sigmaE_new,sigmaT_new,mflu_new,velflu_new] = ...
    f_Rmap_Bif_ElastME(eps_new,phi_new,p_new,N4,e_DatMatSet)
   
   ct = e_DatMatSet.ce;
   BiotM = e_DatMatSet.m_Biot;
   beta = e_DatMatSet.beta;
   PermK = e_DatMatSet.m_PermK;
   
   sigmaE_new = ct*eps_new; %Tensiones efectivas
   sigmaT_new = sigmaE_new - BiotM*p_new; % Tensiones totales
   mflu_new = BiotM'*eps_new + beta*p_new; % Contenido de masa del fluido
   velflu_new = PermK*phi_new; % Velocidad de filtracion del fluido
%    velflu_new = -PermK*phi_new; %VER PORQUE NO ES MENOS (-)
%    velflu_new = -PermK*DerivN*phi_new; % Velocidad de filtracion del fluido
%    mflu_new = BiotM'*M_Id*B*ud;
   
   

   
   