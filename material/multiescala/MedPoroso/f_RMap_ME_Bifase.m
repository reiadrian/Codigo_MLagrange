function [e_TanOp,sigmaE_new,sigmaT_new,velflu_new,mflu_new,mfluY_new,hvar_newMacro] =...
    f_RMap_ME_Bifase(eps_new,phi_new,p_new,hvar_oldMacro,e_DatMatSetMacro,e_VGMacro)

%Se recupera variables micro
xx = e_DatMatSetMacro.xx;
omegaMicro_d = e_DatMatSetMacro.omegaMicro_d;

e_DatSet = e_DatMatSetMacro.e_DatSet;
e_VG = e_DatMatSetMacro.e_VG;
esImplexMacro = e_DatMatSetMacro.esImplex;
e_VarEst_old = hvar_oldMacro.e_VarEst;
u = hvar_oldMacro.u;
c_GdlCond = hvar_oldMacro.c_GdlCond;
Fint         = hvar_oldMacro.Fint;
m_LinCond    = hvar_oldMacro.m_LinCond;
doff         = hvar_oldMacro.doff;
dofl         = hvar_oldMacro.dofl;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%Para calculo de derivadas de variables fluctuantes respecto de
%varaibles macroescala
m_VarFluc_eps0 = hvar_oldMacro.m_VarFluc_eps0;
m_VarFluc_p0 = hvar_oldMacro.m_VarFluc_p0;
m_VarFluc_phi0 = hvar_oldMacro.m_VarFluc_phi0;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

e_VarAux = hvar_oldMacro.e_VarAux;

% VARIABLES GLOBALES
nElemTot = e_VG.nElem;
ntens = e_VG.ntens;
ndoft = e_VG.ndoft;
conshyp = e_VG.conshyp;

% INICIALIZACION DE VARIABLES
Du_step_old = zeros(ndoft,1);  % Vector de incrementos de desplazamientos en el paso de tiempo previo
Fext = zeros(ndoft,1);

%Por si es necesario el paso tiempo a nivel micro.
e_VG.istep = e_VGMacro.istep;
e_VG.Dtime = e_VGMacro.Dtime;
% Guardo conshyp MACRO para calcular velflu en f_MatElem_BifaseMulTSc
e_VG.conshypMacro = e_VGMacro.conshyp;
%Se guarda que se quiere que el modelo constitutivo se comporte elasticamente.
e_VG.elast = e_VGMacro.elast;
%Para impresion y debug se guarda la iteracion macro.
e_VG.iterMacro = e_VGMacro.iter;
%Para imprension se guarda el numero de elemento macro y numero de PG macro.
e_VG.iElemNumMacro = e_VGMacro.iElemNum;
e_VG.iPGMacro = e_VGMacro.iPG;
e_VG.SharedFS = e_VGMacro.SharedFS; 

%Se guarda los datos del matlabpool (se utiliza para imprimir de pasos y iteraciones a nivel micro que no
%convergieron).
e_VG.nLab = e_VGMacro.nLab; 
e_VG.tipoMPool = e_VGMacro.tipoMPool; 

% DESPLAZAMIENTO IMPUESTO
%No seria necesario estas operaciones porque vfix en todos los modelos clasicos de las
%formulaciones multiescala son nulos. Ademas que habria que interpretar como realizar el delta
%psi_value a nivel micro, ya deberia se corresponder con el incremento de tiempo a nivel macro,
%pero lo que implicaria que en cada paso de tiempo se esta resolviendo un RVE distinto.

% Deformacion macro aplicada en la Celda unitaria.
%Se aplica en forma uniforme en todo dominio.
%Deformacion macro por elemento y por punto de gauss, dividida en sets.
c_DefMacro = arrayfun(@(x)repmat(eps_new,[1,x.e_DatElem.npg,x.nElem]),e_DatSet,'UniformOutput',false); % AA:Desbloquee
c_GradPorMacro = arrayfun(@(x)repmat(phi_new,[1,x.e_DatElem.npg,x.nElem]),e_DatSet,'UniformOutput',false);
c_PorMacro = arrayfun(@(x)repmat(p_new,[1,x.e_DatElem.npg,x.nElem]),e_DatSet,'UniformOutput',false);
%Deformacion macro por elemento, dividad en sets.
%    c_DefMacro = arrayfun(@(x)repmat(eps_new,[1,x.nElem]),e_DatSet,'UniformOutput',false);
%    c_GradPorMacro = arrayfun(@(x)repmat(phi_new,[1,x.nElem]),e_DatSet,'UniformOutput',false);
%    c_PorMacro = arrayfun(@(x)repmat(p_new,[1,x.nElem]),e_DatSet,'UniformOutput',false);

% ESQUEMA DE NEWTON-RAPHSON
Du_step_new = zeros(ndoft,1);
if conshyp==15
    [Fext] = f_Fflux(Fext,u,e_DatSet,e_VG); %AA: add function
elseif conshyp==16
    [Fext] = f_Fflux_ML(Fext,u,e_DatSet,e_VG); %AA: add function
end

[u,c_GdlCond,Fint,e_VarEst_new,e_VarAux,Du_step_new,c_CT,KT] = newton_raphsonPM(...
xx,m_LinCond,dofl,doff,u,Du_step_new,c_GdlCond,Du_step_old,Fint,Fext,e_VarEst_old,e_VarAux,e_DatSet,...
c_DefMacro,c_GradPorMacro,c_PorMacro,e_VG);

% OPERADORES TANGENTES HOMOGENEIZADOS
[e_TanOp,m_VarFluc_eps0,m_VarFluc_p0,m_VarFluc_phi0] = f_OpTangHomBifase(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,omegaMicro_d,...
true(nElemTot,1),true(nElemTot,1),xx,e_VG,m_VarFluc_eps0,m_VarFluc_p0,m_VarFluc_phi0);

%Se asume que no se realiza analisis de bifurcacion con el tensor tangente constitutivo homogeneizado, por
%lo que en el caso ser implex, se devuelve nulo el tensor implicito homogeneizado.
% VEEEEER!!!!!!!!!!!!!!
if esImplexMacro
    m_CTHomog = struct('Implex',m_CTHomog,'Impli',zeros(ntens,ntens));
end

% VEEEEER!!!!!!!!!!!!!!   
% CALCULO DE VARIABLES HOMOGENEIZADAS
sigmaE_new = f_HomogArea({e_VarEst_new.sigmaE},ntens,omegaMicro_d,{e_DatSet.m_DetJT_d},e_DatSet,e_VG);
sigmaT_new = f_HomogArea({e_VarEst_new.sigmaT},ntens,omegaMicro_d,{e_DatSet.m_DetJT_d},e_DatSet,e_VG);
mflu_new = f_HomogArea({e_VarEst_new.mflu},1,omegaMicro_d,{e_DatSet.m_DetJT_d},e_DatSet,e_VG);

velflu_new = f_HomogArea({e_VarEst_new.velflu},2,omegaMicro_d,{e_DatSet.m_DetJT_d},e_DatSet,e_VG);
mfluY_new = f_HomogArea({e_VarEst_new.mYcord},2,omegaMicro_d,{e_DatSet.m_DetJT_d},e_DatSet,e_VG);
defHomogFl = f_HomogArea({e_VarEst_new.eps_fluct},ntens,omegaMicro_d,{e_DatSet.m_DetJT_d},e_DatSet,e_VG);

%Verificacion de la media del desplazamiento y de la poropresion
[m_uMedioFl_d,m_uMedioFl_p] = f_MediaDespCU_Bif(u,omegaMicro_d,e_DatSet,e_VG);

fprintf('Elemento %d: PG %d: Norma de la deformacion fluctuante media: %g\n',e_VGMacro.iElemNum,...
e_VGMacro.iPG,norm(defHomogFl))
fprintf('Elemento %d: PG %d: Norma del desplazamiento fluctuante medio: %g\n',e_VGMacro.iElemNum,...
e_VGMacro.iPG,norm(m_uMedioFl_d))
fprintf('Elemento %d: PG %d: Norma de la poropresion fluctuante media: %g\n',e_VGMacro.iElemNum,...
e_VGMacro.iPG,norm(m_uMedioFl_p)) 

hvar_newMacro = struct('u',u,'c_GdlCond',{c_GdlCond},'Fint',Fint,'e_VarEst',e_VarEst_new,...
'e_VarAux',e_VarAux,'m_LinCond',m_LinCond,'doff',doff,'dofl',dofl,'c_DefMacro',{c_DefMacro},...
'c_GradPorMacro',{c_GradPorMacro},'c_PorMacro',{c_PorMacro},...
'm_VarFluc_eps0',m_VarFluc_eps0,'m_VarFluc_p0',m_VarFluc_p0,'m_VarFluc_phi0',m_VarFluc_phi0);


end

%%
function [m_uMedio_d,m_uMedio_p] = f_MediaDespCU_Bif(u,omegaMicro,e_DatSet,e_VG)

   nSet = e_VG.nSet;
   %El desplazamiento es un campo nodal, por lo que para intregrar sobre el dominio se lo lleva a los punto de
   %gauss.
   c_DespElem = cell(nSet,1);
   c_PorpElem = cell(nSet,1);

   %
   for iSet = 1:nSet
      %
      nElem = e_DatSet(iSet).nElem;
      nPG = e_DatSet(iSet).e_DatElem.npg;
      m_DofElem = e_DatSet(iSet).m_DofElem;
      m_FF_d = e_DatSet(iSet).m_FF_d;
      m_FF_p = e_DatSet(iSet).m_FF_p;
      %
      pos_d =  e_DatSet(iSet).e_DatElem.pos_d;
      pos_p =  e_DatSet(iSet).e_DatElem.pos_p;
      m_uElemSet_d = reshape(u(m_DofElem(pos_d,1:nElem)),[],nElem);
      m_uElemSet_p = reshape(u(m_DofElem(pos_p,1:nElem)),[],nElem);
      %
      m_DespPG = zeros(2,nPG,nElem);
      m_PorpPG = zeros(1,nPG,nElem);
      for iElem = 1:nElem
         %squeeze llama reshape, asi que no es mas rapido que esta si se conoce cual es la dimension de la
         %matriz con valor 1.
         %Desplazamientos
         m_DespPG(:,:,iElem) = squeeze(sum(bsxfun(@times,m_FF_d(:,:,:),m_uElemSet_d(:,iElem)'),2));
         %Poropresiones
         m_PorpPG(:,:,iElem) = squeeze(sum(bsxfun(@times,m_FF_p(:,:,:),m_uElemSet_p(:,iElem)'),2));
      end
      %
      c_DespElem{iSet} = m_DespPG;
      c_PorpElem{iSet} = m_PorpPG;
   end
   
   m_uMedio_d = f_HomogArea(c_DespElem,2,omegaMicro,{e_DatSet.m_DetJT_d},e_DatSet,e_VG);
   m_uMedio_p = f_HomogArea(c_PorpElem,1,omegaMicro,{e_DatSet.m_DetJT_p},e_DatSet,e_VG);

end
