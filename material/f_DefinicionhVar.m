function hVarNew = f_DefinicionhVar(conshyp,sihvarpg,nPG,protype)

   %Esta funcion inicializa las variables historicas segun el modelo constitutivo, considerando que el caso
   %multiescala puede ser una array de estructura y para el estandar una matriz (es para llamarse dentro de la
   %funcion del elemento). Esta generacion evita transferencia de datos a los nodos (habria que ver el tiempo
   %adicional agregado por el llamado de esta funcion).
   %sihvarpg = e_DatMatSet.sihvarpg;
   %nPG = e_DatElemSet.npg;
   %conshyp = e_DatMatSet.conshyp;   
   switch conshyp
       case {1,2,4,5,8,10,11,12,13,14,15,16,52,100,110} % AA: add case 14
         hVarNew = zeros(sihvarpg,nPG);
       case {50,55}
           if protype==0
               hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],'e_VarAux',[],...
                   'm_LinCond',[],'doff',[],'dofl',[],'c_DefMacro',[]);
           elseif protype==1
               hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],'e_VarAux',[],...
                   'm_LinCond',[],'doff',[],'dofl',[],'c_DefMacro',[],'c_GradPorMacro',[],'c_PorMacro',[],...
                   'm_VarFluc_eps0',[],'m_VarFluc_p0',[],'m_VarFluc_phi0',[]);
           end
       case 51
         hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],'e_VarAux',[],...
            'm_LinCond',[],'doff',[],'dofl',[],'m_ElemLoc',[],'c_DefMacro',[],'omegaMicroL',[],...
            'lMacro',[],'lMicro',[],'c_NormalesMicro',[],'longFis',[],'facNormMicro',[]);        %,'m_TensProy',[]
       case 53
         hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],'e_VarAux',[],...
            'm_LinCond',[],'doff',[],'dofl',[],'c_DefMacro',[]);        %,'m_TensProy',[]
       case 54
         hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'c_GdlCond',[],'Fint',[],...
             'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],...
             'doff',[],'dofl',[],'m_ElemLoc',[],'c_DefMacro',[],...
             'omegaMicroL',[],'lMacro',[],'lMicro',[],'c_NormalesMicro',[],...
             'longFis',[],'facNormMicro',[] );        %,'m_TensProy',[]        
        otherwise
         error('Matrices Elementales: Variables Historicas: Inicializacion: Modelo constitutivo no definido.')         
   end
   
end