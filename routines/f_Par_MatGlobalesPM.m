%Esta funcion permite calcular la matriz de rigidez y el vector de fuerzas
%internas de cada elemento y global
function [KT,Fint,c_GdlCond,e_VarEst_new,e_VarAux,c_TensorTang,o_Par,...
    c_Cw_eps,c_bw_p,c_kw_phi] = f_Par_MatGlobalesPM(xx,u,du,...
    c_GdlCond,e_DatSet,e_VarEst_new,e_VarEst_old,e_VarAux,DefMacro,GradPorMacro,PorMacro,e_VG)

o_Par = [];
%Para no transferir estas variables que no se utilizan dentro del parfor.
e_VG.smooth_alpha = [];
e_VG.smooth_dalpha = [];
e_VG.smooth_dalpha_old = [];
e_VG.MGlobInv = [];

% Recupera variables globales
ntens = e_VG.ntens;
nSet = e_VG.nSet;
protype = e_VG.protype;
% Inicializaciones
c_Ke = cell(nSet,1);
c_Fint = cell(nSet,1);
c_Fil = cell(nSet,1);
c_Col = cell(nSet,1);
c_FilFza = cell(nSet,1);
c_TensorTang = cell(nSet,1);
c_Csig_eps = cell(nSet,1);
c_bw_eps= cell(nSet,1);
c_Cw_eps3= cell(nSet,1);
c_bsig_p= cell(nSet,1);
c_Mww_p= cell(nSet,1);
c_bw_p3= cell(nSet,1);
c_ksig_phi= cell(nSet,1);
c_kww_phi= cell(nSet,1);
c_kw_phi3= cell(nSet,1);
%
c_Cw_eps= cell(nSet,1);
c_bw_p= cell(nSet,1);
c_kw_phi= cell(nSet,1);
%Loop sobre el set de elementos (tipo y material del elemento)
for iSet = 1:nSet
    % Recuperacion de variables
    %(evita llamar desde las estructura de los sets en el bucle de los elementos, ver que solo
    %recupera punteros)
    nElem = e_DatSet(iSet).nElem;
    m_DofElem = e_DatSet(iSet).m_DofElem;
    eps_old = e_VarEst_old(iSet).eps;
    hvar_old = e_VarEst_old(iSet).hvar;
    aux_var = e_VarAux(iSet).VarAuxGP;
    m_VarAuxElem = e_VarAux(iSet).VarAuxElem;
    e_DatElemSet = e_DatSet(iSet).e_DatElem;
    dofpe = e_DatElemSet.dofpe;
    npg = e_DatElemSet.npg;
    e_DatMatSet = e_DatSet(iSet).e_DatMat;
    m_NumElem = e_DatSet(iSet).m_NumElem;
    m_DefMacroSet = DefMacro{iSet}; %AA: La def macro es cte en todo el elemento aplicada en cada PG
    m_GradPorMacroSet = GradPorMacro{iSet}; 
    m_PorMacroSet = PorMacro{iSet}; 
    conshyp = e_DatMatSet.conshyp;
    if conshyp==14
        m_BT_d = e_DatSet(iSet).m_BT_d; %Matriz B de desplazamientos-deformaciones
        m_DetJT_d = e_DatSet(iSet).m_DetJT_d; %Determinante del Jacobiano en desplazamientos
        m_DerCa_p = e_DatSet(iSet).m_DerCa_p; %Matriz en derivadas cartesianas para poropresiones
        m_DetJT_p = e_DatSet(iSet).m_DetJT_p; %Determinante del Jacobiano en poropresiones
    elseif conshyp==1
        m_BT = e_DatSet(iSet).m_BT; %Matriz B de desplazamientos-deformaciones
        m_DetJT = e_DatSet(iSet).m_DetJT; %Determinante del Jacobiano en desplazamientos
    end

    sigmaE_old = e_VarEst_old(iSet).sigmaE; %Tension efectiva previa
    sigmaT_old = e_VarEst_old(iSet).sigmaT; %Tension total previa

    sigmaE_new   = e_VarEst_new(iSet).sigmaE; %Tension efectiva actual
    sigmaT_new   = e_VarEst_new(iSet).sigmaT; %Tension total actual
    eps_new     = e_VarEst_new(iSet).eps; %Deformaciones actuales
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    velflu_new     = e_VarEst_new(iSet).velflu; %Velocidad de filtracion actual
    mfluY_new     = e_VarEst_new(iSet).mYcord; %VER DE ELIMINAR
    phi_new     = e_VarEst_new(iSet).phi;  %Gradiente de porpresiones actuales
    mflu_new     = e_VarEst_new(iSet).mflu; %Caudal por unidad de masa actuales 
    porpr_new     = e_VarEst_new(iSet).porpr; %Poropresiones actuales en PG 
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    eps_fluct   = e_VarEst_new(iSet).eps_fluct; %Deformaciones fluctuantes 
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    phi_fluct   = e_VarEst_new(iSet).phi_fluct; %Gradiente de poropresiones fluctuantes
    p_fluct   = e_VarEst_new(iSet).p_fluct; %Poropresiones fluctuantes en cada PG
    p_M   = e_VarEst_new(iSet).p_M; %Poropresiones macroescala en cada PG
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    hvar_new    = e_VarEst_new(iSet).hvar;           
    m_VarHistElemNew = e_VarEst_new(iSet).VarHistElem;

    % Inicializaciones
    %La verificacion si esImplex es field de la estructura no seria necesario si a todos los modelos
    %constitutivos se agrega este campo.
    if isfield(e_DatMatSet,'esImplex') && e_DatMatSet.esImplex 
        %No se puede disminuir la matriz m_TensorTang cuando se calcula 
        m_TensorTang = zeros(ntens,ntens,2*npg,nElem);
        %esImplex = 1;
    else
        m_TensorTang = zeros(ntens,ntens,npg,nElem);
        m_Csig_eps= zeros(ntens,ntens,npg,nElem);
        m_bw_eps= zeros(1,ntens,npg,nElem);
        m_Cw_eps3= zeros(2,ntens,npg,nElem);
        m_bsig_p= zeros(ntens,1,npg,nElem);
        m_Mww_p= zeros(1,1,npg,nElem);
        m_bw_p3=zeros(2,1,npg,nElem);
        m_ksig_phi= zeros(ntens,2,npg,nElem);
        m_kww_phi= zeros(1,2,npg,nElem);
        m_kw_phi3= zeros(2,2,npg,nElem);
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        m_Cw_eps= zeros(2,ntens,npg,nElem);
        m_bw_p=zeros(2,1,npg,nElem);
        m_kw_phi= zeros(2,2,npg,nElem);
        %esImplex = 0;
     end
    m_Ke = zeros(dofpe,dofpe,nElem);
    m_Fint = zeros(dofpe,nElem);   

    % Grados de libertad y coordenadas de los nodos de los elementos del set
    dofElemSet = m_DofElem(:);
    m_FilFza = dofElemSet';
    m_Fil = reshape(repmat(reshape(dofElemSet,dofpe,[]),dofpe,1),1,[]);
    m_Col = reshape(repmat(m_FilFza,dofpe,1),1,[]);
    uElemSet  = reshape(u(dofElemSet),[],nElem);
    duElemSet = reshape(du(dofElemSet),[],nElem);

    %Solo para debug, se guarda el numero de set en e_VG (ya que esta variable se pasa a todas las
    %funciones).
    e_VG.iSet = iSet;

    %Se podria poner el parfor de dentro de cada eltype, asi evitar hacer la verificacion de tipo
    %de elemento para cada set.
    %parfor iElem = 1:nElem
    %for iElem = 1:nElem

    %Esta linea no puede estar si esta el parfor activado.
    %e_VG.iElem = iElem;
    %No se puede modificar una variable definida previa al loop, si no es sliced y si no es
    %interpretada como de reduccion, ya que no sabe como interpretarla el MatLab (puede haber
    %superposicion de resultados, al hacer reduccion). Por eso cada Lab debe tener su "copia",
    %para modificarla en forma independiente.
    %Lo que puede ser lento es copiar para cada lab esa e_VG_Aux, que puede ser grande.
    %En realidad como cada procesador tiene su copia local del e_VG (ya que MatLab realiza una
    %copia por cada Lab de todas la variables, excepto las sliced, donde solo copia la parte
    %que le corresponde al procesador), y al "copiarse" esta al e_VG_Aux, lo que unico que se
    %hace es copiarse el puntero, ya que no se esta modificando el e_VG_Aux. Luego lo que se
    %realiza es la modificacion de un valor de un campo (field), donde supuestamente
    %MatLab al cambiar un field de una estructura no hace copia de toda la estructura de nuevo,
    %si solo del field, por lo tanto no tendria que ser mucho mas lenta.
    %Otra seria pasar la variable iElem como argumento en las funciones que llama dentro del
    %parfor.
    switch conshyp
        case {1}
            m_VolElem = e_DatSet(iSet).m_VolElem;
            %                  parfor iElem = 1:nElem
            for iElem = 1:nElem
                e_VG_Aux4 = e_VG;
                e_VG_Aux4.iElemSet = iElem;
                e_VG_Aux4.iElemNum = m_NumElem(iElem);
                %
                [m_Ke(:,:,iElem),m_Fint(:,iElem),sigmaT_new(:,iElem),eps_new(:,iElem),...
                eps_fluct(:,iElem),hvar_new(:,iElem),aux_var(:,iElem),...
                m_TensorTang(:,:,:,iElem)] =...
                f_MatElem_quad_q1(...
                uElemSet(:,iElem),eps_old(:,iElem), hvar_old(:,iElem),aux_var(:,iElem),...
                e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
                m_DefMacroSet(:,iElem),sigmaT_old(:,iElem),...    %m_DefMacroSet(:,:,iElem)
                m_VolElem(iElem),e_VG_Aux4);
            end
        case {14}
            % Fuerza interna y tensor tangente del elemento
            %AA: Cuadrangulo estandar de 8 nodos (Aplicado p/medio bifase)
            m_FF_p = e_DatSet(iSet).m_FF_p;
            conec = e_DatSet(iSet).conec;
            %                  parfor iElem = 1:nElem
            for iElem = 1:nElem
                e_VG_Aux16 = e_VG;
                e_VG_Aux16.iElemSet = iElem;
                e_VG_Aux16.iElemNum = m_NumElem(iElem);
                coord_n = f_CoordElem(xx,conec(iElem,:));
                [m_Ke(:,:,iElem),m_Fint(:,iElem),sigmaE_new(:,iElem),sigmaT_new(:,iElem),...
                eps_new(:,iElem),velflu_new(:,iElem),mflu_new(:,iElem),mfluY_new(:,iElem),eps_fluct(:,iElem),...
                phi_new(:,iElem),phi_fluct(:,iElem),porpr_new(:,iElem),p_fluct(:,iElem),p_M(:,iElem),...
                hvar_new(:,iElem),aux_var(:,iElem),m_Csig_eps(:,:,:,iElem),m_bw_eps(:,:,:,iElem),m_Cw_eps3(:,:,:,iElem),...
                m_bsig_p(:,:,:,iElem),m_Mww_p(:,:,:,iElem),m_bw_p3(:,:,:,iElem), m_ksig_phi(:,:,:,iElem),...
                m_kww_phi(:,:,:,iElem),m_kw_phi3(:,:,:,iElem),m_Cw_eps(:,:,:,iElem),...
                m_bw_p(:,:,:,iElem),m_kw_phi(:,:,:,iElem),m_TensorTang(:,:,:,iElem)] = f_MatElem_BifaseMulTSc(...
                uElemSet(:,iElem),eps_old(:,iElem),coord_n,hvar_old(:,iElem),aux_var(:,iElem),...
                e_DatElemSet,e_DatMatSet,m_BT_d(:,:,:,iElem),m_DetJT_d(:,iElem),...
                m_DerCa_p(:,:,:,iElem),m_DetJT_p(:,iElem),m_FF_p,m_DefMacroSet(:,iElem),...  %m_DefMacroSet(:,:,iElem)
                m_GradPorMacroSet(:,iElem),m_PorMacroSet(:,iElem),...
                sigmaE_old(:,iElem),sigmaT_old(:,iElem),e_VG_Aux16);
            end %for(iElem)
    end

    e_VarEst_new(iSet).sigmaE = sigmaE_new;
    e_VarEst_new(iSet).sigmaT = sigmaT_new;
    e_VarEst_new(iSet).eps = eps_new;
    e_VarEst_new(iSet).velflu = velflu_new;
    e_VarEst_new(iSet).phi = phi_new;
    e_VarEst_new(iSet).mYcord = mfluY_new;
    e_VarEst_new(iSet).mflu = mflu_new;
    e_VarEst_new(iSet).porpr = porpr_new;
    e_VarEst_new(iSet).eps_fluct = eps_fluct;
    e_VarEst_new(iSet).phi_fluct = phi_fluct;
    e_VarEst_new(iSet).p_fluct = p_fluct;
    e_VarEst_new(iSet).p_M = p_M;
    e_VarEst_new(iSet).hvar = hvar_new;
    e_VarEst_new(iSet).VarHistElem = m_VarHistElemNew;
    e_VarAux(iSet).VarAuxGP = aux_var;
    e_VarAux(iSet).VarAuxElem = m_VarAuxElem;
    c_Ke{iSet} = m_Ke(:);
    c_Fint{iSet} = m_Fint(:);
    c_Fil{iSet} = m_Fil;
    c_Col{iSet} = m_Col;
    c_FilFza{iSet} = m_FilFza;
%     if e_DatMatSet.conshyp ~= 50
%         %                m_TensorTang = e_TanOp.m_TensorTang;
%         c_TensorTang{iSet} = m_TensorTang;
%     elseif e_DatMatSet.conshyp == 50 
        % VER DE SACAR LOS RESTANTES OPERADORES TANGENTES!!
        c_TensorTang{iSet} = m_TensorTang;
        c_Csig_eps{iSet} = m_Csig_eps;
        c_bw_eps{iSet}= m_bw_eps;
        c_Cw_eps3{iSet}= m_Cw_eps3;
        c_bsig_p{iSet}= m_bsig_p;
        c_Mww_p{iSet}= m_Mww_p;
        c_bw_p3{iSet}= m_bw_p3;
        c_ksig_phi{iSet}= m_ksig_phi;
        c_kww_phi{iSet}= m_kww_phi;
        c_kw_phi3{iSet}= m_kw_phi3;
        
        c_Cw_eps{iSet}= m_Cw_eps;
        c_bw_p{iSet}= m_bw_p;
        c_kw_phi{iSet}= m_kw_phi;
%     end 
end %for(iSet)

% Ensamble de matriz de fuerzas internas global
Fint = sparse([c_FilFza{:}],1,cat(1,c_Fint{:}));

% Ensamble de matriz de rigidez global
KT = sparse([c_Fil{:}],[c_Col{:}],cat(1,c_Ke{:}));

end
