function [Fext] = f_Fflux4(Fext,u,e_DatSet,e_VG,c_Cw_eps,c_bw_p,c_kw_phi)
% Determina las fuerzas debidado a las poropresiones en cada paso de tiempo
% y los suma al vector global de fuerzas externas

nSet = e_VG.nSet;
Dtime = e_VG.Dtime;

c_Fext = cell(nSet,1);
c_FilFza = cell(nSet,1);

for iSet = 1:nSet 
    e_DatElemSet = e_DatSet(iSet).e_DatElem;
    m_DofElem = e_DatSet(iSet).m_DofElem;
    nElem = e_DatSet(iSet).nElem;
    dofpe = e_DatElemSet.dofpe;
    dofpe_p = e_DatElemSet.dofpe_p;
    dofpe_d = e_DatElemSet.dofpe_d;
    wg = e_DatElemSet.wg;
    nPG = e_DatElemSet.npg;
    pos_p = e_DatElemSet.pos_p;
    pos_d = e_DatElemSet.pos_d;
    pos_p0 = e_DatElemSet.pos_p0;
    
    m_BT_d = e_DatSet(iSet).m_BT_d;
    m_FF_p = e_DatSet(iSet).m_FF_p;
    m_DerCa_p = e_DatSet(iSet).m_DerCa_p;
    m_DetJT_p = e_DatSet(iSet).m_DetJT_p;

    m_Cw_eps= c_Cw_eps{iSet};
    m_bw_p= c_bw_p{iSet};
    m_kw_phi= c_kw_phi{iSet};
    
    % Grados de libertad y coordenadas de los nodos de los elementos del set
    dofElemSet = m_DofElem(:); 
    m_FilFza = dofElemSet';
    uElemSet  = reshape(u(dofElemSet),[],nElem);
    
    m_Fext = zeros(dofpe,nElem);
    for iElem = 1:nElem
        Cw_epsE = m_Cw_eps(:,:,:,iElem);     
        kw_phiE = m_kw_phi(:,:,:,iElem);
        bw_pE = m_bw_p(:,:,:,iElem);

        Kpp = zeros(dofpe_p,dofpe_p);
        Qps = zeros(dofpe_p,dofpe_d);
        m_Dercae_p = m_DerCa_p(:,:,:,iElem);
        ue = uElemSet(:,iElem);
        ue_p = ue(pos_p);
        ue_d = ue(pos_d);
        
        m_pesoPG_p = m_DetJT_p(:,iElem).*wg;        
        m_Be_d=m_BT_d(:,:,:,iElem);
    for iPG = 1:nPG
        N4 = m_FF_p(:,:,iPG);
        B = m_Be_d(:,:,iPG); 
        DerivN = m_Dercae_p(:,:,iPG);
        
        Cw_eps = Cw_epsE(:,:,iPG);
        kw_phi = kw_phiE(:,:,iPG);
        bw_p = bw_pE(:,:,iPG);
        
        Qps = Qps + (-DerivN'*Cw_eps*B)*m_pesoPG_p(iPG);
        Kpp = Kpp + (-DerivN'*kw_phi*DerivN-DerivN'*bw_p*N4)*m_pesoPG_p(iPG); 
    end %for(iPG)

    m_Fext(pos_p,iElem) =  Dtime*Kpp*ue_p+Dtime*Qps*ue_d;
    m_Fext(pos_p0,iElem) =  0.0; %Es necesario?
    
    end %for(iElem)      
    c_Fext{iSet} = m_Fext(:);
    c_FilFza{iSet} = m_FilFza;
end %for(iSet)

Fext = Fext + sparse([c_FilFza{:}],1,cat(1,c_Fext{:}));

    
    

