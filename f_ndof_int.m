%Funcionque permite crear vector de fuerzas para modelos con DOS elementos
%de diferentes gdl. (Elemento para medio poroso saturado y para medio
%solido)
function [ndoft,m_gdl,ndof_pm,ndof_sm,conec_pm,conec_sm] = f_ndof_int(e_DatSet,nSet,ndn_pm,ndn_sm)
if nSet~=2
    error('Lectura de datos: Cargas aplicadas: Tipo de carga no definida')
end
conshyp1 = e_DatSet(1).e_DatMat.conshyp;
    if conshyp1 == 14 
        conec_pm = unique(e_DatSet(1).conec');%Traspongo por si es un unico elemento y quede en formato martriz columna
        conec_sm = unique(e_DatSet(2).conec');
    else
        conec_sm = unique(e_DatSet(1).conec');
        conec_pm = unique(e_DatSet(2).conec');
    end

n_pm = size(conec_pm,1);
for kf = 1:n_pm
    pk = find(conec_sm==conec_pm(kf));
    if ~isempty(pk)
        conec_sm(pk) = [];
    end
end
n_sm = size(conec_sm,1);
ndof_pm = n_pm*ndn_pm;
ndof_sm = n_sm*ndn_sm;
ndoft = ndof_pm+ndof_sm;%Numero de gdl para problemas con elemntos con gdl variable
% m_con_f_sm = zeros(n_sm,4);
% m_con_f_pm = zeros(n_pm,4);
m_con_f_sm = [conec_sm ones(n_sm,1) ones(n_sm,1) zeros(n_sm,1)];
m_con_f_pm = [conec_pm ones(n_pm,1) ones(n_pm,1) ones(n_pm,1)];
m_gdl = [m_con_f_pm; m_con_f_sm];
m_gdl = sortrows(m_gdl,1);
nnodo = size(m_gdl,1);
pos = 0;
for inodo = 1:nnodo
    if m_gdl(inodo,4)==0
        for in_sm=1:ndn_sm
            pos = pos+1;
            m_gdl(inodo,in_sm+1)=pos;
        end
    else
         for in_pm=1:ndn_pm
            pos = pos+1;
            m_gdl(inodo,in_pm+1)=pos;
        end
    end
    
end

