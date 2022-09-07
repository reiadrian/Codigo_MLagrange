%Funcionque permite crear vector de fuerzas para modelos con DOS elementos
%de diferentes gdl. (Elemento para medio poroso saturado y para medio
%solido)
function [ndoft,ndof_pm,ndof_sm,conec_pm,conec_sm] = f_ndof_int(e_DatSet,nSet,ndn_pm,ndn_sm)
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
ndoft = ndof_pm+ndof_sm;  