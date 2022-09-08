function m_DofElem = f_ndof_elem(conec,m_gdl,conshyp,nelem)
if conshyp==1
    m_DofElem = reshape((m_gdl(conec',2:3))',16,nelem);
elseif conshyp==14
    m_DofElem = reshape((m_gdl(conec',2:end))',24,nelem);
end
