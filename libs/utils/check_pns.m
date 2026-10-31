function  [ok, pns_norm, pns_comp, t_axis] = check_pns(seq,sys_type)
%WJG
if strcmp(sys_type,'prisma')
     [ok, pns_norm, pns_comp, t_axis] = seq.calcPNS('MP_GPA_K2309_2250V_951A_AS82.asc'); % 
elseif strcmp(sys_type,'skyra')
    [ok, pns_norm, pns_comp, t_axis] = seq.calcPNS('MP_GPA_K2309_2250V_793A_GC99.asc'); % 
elseif strcmp(sys_type,'Connectome2') || strcmp(sys_type,'C2_simulate_prisma')
    [ok, pns_norm, pns_comp, t_axis] = seq.calcPNS('MP_GradSys_P034_c_CX600.asc'); % 
elseif strcmp(sys_type,'terra')
    [ok, pns_norm, pns_comp, t_axis] = seq.calcPNS('MP_GPA_K2298_2250V_793A_SC72CD_EGA.asc'); % jian added
else
    warning('Unknown system :%s\nSetting asc filename to []', sys_type)
end
if (ok==1)
    disp('PNS check passed')
else
    disp('PNS check failed')
end