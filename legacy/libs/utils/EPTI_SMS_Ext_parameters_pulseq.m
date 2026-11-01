function parameters = EPTI_SMS_Ext_parameters_pulseq(sequence_name,parameters)


seq_obj = mr.Sequence();         % Create a new sequence object
seq_obj.read(sequence_name);

parameters.pe_lines = seq_obj.getDefinition('pe_lines'); 
parameters.PAT_factor = seq_obj.getDefinition('PAT_factor'); 
parameters.CNseg = seq_obj.getDefinition('CNseg'); 
if (isempty(seq_obj.getDefinition('TE')))
     parameters.alTE(1) =0;
else
    parameters.alTE(1) = seq_obj.getDefinition('TE')*1e6;
end
if (isempty(seq_obj.getDefinition('TE2')))
     parameters.alTE(2) =0;
else
    parameters.alTE(2) = seq_obj.getDefinition('TE2')*1e6;
end
parameters.SE_echoindex = seq_obj.getDefinition('center_echo');

grad_shape = seq_obj.getDefinition('TrapezoidGriddingParameters');
parameters.rampUpTime = grad_shape(1);
parameters.rampDownTime = grad_shape(3);
parameters.acqDelayTime = grad_shape(4);
parameters.flatTopTime = grad_shape(2);
parameters.iEffectiveEpiEchoSpacing =  (parameters.rampUpTime+parameters.rampDownTime+parameters.flatTopTime)*1e6;
parameters.dwellTime = seq_obj.adcRasterTime;
parameters.reconNx_ = seq_obj.getDefinition('TargetGriddedSamples'); 
parameters.resolution = seq_obj.getDefinition('Resolution'); 
parameters.DestNx = seq_obj.getDefinition('TargetGriddedSamples'); 
parameters.Rseg =  seq_obj.getDefinition('Rseg'); 
parameters.Nseg =  seq_obj.getDefinition('Nseg'); 
parameters.Rpe =  seq_obj.getDefinition('Rpe'); 
parameters.PF_shift =  seq_obj.getDefinition('PFshift'); 
parameters.Nyseg =  seq_obj.getDefinition('Nyseg'); 
parameters.Ryseg =  seq_obj.getDefinition('Ryseg'); 
parameters.Nzseg =  seq_obj.getDefinition('Nzseg'); 
parameters.Rzseg =  seq_obj.getDefinition('Rzseg'); 
parameters.N_phase_encode = parameters.Rseg *parameters.Nseg ; 

sampling_pattern_GE = seq_obj.getDefinition('pattern_GE')'; 
if isempty(sampling_pattern_GE)
    parameters.sampling_pattern_GE  = [];
else
    sampling_pattern_GE = reshape(sampling_pattern_GE,[],parameters.Nseg);
    parameters.sampling_pattern_GE =sampling_pattern_GE';
end
sampling_pattern_SE = seq_obj.getDefinition('pattern_SE')'; 
if isempty(sampling_pattern_SE)
    parameters.sampling_pattern_SE  = [];
else
    sampling_pattern_SE = reshape(sampling_pattern_SE,[],parameters.Nseg);
    parameters.sampling_pattern_SE =sampling_pattern_SE';
end
if (~isempty(seq_obj.getDefinition('ETL')))
    parameters.ETL = seq_obj.getDefinition('ETL');
end
parameters.GRE_line = seq_obj.getDefinition('GRE_Lines'); 
parameters.SE_line = seq_obj.getDefinition('SE_Lines'); 
if (isempty(parameters.GRE_line))
    parameters.GRE_line=0;
end
if (isempty(parameters.SE_line))
    parameters.SE_line=0;
end
parameters.echo_num = parameters.GRE_line+parameters.SE_line;
parameters.lPhaseEncodingLines = parameters.echo_num;
FOV = seq_obj.getDefinition('FOV')*1000; 
parameters.FOV_1 = FOV(1);
parameters.position = [0,0,0];
parameters.read_dir = [0,0,0];
parameters.N_slices = length(seq_obj.getDefinition('SlicePositions'));
parameters.nslice = length(seq_obj.getDefinition('SlicePositions'));
parameters.SMS = seq_obj.getDefinition('SMS'); 
parameters.SMS_shift = seq_obj.getDefinition('SMS'); 
parameters.TotalDuration = seq_obj.getDefinition('TotalDuration'); 
if isempty(parameters.SMS)
    parameters.SMS=1;
end
slice_order = seq_obj.getDefinition('actual_SlicePositions');

[~,slice_order] = sort(slice_order);
parameters.slice_order = slice_order(end:-1:1);
parameters.PF_shift  = seq_obj.getDefinition('PFshift'); 
if isempty(parameters.PF_shift )
    parameters.PF_shift = 0;
end
