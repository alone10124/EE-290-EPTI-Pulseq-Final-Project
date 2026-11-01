function Matrix_encode = Generate_SMS_encode_Matrix2(ny,nz,ns,ns_sep)
% Generate the SMS encoding matrix for EPTI subspace reconstruction with
% SMS acquisition
% Fuyixue Wang, 2021, MGH 
% Jian Wu 2025
if ns == 1 || ns == 2
    PhaseShiftBase = 2*pi/ns_sep;
    Matrix_encode = ones(ny,nz,ns,1,1,1,ns_sep);
    for slice = 1:ns
        PhaseShift = rem(PhaseShiftBase*(slice-1), 2*pi);
        for count = 1:ns_sep
            Matrix_encode(:,:,slice,:,:,:,count) =  exp(1i*PhaseShift*(count-1));
        end
    end
elseif ns == 3
     PhaseShiftBase = 2*pi/2;
     Matrix_encode = ones(ny,nz,ns,1,1,1,2);
    for slice = 1:ns
        PhaseShift = rem(PhaseShiftBase*(slice-1), 2*pi);
        for count = 1:2
            Matrix_encode(:,:,slice,:,:,:,count) =  exp(1i*PhaseShift*(count-1));
        end
    end
 elseif ns == 4
     PhaseShiftBase = 2*pi/3;
     Matrix_encode = ones(ny,nz,ns,1,1,1,3);
    for slice = 1:ns
        PhaseShift = rem(PhaseShiftBase*(slice-1), 2*pi);
        for count = 1:3
            Matrix_encode(:,:,slice,:,:,:,count) =  exp(1i*PhaseShift*(count-1));
        end
    end
end

Matrix_encode = Matrix_encode(:,:,end:-1:1,:,:,:,end:-1:1);
% if ns_sep == 2
%     Matrix_encode(:,:,1,:,:,:,1) = exp(-1i*pi/2);
%     Matrix_encode(:,:,1,:,:,:,2) = exp(+1i*pi/2);
% elseif ns == 3
%     Matrix_encode(:,:,2,:,:,:,2) = exp(+1i*2*pi/3);
%     Matrix_encode(:,:,2,:,:,:,1) = exp(+1i*4*pi/3);
%     Matrix_encode(:,:,1,:,:,:,2) = exp(+1i*4*pi/3);
%     Matrix_encode(:,:,1,:,:,:,1) = exp(+1i*2*pi/3);
% end

end

