% calculate the area of gz blip of sms
function [kz_acq,deltakz] =  calc_Pattern_SMS_blip(MB,ETL,sliceSep)
kz_acq = zeros(1,ETL);
if (MB==1)
    for loopi = 1:ETL
        kz_acq(loopi) = 0;
    end
    deltakz = 0;
elseif(MB ==2)
    ind = [-1,1];
    for loopi = 1:ETL
        kz_acq(loopi) = ind((mod(loopi,MB)+1));
    end
    deltakz = (1/MB)/sliceSep;
elseif (MB == 3)
    %v1
    % ind = [-2,1,1];
    % for loopi = 1:ETL
    %     kz_acq(loopi) = ind((mod(loopi,MB)+1));
    % end
    % deltakz = (1/MB)/sliceSep;
    %v2
    ind = [-1,1];
    for loopi = 1:ETL
        kz_acq(loopi) = ind((mod(loopi,length(ind))+1));
    end
    deltakz = (1/2)/sliceSep;
elseif (MB == 4)
    ind = [2,-1,-1];
    for loopi = 1:ETL
        kz_acq(loopi) = ind((mod(loopi,length(ind))+1));
    end
    deltakz = (1/3)/sliceSep;
end