function [ky_acq] = calc_Pattern_WJG5(Nseg,Rseg,ETL,Rpe,PFshift)
ky_acq = zeros(Nseg,ETL);

if Nseg>1
    Start_point = (Rseg:Rseg:Nseg*Rseg) - Rseg/2;
else
    Start_point = (Rseg+PFshift)/2;
end
skip = 1;
for shot_num = 1:Nseg
    PE_currect = Start_point(shot_num);
    Sign_PEdirection = -1;
    for i = 1:ETL
            ky_acq(shot_num,i) = PE_currect;
      if i>2
        if(((PE_currect>=Rseg+PFshift-Rpe/2)||(PE_currect<=Rpe/2+PFshift))&&(skip==1)) 
            Sign_PEdirection = Sign_PEdirection*-1;
            PE_currect = PE_currect+Sign_PEdirection*Rpe/2;
            skip=0;
        elseif(((PE_currect>Rseg+PFshift-Rpe)||(PE_currect<=Rpe+PFshift))&&(skip==1)) 
            Sign_PEdirection = Sign_PEdirection*-1;
            PE_currect = PE_currect+Sign_PEdirection*Rpe/2;
            skip=0;
        else
            PE_currect = PE_currect+Sign_PEdirection*Rpe;
            skip=1;
        end
      end
    end
end
