function [ky_acq,ind] = calc_Pattern_WJG4(Nseg,Rseg,ETL,Rpe,Rpe2,PFshift)
ky_acq = zeros(Nseg,ETL);

if Nseg>1
    Start_point = (Rseg:Rseg:Nseg*Rseg) - Rseg/2;
else
    Start_point = (Rseg+PFshift)/2;
end
ceter_index = [];
skip = 1;
skip2 = 0;
for shot_num = 1:Nseg
    PE_currect = Start_point(shot_num);
    Sign_PEdirection = -1;
    for i = 1:ETL
            ky_acq(shot_num,i) = PE_currect;
      if i>2
        if(PE_currect<=Start_point+3)&&((PE_currect>=Start_point-3))&&(Sign_PEdirection==1)
            PE_currect = PE_currect+Sign_PEdirection*Rpe2;
            ceter_index = [ceter_index,i];
        elseif(PE_currect<=Start_point+4)&&((PE_currect>=Start_point-2))&&(Sign_PEdirection==-1)
            PE_currect = PE_currect+Sign_PEdirection*Rpe2;
            ceter_index = [ceter_index,i];
        elseif(PE_currect>=Start_point-9)&&(PE_currect<=Start_point-3)&&(Sign_PEdirection==1)&&(skip2==1)
            PE_currect=Start_point-Rpe2-1;
            skip2=0;
        elseif(PE_currect<=Start_point+10)&&(PE_currect>=Start_point+4)&&(Sign_PEdirection==-1)&&(skip2==1)
            PE_currect=Start_point+2*Rpe2;
            skip2=0;
        elseif(((PE_currect>=Rseg+PFshift-Rpe/2)||(PE_currect<=Rpe/2+PFshift))&&(skip==1)) 
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
            skip2=1;
        end
      end
    end
end
ind = find(diff(diff(ceter_index))<0);

ind = ceter_index(ind+1);
ind(4) = ind(3);
ind(3) = ind(2);
