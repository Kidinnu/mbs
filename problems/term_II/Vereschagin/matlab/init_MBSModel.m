function Model=init_MBSModel(Model)
%
%
nq = sum(Model.na);

Model.iq  = cell(Model.n,1);
Model.idq = cell(Model.n,1);
istart=1;

for i=1:Model.n
    % Срезы для координат    
    Model.iq{i} =(istart:istart+Model.na(i)-1)';
    % Срезы для скоростей
    Model.idq{i}=Model.iq{i}+nq;    
    istart=istart+Model.na(i);
end
% Срез для дополнительных д.у.
if Model.n0 > 0
    Model.iauxq=((istart-1)*2+1:(istart-1)*2+Model.n0)';
else
    Model.iauxq=NaN;
end


