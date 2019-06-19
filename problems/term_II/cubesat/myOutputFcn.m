function status = myOutputFcn(t,q,flag,params)

if strcmp(flag,'init') && (t(1) == 0)
   params.RaRb = [];   
end

if ~strcmp(flag,'done') && ~strcmp(flag,'init')
   [~, Ra, Rb, zB, yA, zA, yB, Fp ] = dqdt(t, q, params);  
   fprintf('t=%6.3f s: Ra=%+5.3f N,  Rb=%+5.3f N\n',t, Ra, Rb);
end

% continue...
status = 0;

end

