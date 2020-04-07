function [Value, isterminal, direction] = PhiEvent(t, q,p)
    % Если угол поврота превышает phimax...
    Value = q(6)-deg2rad(p.phimax);     
    % при увеличении угла, 
    direction = 1;                       
    % то остановить интегрирование
    isterminal = 1;                     
end

