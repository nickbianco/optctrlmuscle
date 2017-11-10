function [Fpe] = HillModel_PassiveMuscleFiber(lMtilde, e0, kpe)

t50 = exp(kpe .* (0.2 - 0.10e1) ./ e0);
pp1 = (t50 - 0.10e1); 
t7 = exp(kpe); 
pp2 = (t7 - 0.10e1);
t5 = exp(kpe .* (lMtilde - 0.10e1) ./ e0);
Fpe = ((t5 - 0.10e1) - pp1) ./ pp2;


end