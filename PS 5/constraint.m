function [c, ceq] = constraint(h)
    
    c = [];
    ceq = h'*h - 1;
    

end