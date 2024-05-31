function obj = lr_target(h, cholirf, tfp_pos)

    h1 = [0; h];    % Impose that the initial response (C0Sh = ISh = Sh) delivers a zero response of TFP
    h1 = h1./norm(h1); % Impose that h has unit length

    lr = sum(cholirf,3)*h1;     % The LR response is C(1)Sh
    obj = -lr(tfp_pos);         % The routine is for minimization and we want a maximum

end