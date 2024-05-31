function obj = lr_target_level(h, cholirf, tfp_pos)

    h = h./norm(h); % Impose that h has unit length
    h1 = [0; h];    % Impose that the initial response (C0Sh = ISh = Sh) delivers a zero response of TFP
    
    lr = squeeze(cholirf(:,:,end))*h1;     % The LR response is C(1)Sh
    obj = -lr(tfp_pos);         % The routine is for minimization and we want a maximum

end