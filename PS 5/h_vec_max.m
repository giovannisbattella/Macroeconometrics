function max_effect=h_vec_max(theta,D)

theta_13=theta(1);
theta_14=theta(2);

h=[0;-cos(theta_13)*cos(theta_14); -sin(theta_13)*cos(theta_14); -sin(theta_14)];

% Recall that we want the maximum response on the level, so we need to take
% cumsum
D_last=cumsum(D(:,:,:),3);

% The target function is C_40 
max_effect=squeeze(D_last(1,:,41))*h;


end