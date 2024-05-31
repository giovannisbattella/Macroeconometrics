function irf = partial_irf(cholirf,h)

    hor = size(cholirf,3);
    irf = zeros(size(cholirf,1), hor);

    for i = 1:hor
        irf(:,i) = cholirf(:,:,i)*h;
    end

end