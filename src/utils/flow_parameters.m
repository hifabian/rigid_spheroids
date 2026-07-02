function [gxz, wy, gyz, wx] = flow_parameters(Du)
    if isnumeric(Du) && ismatrix(Du)
        gxz = Du(3,1)+Du(1,3);
        wy  = -Du(1,3);
        gyz = Du(3,2)+Du(2,3);
        wx  = Du(2,3);

        % TODO warning if any other values are nonzero!
    else
        gxz = @(Du) Du(3,1)+Du(1,3);
        wy  = @(Du) -Du(1,3);
        gyz = @(Du) Du(3,2)-Du(2,3);
        wx  = @(Du) Du(2,3);
    end
end