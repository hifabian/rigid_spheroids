function [gxz, wy, gyz, wx] = flow_parameters(Du)
    if isnumeric(Du) && ismatrix(Du)
        gxz = Du(3,1)+Du(1,3);
        wy  = -Du(1,3);
        gyz = Du(3,2)+Du(2,3);
        wx  = Du(2,3);

        assert(Du(1,1) == 0, "Non-zero u_{x,x} not supported!")
        assert(Du(1,2) == 0, "Non-zero u_{y,x} not supported!")
        assert(Du(2,1) == 0, "Non-zero u_{x,y} not supported!")
        assert(Du(2,2) == 0, "Non-zero u_{y,y} not supported!")
        assert(Du(3,3) == 0, "Non-zero u_{z,z} not supported!")
    else
        gxz = @(Du) Du(3,1)+Du(1,3);
        wy  = @(Du) -Du(1,3);
        gyz = @(Du) Du(3,2)-Du(2,3);
        wx  = @(Du) Du(2,3);
    end
end