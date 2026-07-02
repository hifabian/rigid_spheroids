function [f, Fjac] = assemble_unsteady(Du, Dr, beta, L2, Gxz, iLy, Gyz, iLx)
% Assembles the discretized linear operator with Lagrange multiplier

    L = -Dr*L2;
    Sxz = beta*Gxz+0.5*(1-beta)*iLy;
    Syz = beta*Gyz-0.5*(1-beta)*iLx;

    % Extract parameters from flow gradient
    [gxz, wy, gyz, wx] = flow_parameters(Du);

    % Assemble function and Jacobian
    f = @(t,y) L*y;
    Fjac = @(t,y) L;
    if isnumeric(Du) && ismatrix(Du)
        f = @(t,y) f(t,y) - gxz*(Sxz*y) - wy*(iLy*y) ...
                          - gyz*(Syz*y) - wx*(iLx*y);
        Fjac = @(t,y) Fjac(t,y) - gxz*Sxz - wy*iLy ...
                                - gyz*Syz - wx*iLx;
    else
        f = @(t,y) f(t,y) - gxz(Du(t))*(Sxz*y) - wy(Du(t))*(iLy*y) ...
                          - gyz(Du(t))*(Syz*y) - wx(Du(t))*(iLx*y);
        Fjac = @(t,y) Fjac(t,y) - gxz(Du(t))*Sxz - wy(Du(t))*iLy ...
                                - gyz(Du(t))*Syz - wx(Du(t))*iLx;
    end

end