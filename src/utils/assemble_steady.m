function A = assemble_steady(Du, Dr, beta, L2, Gxz, iLy, Gyz, iLx)
% Assembles the discretized linear operator with Lagrange multiplier

    % Extract parameters from flow gradient
    [gxz, wy, gyz, wx] = flow_parameters(Du);

    % For Lagrange multiplier
    c = sparse(1,1,(4*pi)^0.5, size(L2,1), 1);

    % Assemble matrix
    A = [0, c'; ...
         c, -Dr*L2 ...
            - gxz*(beta*Gxz+0.5*(1-beta)*iLy) - wy*iLy ...
            - gyz*(beta*Gyz-0.5*(1-beta)*iLx) - wx*iLx];
end