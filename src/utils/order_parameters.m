function [Sy, Sz, Sx, S, R, ExtChi, ExtTheta] = order_parameters(Q)
% Sy:       Projected order parameter in xz
% Sz:       Projected order parameter in xy
% Sx:       Projected order parameter in yz
% S:        Order parameter defined as l1-l3 (with eigenvalues l1>=l2>=l3)
%               (nematic orderinging)
% R:        Order parameter defined as l2-l3 
%               (biaxiality)
% ExtChi:   Azumithal angle of largest eigenvector
% ExtTheta: Polar angle of largest eigenvector
%
% Disclaimer: Much of this was written with help of Claude, though from my
% testing it appears correct. It does mainly some weird cos-arcos trick to
% avoid third roots, which could be numerically unstable.


    ibuf = repmat({':'}, 1, numel(size(Q))-1); % python >>> matlab
    Q1 = Q(ibuf{:},1); Q2 = Q(ibuf{:},2); Q3 = Q(ibuf{:},3);
    Q4 = Q(ibuf{:},4); Q5 = Q(ibuf{:},5); Q6 = Q(ibuf{:},6);

    % Projected order parameters
    Sy = sqrt((Q1-Q3).^2 + (2*Q6).^2);
    Sz = sqrt((Q1-Q2).^2 + (2*Q4).^2);
    Sx = sqrt((Q2-Q3).^2 + (2*Q5).^2);

    %% Exact 3x3 eigendecomposition for the leading eigenvector

    % traceless only
    trQ = Q1 + Q2 + Q3;
    Mxx = Q1 - trQ/3; Myy = Q2 - trQ/3; Mzz = Q3 - trQ/3;
    Mxy = Q4; Myz = Q5; Mzx = Q6;

    % Recall Cayley-Hamilton: A^3 - 0*A^2 + (0^2-tr(A^2))*A - det(A)*I = 0
    % --> A^3 + (-tr(A^2))*A + (-det(A))*I = 0
    % c.f. cubic equation: x^3 + p*t + q = 0  --> Cardano formula

    % largest eigenvalue via the Cardano formula:
    %   lam1 = u + v
    %   u = (-q/2 + D)^(1/3)    and     v = (-q/2 - D)^(1/3)
    %   D = Sqrt[ q^2/4 + p^3/27 ]
    trM2 = Mxx.^2 + Myy.^2 + Mzz.^2 + 2*(Mxy.^2 + Myz.^2 + Mzx.^2);
    p = sqrt(max(trM2,0)/6);

    detM = Mxx.*(Myy.*Mzz - Myz.^2) ...
         - Mxy.*(Mxy.*Mzz - Myz.*Mzx) ...
         + Mzx.*(Mxy.*Myz - Myy.*Mzx);

    pSafe = p; pSafe(p==0) = 1; % avoid 0/0 at isotropic points
    r = detM ./ (2*pSafe.^3);
    r = min(max(r,-1),1);       % clip numerical noise
    phi = acos(r)/3;            % --> r = cos(3*phi)

    lam1 = 2*p.*cos(phi);       % 
    lam1(p==0) = 0;
    lam3 = 2*p.*cos(phi+2*pi/3);
    lam2 = -lam1 - lam3;

    S = lam1-lam3;
    R = lam2-lam3;

    % eigenvector = null vector of (M - lam1*I), via cross product of two
    % rows; try all three row pairs and keep the largest-norm result
    % (robust against any single pair being parallel)
    Axx = Mxx - lam1; Ayy = Myy - lam1; Azz = Mzz - lam1;

    [vx,vy,vz]    = crossRows(Axx,Mxy,Mzx,  Mxy,Ayy,Myz);
    [vx2,vy2,vz2] = crossRows(Mxy,Ayy,Myz,  Mzx,Myz,Azz);
    [vx3,vy3,vz3] = crossRows(Mzx,Myz,Azz,  Axx,Mxy,Mzx);

    n1 = vx.^2+vy.^2+vz.^2; n2 = vx2.^2+vy2.^2+vz2.^2; n3 = vx3.^2+vy3.^2+vz3.^2;
    use2 = n2>n1 & n2>=n3; use3 = n3>n1 & n3>n2;
    vx(use2)=vx2(use2); vy(use2)=vy2(use2); vz(use2)=vz2(use2);
    vx(use3)=vx3(use3); vy(use3)=vy3(use3); vz(use3)=vz3(use3);

    vnorm = sqrt(vx.^2+vy.^2+vz.^2);
    vnorm(vnorm==0) = 1;                    % isotropic points: direction undefined
    vx = vx./vnorm; vy = vy./vnorm; vz = vz./vnorm;

    % fix n <-> -n gauge freedom (director symmetry): pick nz >= 0
    flip = vz < 0;
    vx(flip) = -vx(flip); vy(flip) = -vy(flip); vz(flip) = -vz(flip);

    ExtTheta = acos(min(max(vz,-1),1));   % polar angle from +z, in [0, pi]
    ExtChi   = atan2(vy, vx);             % azimuthal angle in xy, in (-pi, pi]
end

function [cx,cy,cz] = crossRows(a1,a2,a3, b1,b2,b3)
    cx = a2.*b3 - a3.*b2;
    cy = a3.*b1 - a1.*b3;
    cz = a1.*b2 - a2.*b1;
end