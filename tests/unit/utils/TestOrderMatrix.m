classdef TestOrderMatrix < matlab.unittest.TestCase
    methods (Test)
        function testIsotropic(tc)
            % Isotropic: only b_00 nonzero -> Q = (1/3)*I, off-diag = 0
            Lmax = 24;
            N = idx(Lmax, Lmax, Lmax);
            psi = zeros(N,1);
            psi(idx(0,0,Lmax)) = 1/sqrt(4*pi);
            Q = order_matrix(psi, 'type', 'xz');
            % Diagonal entries should all equal 1/3
            tc.verifyEqual(Q(1), Q(2), 'AbsTol', 1e-12);
            tc.verifyEqual(Q(2), Q(3), 'AbsTol', 1e-12);
            % Off-diagonal entries should vanish
            tc.verifyEqual(Q(4), 0, 'AbsTol', 1e-12);
            tc.verifyEqual(Q(5), 0, 'AbsTol', 1e-12);
            tc.verifyEqual(Q(6), 0, 'AbsTol', 1e-12);
        end

        function testTrace(tc)
            % Trace of Q = <ux^2>+<uy^2>+<uz^2> = 0 always
            Lmax = 24;
            N = idx(Lmax, Lmax, Lmax);
            psi = zeros(N,1);
            psi(idx(0,0,Lmax)) = 1/sqrt(4*pi);
            psi(idx(2,0,Lmax)) = 0.1;
            Q = order_matrix(psi, 'type', 'xz');
            tc.verifyEqual(Q(1)+Q(2)+Q(3), 0, 'AbsTol', 1e-12);
        end

        function testB20AlignedWithZ(tc)
            % Pure b_20 > 0 -> rods align with z -> qz2 largest
            Lmax = 24;
            N = idx(Lmax, Lmax, Lmax);
            psi = zeros(N,1);
            psi(idx(0,0,Lmax)) = 1/sqrt(4*pi);
            psi(idx(2,0,Lmax)) = 0.1;
            Q = order_matrix(psi, 'type', 'xz');
            tc.verifyGreaterThan(Q(3), Q(1));  % qz2 > qx2
            tc.verifyGreaterThan(Q(3), Q(2));  % qz2 > qy2
        end
        
        function testXYSymmetryUnderB22(tc)
            % b_22 breaks x/y symmetry but not trace
            Lmax = 24;
            N = idx(Lmax, Lmax, Lmax);
            psi = zeros(N,1);
            psi(idx(0,0,Lmax)) = 1/sqrt(4*pi);
            psi(idx(2,2,Lmax)) = 0.1;
            Q = order_matrix(psi, 'type', 'xz');
            tc.verifyNotEqual(Q(1), Q(2));  % qx2 != qy2
            tc.verifyEqual(Q(1)+Q(2)+Q(3), 0, 'AbsTol', 1e-12);
        end
    end
end