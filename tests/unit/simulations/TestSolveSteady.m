classdef TestSolveSteady < matlab.unittest.TestCase
    methods (Test)
        function testZeroFlowIsIsotropic(tc)
            % No flow -> isotropic distribution -> Sy=Sz=Sx=0
            psi = solve_steady(4, 0.5, 0, 0);
            Q = order_matrix(psi, 'type', 'xz');
            [Sy, Sz, Sx] = order_parameters(Q);
            tc.verifyEqual(Sy, 0, 'AbsTol', 1e-10);
            tc.verifyEqual(Sz, 0, 'AbsTol', 1e-10);
            tc.verifyEqual(Sx, 0, 'AbsTol', 1e-10);
        end

        function testNormalization(tc)
            % psi integrates to 1: b_00 = 1/sqrt(4*pi)
            Lmax = 24;
            psi = solve_steady(Lmax, 0.5, 1, 1);
            tc.verifyEqual(psi(idx(0,0,Lmax)), 1/sqrt(4*pi), 'AbsTol', 1e-10);
        end

        function testXZvsYZOrderParamter(tc)
            % solve_steady(sxz, 0) and solve_steady(0, syz) should give
            % same order parameter magnitude by reflection symmetry
            psi_xz = solve_steady(24, 0.5, 1.0, 0.0);
            psi_yz = solve_steady(24, 0.5, 0.0, 1.0);
            Q_xz = order_matrix(psi_xz, 'type', 'xz');
            Q_yz = order_matrix(psi_yz, 'type', 'xz');
            [Sy_xz, Sz_xz, Sx_xz] = order_parameters(Q_xz);
            [Sy_yz, Sz_yz, Sx_yz] = order_parameters(Q_yz);
            tc.verifyEqual(Sz_xz, Sz_yz, 'AbsTol', 1e-10);
            tc.verifyEqual(Sy_xz, Sx_yz, 'AbsTol', 1e-10);
            tc.verifyEqual(Sx_xz, Sy_yz, 'AbsTol', 1e-10);
        end

        function testXZvsYZRelfection(tc)
            % (y=-x)-reflection between xz and yz shear 
            Lmax = 24;
            psi_xz = solve_steady(Lmax, 0.5, 10.0, 0.0);
            psi_yz = solve_steady(Lmax, 0.5, 0.0, 10.0);
            N = size(psi_xz, 1);
            T = speye(N);
            for l = 0:2:Lmax
                for m = 1:l
                    ii = idx(l,  m, Lmax);
                    jj = idx(l, -m, Lmax);
                    c = cos(m*pi/2); s = sin(m*pi/2);
                    % 2x2 reflection block for (+m, -m) pair
                    T(ii,ii) = c; T(ii,jj) =  s;
                    T(jj,ii) = s; T(jj,jj) = -c;
                end
            end
            psi_yz_rotated = T * psi_xz;
            tc.verifyEqual(psi_yz_rotated, psi_yz, 'AbsTol', 1e-10);
        end

        function testXZvsMixedReflection(tc)
            % reflection along tan(theta)=(syz/sxy) and with
            % (sxz^2+syz^2)^0.5 shear-rate
            sxz = 100.0; syz = 20.0;
            theta = atan2(syz,sxz);
            Lmax = 24;
            psi_xz = solve_steady(Lmax, 0.5, (sxz^2+syz^2)^0.5, 0.0);
            psi_mz = solve_steady(Lmax, 0.5, sxz, syz);
            N = size(psi_xz, 1);
            T = speye(N);
            for l = 0:2:Lmax
                for m = 1:l
                    ii = idx(l,  m, Lmax);
                    jj = idx(l, -m, Lmax);
                    c = cos(m*theta); s = sin(m*theta);
                    % 2x2 reflection block for (+m, -m) pair
                    T(ii,ii) = c; T(ii,jj) =  s;
                    T(jj,ii) = s; T(jj,jj) = -c;
                end
            end
            psi_mz_rotated = T * psi_xz;
            tc.verifyEqual(psi_mz_rotated, psi_mz, 'AbsTol', 1e-10);
        end
    end
end