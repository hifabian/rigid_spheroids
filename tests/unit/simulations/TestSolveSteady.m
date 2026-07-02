classdef TestSolveSteady < matlab.unittest.TestCase
    properties
        Lmax = 12;
        beta = 0.8;
    end

    methods (Test)

        function testZeroFlowIsIsotropic(tc)
            % No flow -> isotropic distribution -> Sy=Sz=Sx=0
            Du = [0,0,0;0,0,0;0,0,0];

            psi = solve_steady(tc.Lmax, tc.beta, Du, ...
                'store', false);
            Q = order_matrix(psi);
            [Sy, Sz, Sx] = order_parameters(Q);

            tc.verifyEqual(Sy, 0, 'AbsTol', 1e-10);
            tc.verifyEqual(Sz, 0, 'AbsTol', 1e-10);
            tc.verifyEqual(Sx, 0, 'AbsTol', 1e-10);
        end


        function testNormalization(tc)
            % b_00 = 1/sqrt(4*pi)
            Du = [0,0,0;0,0,0;1,1,0];

            psi = solve_steady(tc.Lmax, tc.beta, Du, ...
                'store', false);
            
            tc.verifyEqual(psi(idx(0,0,tc.Lmax)), 1/sqrt(4*pi), ...
                'RelTol', 1e-8, 'AbsTol', 1e-10);
        end


        function testXZvsYZRelfection(tc)
            % (y=-x)-reflection between xz and yz shear 
            Du_xz = [0,0,0;0,0,0;10.0,0,0];
            psi_xz = solve_steady(tc.Lmax, tc.beta, Du_xz, ...
                'store', false);
            
            Du_yz = [0,0,0;0,0,0;0,10.0,0];
            psi_yz = solve_steady(tc.Lmax, tc.beta, Du_yz, ...
                'store', false);

            N = size(psi_xz, 1); T = speye(N);
            for l = 0:2:tc.Lmax
                for m = 1:l
                    ii = idx(l,  m, tc.Lmax); jj = idx(l, -m, tc.Lmax);
                    c = cos(m*pi/2); s = sin(m*pi/2);
                    % 2x2 reflection block for (+m, -m) pair
                    T(ii,ii) = c; T(ii,jj) =  s;
                    T(jj,ii) = s; T(jj,jj) = -c;
                end
            end
            psi_yzr = T * psi_xz;

            tc.verifyEqual(psi_yzr, psi_yz, 'RelTol',1e-8,'AbsTol',1e-10);
        end


        function testXZvsMixedReflection(tc)
            % reflection along tan(theta)=(syz/sxy) and with
            % (sxz^2+syz^2)^0.5 shear-rate
            sxz = 100.0; syz = 20.0;
            smz = (sxz^2+syz^2)^0.5;
            theta = atan2(syz,sxz);

            Du_xz = [0,0,0;0,0,0;smz,0,0];
            psi_xz = solve_steady(tc.Lmax, tc.beta, Du_xz, ...
                'store', false);

            Du_mz = [0,0,0;0,0,0;sxz,syz,0];
            psi_mz = solve_steady(tc.Lmax, tc.beta, Du_mz, ...
                'store', false);
            
            N = size(psi_xz, 1); T = speye(N);
            for l = 0:2:tc.Lmax
                for m = 1:l
                    ii = idx(l,  m, tc.Lmax); jj = idx(l, -m, tc.Lmax);
                    c = cos(m*theta); s = sin(m*theta);
                    % 2x2 reflection block for (+m, -m) pair
                    T(ii,ii) = c; T(ii,jj) =  s;
                    T(jj,ii) = s; T(jj,jj) = -c;
                end
            end
            psi_mzr = T * psi_xz;

            tc.verifyEqual(psi_mzr, psi_mz, 'RelTol',1e-8,'AbsTol',1e-10);
        end


        function testXZvsYZRelfectionAdaptive(tc)
            % (y=-x)-reflection between xz and yz shear
            Du_xz = [0,0,0;0,0,0;10.0,0,0];
            psi_xz = solve_steady(64, tc.beta, Du_xz, ...
                'Ladaptive', true, 'store', false);

            Du_yz = [0,0,0;0,0,0;0,10.0,0];
            psi_yz = solve_steady(64, tc.beta, Du_yz, ...
                'Ladaptive', true, 'store', false);

            N = size(psi_xz, 1); T = speye(N); [Lmax, ~] = lmdx(N);
            for l = 0:2:Lmax
                for m = 1:l
                    ii = idx(l,  m, Lmax); jj = idx(l, -m, Lmax);
                    c = cos(m*pi/2); s = sin(m*pi/2);
                    % 2x2 reflection block for (+m, -m) pair
                    T(ii,ii) = c; T(ii,jj) =  s;
                    T(jj,ii) = s; T(jj,jj) = -c;
                end
            end
            psi_yzr = T * psi_xz;

            tc.verifyEqual(psi_yzr, psi_yz, 'RelTol',1e-8,'AbsTol',1e-10);
        end


        function testXZvsMixedReflectionAdaptive(tc)
            % reflection along tan(theta)=(syz/sxy) and with
            % (sxz^2+syz^2)^0.5 shear-rate
            sxz = 100.0; syz = 20.0;
            smz = (sxz^2+syz^2)^0.5;
            theta = atan2(syz,sxz);

            Du_xz = [0,0,0;0,0,0;smz,0,0];
            psi_xz = solve_steady(64, tc.beta, Du_xz, ...
                'Ladaptive', true, 'store', false);

            Du_mz = [0,0,0;0,0,0;sxz,syz,0];
            psi_mz = solve_steady(64, tc.beta, Du_mz, ...
                'Ladaptive', true, 'store', false);

            N = size(psi_xz, 1); T = speye(N); [Lmax, ~] = lmdx(N);
            for l = 0:2:Lmax
                for m = 1:l
                    ii = idx(l,  m, Lmax); jj = idx(l, -m, Lmax);
                    c = cos(m*theta); s = sin(m*theta);
                    % 2x2 reflection block for (+m, -m) pair
                    T(ii,ii) = c; T(ii,jj) =  s;
                    T(jj,ii) = s; T(jj,jj) = -c;
                end
            end
            psi_mzr = T * psi_xz;

            tc.verifyEqual(psi_mzr, psi_mz, 'RelTol',1e-8,'AbsTol',1e-10);
        end

        
        function testAdaptive(tc)
            % reflection along tan(theta)=(syz/sxy) and with
            % (sxz^2+syz^2)^0.5 shear-rate
            sxz = 100.0; syz = 20.0;
            Du = [0,0,0;0,0,0;sxz,syz,0];

            psi_direct = solve_steady(64, tc.beta, Du, ...
                'Ladaptive', false, 'store', false);
            psi_adapt = solve_steady(64, tc.beta, Du, ...
                'Ladaptive', true, 'store', false);

            Q_d = order_matrix(psi_direct, 'type', 'xz');
            Q_a = order_matrix(psi_adapt, 'type', 'xz');
            [Sy_d, Sz_d, Sx_d] = order_parameters(Q_d);
            [Sy_a, Sz_a, Sx_a] = order_parameters(Q_a);

            tc.verifyEqual(Sz_d, Sz_a, 'RelTol',1e-8,'AbsTol',1e-10);
            tc.verifyEqual(Sy_d, Sy_a, 'RelTol',1e-8,'AbsTol',1e-10);
            tc.verifyEqual(Sx_d, Sx_a, 'RelTol',1e-8,'AbsTol',1e-10);
        end

    end
end