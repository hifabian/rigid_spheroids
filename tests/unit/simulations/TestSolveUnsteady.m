classdef TestSolveUnsteady < matlab.unittest.TestCase
    properties
        Lmax = 8;
        beta = 0.5;
        Dr = 1.0;
        psi0;
    end

    methods (TestMethodSetup)
        function setup(tc)
            N = idx(tc.Lmax, tc.Lmax, tc.Lmax);
            tc.psi0 = zeros(N, 1);
            tc.psi0(idx(0,0,tc.Lmax)) = 1/sqrt(4*pi);
            tc.psi0(idx(2,0,tc.Lmax)) = 0.1;
            tc.psi0(idx(2,2,tc.Lmax)) = 0.05;
        end
    end

    methods (Test)
        function testInitialCondition(tc)
            % t=0 should recover psi0 exactly
            t = [0, 1];
            Du = [0,0,0;0,0,0;0,0,0];

            [~, psi] = solve_unsteady(t,tc.psi0,tc.Lmax,tc.beta,Du,tc.Dr);

            tc.verifyEqual(psi(1,:)', tc.psi0, 'AbsTol', 1e-12);
        end


        function testNormalizationConserved(tc)
            % b_00 = 1/sqrt(4pi) must be constant at all times
            t = linspace(0, 2, 20);
            Du = [0,0,0;0,0,0;0,0,0];

            [~, psi] = solve_unsteady(t,tc.psi0,tc.Lmax,tc.beta,Du,tc.Dr);
            b00 = psi(:, idx(0,0,tc.Lmax));

            tc.verifyEqual(b00, (1/sqrt(4*pi))*ones(size(b00)), ...
                'AbsTol', 1e-10);
        end


        function testZeroFlowDecaysToIsotropic(tc)
            % Without flow, all l>0 modes decay to zero
            t = linspace(0, 20/tc.Dr, 100);
            Du = [0,0,0;0,0,0;0,0,0];

            [~, psi] = solve_unsteady(t,tc.psi0,tc.Lmax,tc.beta,Du,tc.Dr);
            psi_end = psi(end,:)';

            psi_iso = zeros(size(psi_end));
            psi_iso(idx(0,0,tc.Lmax)) = 1/sqrt(4*pi);

            tc.verifyEqual(psi_end, psi_iso, 'AbsTol', 1e-8);
        end


        function testZeroFlowDecayRate(tc)
            % l=2 modes decay as exp(-l*(l+1)*Dr*t) = exp(-6*Dr*t)
            t = linspace(0, 2/tc.Dr, 20);
            Du = [0,0,0;0,0,0;0,0,0];
            
            [~, psi] = solve_unsteady(t,tc.psi0,tc.Lmax,tc.beta,Du,tc.Dr);
            b20 = psi(:, idx(2,0,tc.Lmax));
            b22 = psi(:, idx(2,2,tc.Lmax));

            b20_expected = tc.psi0(idx(2,0,tc.Lmax)) * exp(-6*tc.Dr*t');
            b22_expected = tc.psi0(idx(2,2,tc.Lmax)) * exp(-6*tc.Dr*t');

            tc.verifyEqual(b20, b20_expected, 'AbsTol',1e-4,'RelTol',1e-3);
            tc.verifyEqual(b22, b22_expected, 'AbsTol',1e-4,'RelTol',1e-3);
        end


        function testLongTimeLimitMatchesSteady(tc)
            % Long-time solution should match solve_steady
            t = linspace(0, 50/tc.Dr, 200);
            Du = [0,0,2.0;0,0,8.0;4.0,6.0,0];

            [~, psi] = solve_unsteady(t,tc.psi0,tc.Lmax,tc.beta,Du,tc.Dr);
            psi_unsteady = psi(end,:)';
            psi_steady = solve_steady(tc.Lmax, tc.beta, Du, ...
                'store', false);

            tc.verifyEqual(psi_unsteady, psi_steady, ...
                'AbsTol', 1e-12, 'RelTol', 1e-8);
        end


        function testConstantFunctionMatchesConstant(tc)
            % Constant time function should match constant for Du
            gxz = 2.0; gyz = 5.0;
            t = linspace(0, 2/tc.Dr, 200);
            Dus = [0,0,0;0,0,0;gxz,gyz,0];
            Dut = @(t) [0,0,0;0,0,0;gxz,gyz,0];

            [~, psis] = solve_unsteady(t,tc.psi0,tc.Lmax,tc.beta,Dus,tc.Dr);
            [~, psit] = solve_unsteady(t,tc.psi0,tc.Lmax,tc.beta,Dut,tc.Dr);

            tc.verifyEqual(psis, psit, 'AbsTol', 1e-12, 'RelTol', 1e-8);
        end


        function testXZvsYZOrderParamter(tc)
            % By reflection symmetry, gxz and gyz flows should give
            % same order parameter magnitude
            t = linspace(0, 5/tc.Dr, 50);
            N = idx(tc.Lmax, tc.Lmax, tc.Lmax);
            psi0_iso = zeros(N,1);
            psi0_iso(idx(0,0,tc.Lmax)) = 1/sqrt(4*pi);
            Du_xz = [0,0,0;0,0,0;20.0,0,0];
            Du_yz = [0,0,0;0,0,0;0,20.0,0];

            [~, psi_xz] = solve_unsteady(t, psi0_iso, tc.Lmax, tc.beta, ...
                Du_xz, tc.Dr);
            [~, psi_yz] = solve_unsteady(t, psi0_iso, tc.Lmax, tc.beta, ...
                Du_yz, tc.Dr);

            Q_xz = order_matrix(psi_xz, 'type', 'xz');
            Q_yz = order_matrix(psi_yz, 'type', 'xz');
            [Sy_xz, Sz_xz, Sx_xz] = order_parameters(Q_xz);
            [Sy_yz, Sz_yz, Sx_yz] = order_parameters(Q_yz);

            tc.verifyEqual(Sz_xz, Sz_yz, 'AbsTol', 1e-12, 'RelTol', 1e-8);
            tc.verifyEqual(Sx_xz, Sy_yz, 'AbsTol', 1e-12, 'RelTol', 1e-8);
            tc.verifyEqual(Sy_xz, Sx_yz, 'AbsTol', 1e-12, 'RelTol', 1e-8);
        end


        function testXZvsYZReflection(tc)
            % By reflection symmetry, gxz and gyz flows should give
            % same order parameter magnitude
            t = linspace(0, 5/tc.Dr, 50);
            N = idx(tc.Lmax, tc.Lmax, tc.Lmax);
            psi0_iso = zeros(N,1);
            psi0_iso(idx(0,0,tc.Lmax)) = 1/sqrt(4*pi);
            Du_xz = [0,0,0;0,0,0;20.0,0,0];
            Du_yz = [0,0,0;0,0,0;0,20.0,0];

            [~, psi_xz] = solve_unsteady(t, psi0_iso, tc.Lmax, tc.beta, ...
                Du_xz, tc.Dr);
            [~, psi_yz] = solve_unsteady(t, psi0_iso, tc.Lmax, tc.beta, ...
                Du_yz, tc.Dr);

            T = speye(N);
            for l = 0:2:tc.Lmax
                for m = 1:l
                    ii = idx(l,  m, tc.Lmax);
                    jj = idx(l, -m, tc.Lmax);
                    c = cos(m*pi/2); s = sin(m*pi/2);
                    % 2x2 reflection block for (+m, -m) pair
                    T(ii,ii) = c; T(ii,jj) =  s;
                    T(jj,ii) = s; T(jj,jj) = -c;
                end
            end
            psi_yzr = (T * psi_xz')';
            tc.verifyEqual(psi_yzr, psi_yz, 'AbsTol',1e-12,'RelTol',1e-8);
        end


        function testXZvsMixedRelfection(tc)
            % reflection along tan(theta)=(syz/sxy) and with
            % (sxz^2+syz^2)^0.5 shear-rate
            sxz = 100.0; syz = 20.0;
            smz = (sxz^2+syz^2)^0.5;
            theta = atan2(syz,sxz);
            t = linspace(0, 5/tc.Dr, 50);
            N = idx(tc.Lmax, tc.Lmax, tc.Lmax);
            psi0_iso = zeros(N,1);
            psi0_iso(idx(0,0,tc.Lmax)) = 1/sqrt(4*pi);

            Du_xz = [0,0,0;0,0,0;smz,0,0];
            [~, psi_xz] = solve_unsteady(t, psi0_iso, tc.Lmax, tc.beta, ...
                Du_xz, tc.Dr);

            Du_mz = [0,0,0;0,0,0;sxz,syz,0];
            [~, psi_mz] = solve_unsteady(t, psi0_iso, tc.Lmax, tc.beta, ...
                Du_mz, tc.Dr);

            T = speye(N);
            for l = 0:2:tc.Lmax
                for m = 1:l
                    ii = idx(l,  m, tc.Lmax);
                    jj = idx(l, -m, tc.Lmax);
                    c = cos(m*theta); s = sin(m*theta);
                    % 2x2 reflection block for (+m, -m) pair
                    T(ii,ii) = c; T(ii,jj) =  s;
                    T(jj,ii) = s; T(jj,jj) = -c;
                end
            end
            psi_mzr = (T * psi_xz')';
            tc.verifyEqual(psi_mzr, psi_mz, 'AbsTol',1e-4,'RelTol',1e-3);
        end

    end
end