classdef TestFpUnsteady < matlab.unittest.TestCase
    properties
        Lmax = 24;
        beta = 0.5;
        l0 = 1e-5;
        Dr = 1;
    end

    methods (Test)

        function testMonodisperseMatchesSolveUnsteady(tc)
            sxz0 = 0; syz0 = 0;  % isotropic initial condition
            Du0 = [0,0,0;0,0,0;sxz0,syz0,0];
            init = fp_init(Du0, tc.l0, 1, tc.Dr, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            sxz = 2.0; syz = 1.0;
            Du = [0,0,0;0,0,0;sxz,syz,0];
            T = 5/init.Dr;
            t = linspace(0, T, 30);

            result = fp_unsteady(init, t, Du, ...
                'verbose', false, 'store', false);
            [~, psi_direct] = solve_unsteady(t, init.psi0{1}, init.Dr, ...
                tc.beta, Du);
            Q_direct = order_matrix(psi_direct);

            tc.verifyEqual(result.Q, Q_direct, ...
                'AbsTol', 1e-8, 'RelTol', 1e-6);
        end

        function testDegeneratePolydisperseMatchesMonodisperse(tc)
            % A "polydisperse" grid of repeated identical lengths should
            % give the same transient order parameters as the
            % monodisperse case.
            sxz0 = 0; syz0 = 0;
            Du0 = [0,0,0;0,0,0;sxz0,syz0,0];

            init_mono = fp_init(Du0, tc.l0, 1, tc.Dr, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            q_poly = [tc.l0-1e-12,tc.l0,tc.l0+1e-12]; % Close enough
            w_poly = [1;1;1];
            w_poly = w_poly / sum(w_poly);
            Dr_poly = ones(3)*tc.Dr;
            init_poly = fp_init(Du0, q_poly, w_poly, Dr_poly, ...
                tc.beta, 'Lmax', tc.Lmax, 'store', false);

            sxz = 1.5; syz = 0.5;
            Du = [0,0,0;0,0,0;sxz,syz,0];
            T = 5/init_mono.Dr;
            t = linspace(0, T, 25);

            result_mono = fp_unsteady(init_mono, t, Du, ...
                'store', false);
            result_poly = fp_unsteady(init_poly, t, Du, ...
                'store', false);

            tc.verifyEqual(result_poly.Q, result_mono.Q, ...
                'AbsTol', 1e-8, 'RelTol', 1e-6);
        end

        function testZeroFlowDecaysToIsotropic(tc)
            % Start anisotropic (via nonzero sxz0 in fp_init), then
            % relax under zero flow -> order parameters decay to zero.
            Du0 = [0,0,10.0;0,0,2.0;5.0,8.0,0];
            init = fp_init(Du0, tc.l0, 1, tc.Dr, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);
            T = 20/init.Dr;
            t = linspace(0, T, 60);

            Du = [0,0,0;0,0,0;0,0,0];
            result = fp_unsteady(init, t, Du, 'store', false);

            tc.verifyEqual(result.Q(end,:), zeros(size(result.Q(end,:))), ...
                'AbsTol', 1e-8);
        end

        function testLongTimeLimitMatchesFpSteady(tc)
            % At long times under constant flow, the transient solution
            % should relax to the steady-state result from fp_steady
            % for the same (dimensional) shear rates and rod length.
            Du0 = [0,0,0;0,0,0;0,0,0];
            init = fp_init(Du0, tc.l0, 1, tc.Dr, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            Du = [0,0,2.0;0,0,8.0;4.0,6.0,0];
            T = 60/init.Dr;
            t = linspace(0, T, 200);

            result_unsteady = fp_unsteady(init, t, Du, ...
                'store', false);
            result_steady = fp_steady(Du, tc.l0, 1, tc.Dr, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            tc.verifyEqual(result_unsteady.Q(end,:), result_steady.Q, ...
                'AbsTol', 1e-10, 'RelTol', 1e-8);
        end

        function testTrueLengthDependenceChangesResult(tc)
            % A genuinely polydisperse population (different lengths,
            % different Dr) should not generally reduce to the
            % monodisperse transient result.
            sxz0 = 0; syz0 = 0;
            Du0 = [0,0,sxz0;0,0,syz0;0,0,0];
            init_mono = fp_init(Du0, tc.l0, 1, tc.Dr, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            q_poly = tc.l0*[0.5, 1.5];
            w_poly = [1.0; 1.0];
            w_poly = w_poly / sum(w_poly);
            Dr_poly = [0.5, 1.5];
            init_poly = fp_init(Du0, q_poly, w_poly, Dr_poly, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            sxz = 2.0; syz = 0.0;
            Du = [0,0,sxz;0,0,syz;0,0,0];
            T = 5/init_mono.Dr;
            t = linspace(0, T, 20);

            result_mono = fp_unsteady(init_mono, t, Du, ...
                'store', false);
            result_poly = fp_unsteady(init_poly, t, Du, ...
                'store', false);

            tc.verifyNotEqual(result_poly.Q, result_mono.Q);
        end


        function testConstantFunctionMatchesConstant(tc)
            % Constant time function should match constant for Du
            Du0 = [0,0,0;0,0,0;0,0,0];
            init = fp_init(Du0, tc.l0, 1, tc.Dr, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            gxz = 2.0; gyz = 5.0;
            t = linspace(0, 2/tc.Dr, 200);
            Dus = [0,0,0;0,0,0;gxz,gyz,0];
            Dut = @(t) [0,0,0;0,0,0;gxz,gyz,0];

            psis = fp_unsteady(init, t, Dus, 'store', false);
            psit = fp_unsteady(init, t, Dut, 'store', false);

            tc.verifyEqual(psis.Q, psit.Q, 'AbsTol', 1e-12, 'RelTol', 1e-8);
        end

    end
end