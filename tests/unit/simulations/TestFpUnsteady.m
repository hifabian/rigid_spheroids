classdef TestFpUnsteady < matlab.unittest.TestCase
    properties
        Lmax = 24;
        beta = 0.5;
        l0 = 1e-9;
    end

    methods (Test)

        function testMonodisperseMatchesSolveUnsteady(tc)
            sxz0 = 0; syz0 = 0;  % isotropic initial condition
            init = fp_init(sxz0, syz0, tc.l0, 1, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            sxz = 2.0; syz = 1.0;
            T = 5/init.Dr;
            t = linspace(0, T, 30);

            result = fp_unsteady(init, t, sxz, syz, ...
                'verbose', false, 'store', false);

            [~, psi_direct] = solve_unsteady(t, init.psi0{1}, tc.Lmax, ...
                tc.beta, sxz, syz, init.Dr);
            Q_direct = order_matrix(psi_direct, 'type', 'xz');
            [Sy, Sz, Sx, ExtChi, ExtTheta] = order_parameters(Q_direct);

            tc.verifyEqual(result.Sy, Sy', 'AbsTol', 1e-8, 'RelTol', 1e-6);
            tc.verifyEqual(result.Sz, Sz', 'AbsTol', 1e-8, 'RelTol', 1e-6);
            tc.verifyEqual(result.Sx, Sx', 'AbsTol', 1e-8, 'RelTol', 1e-6);
            tc.verifyEqual(result.ExtChi, ExtChi', 'AbsTol', 1e-8, 'RelTol', 1e-6);
            tc.verifyEqual(result.ExtTheta, ExtTheta', 'AbsTol', 1e-8, 'RelTol', 1e-6);
        end

        function testDegeneratePolydisperseMatchesMonodisperse(tc)
            % A "polydisperse" grid of repeated identical lengths should
            % give the same transient order parameters as the
            % monodisperse case.
            sxz0 = 0; syz0 = 0;
            init_mono = fp_init(sxz0, syz0, tc.l0, 1, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            lv_poly = [tc.l0-1e-12,tc.l0,tc.l0+1e-12]; % Close enough
            fv_poly = [1,1,1];
            fv_poly = fv_poly / trapz(lv_poly, fv_poly);
            init_poly = fp_init(sxz0, syz0, lv_poly, fv_poly, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            sxz = 1.5; syz = 0.5;
            T = 5/init_mono.Dr;
            t = linspace(0, T, 25);

            result_mono = fp_unsteady(init_mono, t, sxz, syz, ...
                'store', false);
            result_poly = fp_unsteady(init_poly, t, sxz, syz, ...
                'store', false);

            tc.verifyEqual(result_poly.Sy, result_mono.Sy, ...
                'AbsTol', 1e-8, 'RelTol', 1e-6);
            tc.verifyEqual(result_poly.Sz, result_mono.Sz, ...
                'AbsTol', 1e-8, 'RelTol', 1e-6);
            tc.verifyEqual(result_poly.Sx, result_mono.Sx, ...
                'AbsTol', 1e-8, 'RelTol', 1e-6);
            tc.verifyEqual(result_poly.ExtChi, result_mono.ExtChi, ...
                'AbsTol', 1e-8, 'RelTol', 1e-6);
            tc.verifyEqual(result_poly.ExtTheta, result_mono.ExtTheta, ...
                'AbsTol', 1e-8, 'RelTol', 1e-6);
        end

        function testZeroFlowDecaysToIsotropic(tc)
            % Start anisotropic (via nonzero sxz0 in fp_init), then
            % relax under zero flow -> order parameters decay to zero.
            init = fp_init(5.0, 0.0, tc.l0, 1, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);
            T = 20/init.Dr;
            t = linspace(0, T, 60);

            result = fp_unsteady(init, t, 0, 0, 'store', false);

            tc.verifyEqual(result.Sy(end), 0, 'AbsTol', 1e-6);
            tc.verifyEqual(result.Sz(end), 0, 'AbsTol', 1e-6);
            tc.verifyEqual(result.Sx(end), 0, 'AbsTol', 1e-6);
        end

        function testLongTimeLimitMatchesFpSteady(tc)
            % At long times under constant flow, the transient solution
            % should relax to the steady-state result from fp_steady
            % for the same (dimensional) shear rates and rod length.
            sxz0 = 0; syz0 = 0;
            init = fp_init(sxz0, syz0, tc.l0, 1, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            sxz = 3.0; syz = 1.0;
            T = 60/init.Dr;
            t = linspace(0, T, 200);

            result_unsteady = fp_unsteady(init, t, sxz, syz, ...
                'store', false);
            result_steady = fp_steady(sxz, syz, tc.l0, 1, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            tc.verifyEqual(result_unsteady.Sy(end), result_steady.Sy, ...
                'AbsTol', 1e-6, 'RelTol', 1e-4);
            tc.verifyEqual(result_unsteady.Sz(end), result_steady.Sz, ...
                'AbsTol', 1e-6, 'RelTol', 1e-4);
            tc.verifyEqual(result_unsteady.Sx(end), result_steady.Sx, ...
                'AbsTol', 1e-6, 'RelTol', 1e-4);
        end

        function testTrueLengthDependenceChangesResult(tc)
            % A genuinely polydisperse population (different lengths,
            % different Dr) should not generally reduce to the
            % monodisperse transient result.
            sxz0 = 0; syz0 = 0;
            init_mono = fp_init(sxz0, syz0, tc.l0, 1, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            lv_poly = tc.l0*[0.5, 1.5];
            fv_poly = [1.0, 1.0];
            init_poly = fp_init(sxz0, syz0, lv_poly, fv_poly, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            sxz = 2.0; syz = 0.0;
            T = 5/init_mono.Dr;
            t = linspace(0, T, 20);

            result_mono = fp_unsteady(init_mono, t, sxz, syz, ...
                'store', false);
            result_poly = fp_unsteady(init_poly, t, sxz, syz, ...
                'store', false);

            tc.verifyNotEqual(result_poly.Sy, result_mono.Sy);
        end

    end
end