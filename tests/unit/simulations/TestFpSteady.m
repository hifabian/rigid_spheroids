classdef TestFpSteady < matlab.unittest.TestCase
    properties
        Lmax = 24;
        beta = 0.5;
        l0 = 1e-5;   % rod length (m), arbitrary but fixed
    end

    methods (Test)

        function testMonodisperseMatchesSolveSteady(tc)
            % fp_steady with scalar lv should exactly reproduce
            % solve_steady once shear rates are converted to Peclet
            % numbers using the resulting Dr.
            sxz = 3.0; syz = 1.5;

            result = fp_steady(sxz, syz, tc.l0, 1, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            Dr = result.Dr;
            psi_direct = solve_steady(tc.Lmax, tc.beta, sxz/Dr, syz/Dr, ...
                'store', false);
            Q_direct = order_matrix(psi_direct, 'type', 'xz');
            [Sy, Sz, Sx, ExtChi, ExtTheta] = order_parameters(Q_direct);

            tc.verifyEqual(result.Sy, Sy, 'RelTol', 1e-8, 'AbsTol', 1e-10);
            tc.verifyEqual(result.Sz, Sz, 'RelTol', 1e-8, 'AbsTol', 1e-10);
            tc.verifyEqual(result.Sx, Sx, 'RelTol', 1e-8, 'AbsTol', 1e-10);
            tc.verifyEqual(result.ExtChi, ExtChi, 'RelTol', 1e-8, ...
                'AbsTol', 1e-10);
            tc.verifyEqual(result.ExtTheta, ExtTheta, 'RelTol', 1e-8, ...
                'AbsTol', 1e-10);
        end

        function testDegeneratePolydisperseMatchesMonodisperse(tc)
            sxz = 10.0; syz = 2.0;

            result_mono = fp_steady(sxz, syz, tc.l0, 1, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            lv_poly = [tc.l0-1e-10,tc.l0,tc.l0+1e-10]; % Close enough
            fv_poly = [1,1,1];
            fv_poly = fv_poly / trapz(lv_poly, fv_poly);
            result_poly = fp_steady(sxz, syz, lv_poly, fv_poly, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            tc.verifyEqual(result_poly.Sy, result_mono.Sy, ...
                'AbsTol', 1e-10, 'RelTol', 1e-8);
            tc.verifyEqual(result_poly.Sz, result_mono.Sz, ...
                'AbsTol', 1e-10, 'RelTol', 1e-8);
            tc.verifyEqual(result_poly.Sx, result_mono.Sx, ...
                'AbsTol', 1e-10, 'RelTol', 1e-8);
            tc.verifyEqual(result_poly.ExtChi, result_mono.ExtChi, ...
                'AbsTol', 1e-10, 'RelTol', 1e-8);
            tc.verifyEqual(result_poly.ExtTheta, result_mono.ExtTheta, ...
                'AbsTol', 1e-10, 'RelTol', 1e-8);
        end

        function testZeroFlowIsIsotropic(tc)
            % No flow -> isotropic distribution -> Sy=Sz=Sx=0
            % (same physical check as TestSolveSteady.testZeroFlowIsIsotropic)
            result = fp_steady(0, 0, tc.l0, 1, tc.beta, 'Lmax', tc.Lmax, ...
                'store', false);
            tc.verifyEqual(result.Sy, 0, 'AbsTol', 1e-10);
            tc.verifyEqual(result.Sz, 0, 'AbsTol', 1e-10);
            tc.verifyEqual(result.Sx, 0, 'AbsTol', 1e-10);
        end

        function testMultipleShearRatesMatchSolveSteadyEachXZ(tc)
            % fp_steady can be called with a vector of sxz (or syz) at
            % once; each entry should match a separate solve_steady call.
            sxz_vec = [1.0, 2.0, 4.0];
            syz0 = 0.0;

            result = fp_steady(sxz_vec, syz0, tc.l0, 1, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);
            Dr = result.Dr;

            for i = 1:length(sxz_vec)
                psi = solve_steady(tc.Lmax, tc.beta, ...
                    sxz_vec(i)/Dr, syz0/Dr, 'store', false);
                Q = order_matrix(psi, 'type', 'xz');
                [Sy, Sz, Sx] = order_parameters(Q);
                tc.verifyEqual(result.Sy(i), Sy, 'AbsTol', 1e-10);
                tc.verifyEqual(result.Sz(i), Sz, 'AbsTol', 1e-10);
                tc.verifyEqual(result.Sx(i), Sx, 'AbsTol', 1e-10);
            end
        end

        function testMultipleShearRatesMatchSolveSteadyEachYZ(tc)
            % fp_steady can be called with a vector of sxz (or syz) at
            % once; each entry should match a separate solve_steady call.
            syz_vec = [1.0, 2.0, 4.0];
            sxz0 = 0.0;

            result = fp_steady(syz_vec, sxz0, tc.l0, 1, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);
            Dr = result.Dr;

            for i = 1:length(syz_vec)
                psi = solve_steady(tc.Lmax, tc.beta, ...
                    syz_vec(i)/Dr, sxz0/Dr, 'store', false);
                Q = order_matrix(psi, 'type', 'xz');
                [Sy, Sz, Sx] = order_parameters(Q);
                tc.verifyEqual(result.Sy(i), Sy, 'AbsTol', 1e-10);
                tc.verifyEqual(result.Sz(i), Sz, 'AbsTol', 1e-10);
                tc.verifyEqual(result.Sx(i), Sx, 'AbsTol', 1e-10);
            end
        end

        function testTrueLengthDependenceChangesResult(tc)
            % Sanity check: a polydisperse population with genuinely
            % different lengths (and thus different Dr) should NOT
            % generally reduce to the single-length monodisperse result.
            sxz = 5.0; syz = 0.0;

            result_mono = fp_steady(sxz, syz, tc.l0, 1, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            lv_poly = tc.l0*[0.5, 1.0, 1.5];
            fv_poly = [0.5, 1.0, 0.5];  % non-constant weighting
            result_poly = fp_steady(sxz, syz, lv_poly, fv_poly, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            tc.verifyNotEqual(result_poly.Sy, result_mono.Sy);
        end

    end
end