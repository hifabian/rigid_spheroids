classdef TestFpInit < matlab.unittest.TestCase
    properties
        Lmax = 12;
        beta = 0.5;
        l0 = 1e-5;
        Dr = 1;
    end

    methods (Test)

        function testMonodisperseMatchesSolveSteady(tc)
            sxz0 = 3.0; syz0 = 1.5;
            Du0 = [0,0,0;0,0,0;sxz0,syz0,0];

            result = fp_init(Du0, tc.l0, 1, tc.Dr, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);
            psi_direct = solve_steady(tc.Lmax, tc.beta, Du0/result.Dr, ...
                'store', false);

            tc.verifyEqual(result.psi0{1}, psi_direct, 'AbsTol', 1e-10);
        end

        function testDegeneratePolydisperseMatchesMonodisperse(tc)
            % Repeating the same length for every entry of a
            % "polydisperse" grid should give identical psi0{j} for
            % every j, equal to the monodisperse psi0.
            sxz0 = 2.0; syz0 = 0.5;
            Du0 = [0,0,0;0,0,0;sxz0,syz0,0];

            result_mono = fp_init(Du0, tc.l0, 1, tc.Dr, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            lv_poly = tc.l0*ones(1,4);
            fv_poly = ones(1,4);
            Dr_poly = ones(1,4)*tc.Dr;
            result_poly = fp_init(Du0, lv_poly, fv_poly, Dr_poly, ...
                tc.beta, 'Lmax', tc.Lmax, 'store', false);

            for j = 1:length(lv_poly)
                tc.verifyEqual(result_poly.psi0{j}, result_mono.psi0{1}, ...
                    'AbsTol', 1e-10);
            end
        end

        function testZeroFlowIsIsotropicCoefficients(tc)
            % No flow -> isotropic distribution: only the (0,0)
            % coefficient is nonzero, equal to 1/sqrt(4*pi).
            Du = [0,0,0;0,0,0;0,0,0];

            result = fp_init(Du, tc.l0, 1, tc.Dr, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);
            psi0 = result.psi0{1};

            psi_iso = zeros(size(psi0));
            psi_iso(idx(0,0,tc.Lmax)) = 1/sqrt(4*pi);

            tc.verifyEqual(psi0, psi_iso, 'AbsTol', 1e-10);
        end

        function testNormalization(tc)
            % b_00 = 1/sqrt(4*pi) regardless of flow (probability
            % conservation)
            Du = [0,0,0;0,0,0;4.0,2.0,0];

            result = fp_init(Du, tc.l0, 1, tc.Dr, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            psi0 = result.psi0{1};
            tc.verifyEqual(psi0(idx(0,0,tc.Lmax)), 1/sqrt(4*pi), ...
                'AbsTol', 1e-10);
        end

        function testDifferentLengthsGiveDifferentDr(tc)
            % Sanity check that Dr actually depends on length, so a
            % genuinely polydisperse population is not degenerate.
            lv_poly = tc.l0*[0.5, 1.0, 2.0];
            Du = [0,0,0;0,0,0;1.0,0,0];
            Dr_poly = [1, 2, 3];

            result = fp_init(Du, lv_poly, ones(1,3), Dr_poly, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            tc.verifyNotEqual(result.Dr(1), result.Dr(2));
            tc.verifyNotEqual(result.Dr(2), result.Dr(3));
        end

        function testTrueLengthDependenceChangesResult(tc)
            % Genuinely different lengths should not generally reduce to
            % the same psi0 as a single fixed length.
            sxz0 = 5.0; syz0 = 2.0;
            Du0 = [0,0,0;0,0,0;sxz0,syz0,0];

            result_mono = fp_init(Du0, tc.l0, 1, tc.Dr, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            lv_poly = tc.l0*[0.5, 1.5];
            fv_poly = [1.0, 1.0];
            Dr_poly = [0.5, 1.5];
            result_poly = fp_init(Du0, lv_poly, fv_poly, Dr_poly, ...
                tc.beta, 'Lmax', tc.Lmax, 'store', false);

            tc.verifyNotEqual(result_poly.psi0{1}, result_mono.psi0{1});
        end

    end
end