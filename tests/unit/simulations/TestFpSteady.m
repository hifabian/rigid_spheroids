classdef TestFpSteady < matlab.unittest.TestCase
    properties
        Lmax = 24;
        beta = 0.5;
        l0 = 1e-5;   % rod length (m), arbitrary but fixed
        Dr = 1.0;
    end

    methods (Test)

        function testMonodisperseMatchesSolveSteady(tc)
            % fp_steady with scalar lv should exactly reproduce
            % solve_steady once shear rates are converted to Peclet
            % numbers using the resulting Dr.
            sxz = 3.0; syz = 1.5;
            Du = [0,0,0;0,0,0;sxz,syz,0];

            result = fp_steady(Du, tc.l0, 1, tc.Dr, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);
            psi_direct = solve_steady(tc.Lmax, tc.beta, Du/tc.Dr, ...
                'store', false);
            Q_direct = order_matrix(psi_direct);

            tc.verifyEqual(result.Q, Q_direct, ...
                'RelTol', 1e-8, 'AbsTol', 1e-10);
        end


        function testDegeneratePolydisperseMatchesMonodisperse(tc)
            sxz = 10.0; syz = 2.0;
            Du = [0,0,0;0,0,0;sxz,syz,0];
            lv_poly = [tc.l0-1e-10,tc.l0,tc.l0+1e-10]; % Close enough
            fv_poly = [1,1,1];
            fv_poly = fv_poly / trapz(lv_poly, fv_poly);
            Dr_poly = ones(size(fv_poly))*tc.Dr;

            result_mono = fp_steady(Du, tc.l0, 1, tc.Dr, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);
            result_poly = fp_steady(Du, lv_poly, fv_poly, Dr_poly, ...
                tc.beta, 'Lmax', tc.Lmax, 'store', false);

            tc.verifyEqual(result_poly.Q, result_mono.Q, ...
                'AbsTol', 1e-10, 'RelTol', 1e-8);
        end


        function testZeroFlowIsIsotropic(tc)
            % No flow -> isotropic distribution -> Q == 0
            Du = [0,0,0;0,0,0;0,0,0];

            result = fp_steady(Du, tc.l0, 1, tc.Dr, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            tc.verifyEqual(result.Q, zeros(size(result.Q)), ...
                'AbsTol', 1e-10);
        end


        function testMultipleShearRatesMatchSolveSteadyEachXZ(tc)
            % fp_steady can be called with a vector of sxz (or syz) at
            % once; each entry should match a separate solve_steady call.
            Du = repmat([0,0,0;0,0,0;0,0,0], 1, 1, 3);
            Du(3,1,:) = [1.0, 2.0, 4.0];

            result = fp_steady(Du, tc.l0, 1, tc.Dr, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            for i = 1:length(Du(3,1))
                psi = solve_steady(tc.Lmax, tc.beta, Du(:,:,i)/tc.Dr, ...
                    'store', false);
                Q = order_matrix(psi);

                tc.verifyEqual(result.Q(i,:), Q, ...
                    'AbsTol', 1e-10, 'RelTol', 1e-8);
            end
        end


        function testMultipleShearRatesMatchSolveSteadyEachYZ(tc)
            % fp_steady can be called with a vector of sxz (or syz) at
            % once; each entry should match a separate solve_steady call.
            Du = repmat([0,0,0;0,0,0;0,0,0], 1, 1, 3);
            Du(3,2,:) = [1.0, 2.0, 4.0];

            result = fp_steady(Du, tc.l0, 1, tc.Dr, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            for i = 1:length(Du(3,2))
                psi = solve_steady(tc.Lmax, tc.beta, Du(:,:,i)/tc.Dr, ...
                    'store', false);
                Q = order_matrix(psi);

                tc.verifyEqual(result.Q(i,:), Q, ...
                    'AbsTol', 1e-10, 'RelTol', 1e-8);
            end
        end


        function testTrueLengthDependenceChangesResult(tc)
            % Sanity check: a polydisperse population with genuinely
            % different lengths (and thus different Dr) should NOT
            % generally reduce to the single-length monodisperse result.
            sxz = 5.0; syz = 0.0;
            Du = [0,0,0;0,0,0;sxz,syz,0];

            result_mono = fp_steady(Du, tc.l0, 1, tc.Dr, tc.beta, ...
                'Lmax', tc.Lmax, 'store', false);

            lv_poly = tc.l0*[0.5, 1.0, 1.5];
            fv_poly = [0.5, 1.0, 0.5];  % non-constant weighting
            Dr_poly = ones(size(fv_poly))*tc.Dr;
            result_poly = fp_steady(Du, lv_poly, fv_poly, Dr_poly, ...
                tc.beta, 'Lmax', tc.Lmax, 'store', false);

            tc.verifyNotEqual(result_poly.Q, result_mono.Q);
        end

    end
end