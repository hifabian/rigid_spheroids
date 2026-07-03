classdef TestOrderParameters < matlab.unittest.TestCase
    methods (Test)
        function testIsotropicGivesZeroOrderParameters(tc)
            Q = [0, 0, 0, 0, 0, 0];
            [Sy, Sz, Sx, S, R] = order_parameters(Q);
            % No angle check because geometric multiplicity is 3
            tc.verifyEqual(Sy, 0, 'AbsTol', 1e-12);
            tc.verifyEqual(Sz, 0, 'AbsTol', 1e-12);
            tc.verifyEqual(Sx, 0, 'AbsTol', 1e-12);
            tc.verifyEqual(S,  0, 'AbsTol', 1e-12);
            tc.verifyEqual(R,  0, 'AbsTol', 1e-12);
        end

        function testIsotropicNonZeroTrace(tc)
            % Non-zero but isotropic trace should still give S=R=0
            % since the function subtracts trQ/3.
            a = 3.7;
            Q = [a, a, a, 0, 0, 0];
            % No angle check because geometric multiplicity is 3
            [Sy, Sz, Sx, S, R] = order_parameters(Q);
            tc.verifyEqual(Sy, 0, 'AbsTol', 1e-12);
            tc.verifyEqual(Sz, 0, 'AbsTol', 1e-12);
            tc.verifyEqual(Sx, 0, 'AbsTol', 1e-12);
            tc.verifyEqual(S,  0, 'AbsTol', 1e-12);
            tc.verifyEqual(R,  0, 'AbsTol', 1e-12);
        end

        function testSRAgainstEig_Diagonal(tc)
            % Purely diagonal traceless tensor: easy to verify by hand
            % and against eig.
            Q = [0.5, -0.3, -0.2, 0, 0, 0];
            [~,~,~, S, R] = order_parameters(Q);
            [S_ref, R_ref] = tc.SR_from_eig(Q);
            tc.verifyEqual(S, S_ref, 'AbsTol', 1e-12);
            tc.verifyEqual(R, R_ref, 'AbsTol', 1e-12);
        end

        function testSRAgainstEig_FullTensor(tc)
            % Full off-diagonal tensor
            Q = [0.3, -0.1, -0.2, 0.15, -0.05, 0.08];
            [~,~,~, S, R] = order_parameters(Q);
            [S_ref, R_ref] = tc.SR_from_eig(Q);
            tc.verifyEqual(S, S_ref, 'AbsTol', 1e-12);
            tc.verifyEqual(R, R_ref, 'AbsTol', 1e-12);
        end

        function testSRAgainstEig_NonTraceless(tc)
            % Non-traceless input: function subtracts trQ/3 internally.
            Q = [1.3, 0.7, 0.8, 0.15, -0.05, 0.08];
            [~,~,~, S, R] = order_parameters(Q);
            [S_ref, R_ref] = tc.SR_from_eig(Q);
            tc.verifyEqual(S, S_ref, 'AbsTol', 1e-12);
            tc.verifyEqual(R, R_ref, 'AbsTol', 1e-12);
        end

        function testSRAgainstEig_Uniaxial(tc)
            % Uniaxial tensor: R should be 0, S = l1 - l3.
            % Construct via known director n = [0,0,1], S0 = 0.6:
            %   M = S0*(n⊗n - I/3) = S0*diag([-1/3,-1/3,2/3])
            S0 = 0.6;
            Q = S0*[-1/3, -1/3, 2/3, 0, 0, 0];
            [~,~,~, S, R] = order_parameters(Q);
            tc.verifyEqual(S, S0,  'AbsTol', 1e-12);
            tc.verifyEqual(R, 0,   'AbsTol', 1e-12);
        end

        function testSRAgainstEig_Biaxial(tc)
            % Biaxial tensor: R > 0.
            Q = [0.4, 0.1, -0.5, 0, 0, 0];
            [~,~,~, S, R] = order_parameters(Q);
            [S_ref, R_ref] = tc.SR_from_eig(Q);
            tc.verifyEqual(S, S_ref, 'AbsTol', 1e-12);
            tc.verifyEqual(R, R_ref, 'AbsTol', 1e-12);
            tc.verifyGreaterThan(R, 0);
        end

        function testProjectedOrderParameters_Diagonal(tc)
            % For diagonal Q with known eigenvalues, Sy/Sz/Sx reduce to
            % simple differences.
            Q = [0.5, -0.3, -0.2, 0, 0, 0];
            [Sy, Sz, Sx] = order_parameters(Q);
            % Sy = |Q1-Q3|, Sz = |Q1-Q2|, Sx = |Q2-Q3|
            tc.verifyEqual(Sy, abs(Q(1)-Q(3)), 'AbsTol', 1e-12);
            tc.verifyEqual(Sz, abs(Q(1)-Q(2)), 'AbsTol', 1e-12);
            tc.verifyEqual(Sx, abs(Q(2)-Q(3)), 'AbsTol', 1e-12);
        end

        function testProjectedOrderParametersNonNegative(tc)
            Q = [0.3, -0.1, -0.2, 0.15, -0.05, 0.08];
            [Sy, Sz, Sx] = order_parameters(Q);
            tc.verifyGreaterThanOrEqual(Sy, 0);
            tc.verifyGreaterThanOrEqual(Sz, 0);
            tc.verifyGreaterThanOrEqual(Sx, 0);
        end

        function testAnglesAgainstEig_Diagonal_ZAligned(tc)
            % Director along z: ExtTheta=0, ExtChi undefined (any value ok)
            S0 = 0.6;
            Q = S0*[-1/3, -1/3, 2/3, 0, 0, 0];
            [~,~,~,~,~, ~, ExtTheta] = order_parameters(Q);
            tc.verifyEqual(ExtTheta, 0, 'AbsTol', 1e-12);
        end

        function testAnglesAgainstEig_XAligned(tc)
            % Director along x: ExtTheta=pi/2, ExtChi=0
            S0 = 0.6;
            Q = S0*[2/3, -1/3, -1/3, 0, 0, 0];
            [~,~,~,~,~, ExtChi, ExtTheta] = order_parameters(Q);
            tc.verifyEqual(ExtTheta, pi/2, 'AbsTol', 1e-12);
            tc.verifyEqual(ExtChi,   0,    'AbsTol', 1e-12);
        end

        function testAnglesAgainstEig_YAligned(tc)
            % Director along y: ExtTheta=pi/2, ExtChi=pi/2
            S0 = 0.6;
            Q = S0*[-1/3, 2/3, -1/3, 0, 0, 0];
            [~,~,~,~,~, ExtChi, ExtTheta] = order_parameters(Q);
            tc.verifyEqual(ExtTheta, pi/2, 'AbsTol', 1e-12);
            tc.verifyEqual(ExtChi,   pi/2, 'AbsTol', 1e-12);
        end

        function testAnglesAgainstEig_FullTensor(tc)
            % General tensor: compare angles to eig-derived eigenvector.
            Q = [0.3, -0.1, -0.2, 0.15, -0.05, 0.08];
            [~,~,~,~,~, ExtChi, ExtTheta] = order_parameters(Q);
            [chi_ref, theta_ref] = tc.angles_from_eig(Q);
            tc.verifyEqual(ExtChi,   chi_ref,   'AbsTol', 1e-10);
            tc.verifyEqual(ExtTheta, theta_ref, 'AbsTol', 1e-10);
        end

        function testAnglesAgainstEig_TiltedDirector(tc)
            % Construct a tensor with known tilted director and verify.
            theta0 = pi/5;
            chi0   = pi/7;
            n = [sin(theta0)*cos(chi0); sin(theta0)*sin(chi0); cos(theta0)];
            S0 = 0.7;
            M  = S0*(n*n' - eye(3)/3);
            Q  = [M(1,1), M(2,2), M(3,3), M(1,2), M(2,3), M(1,3)];
            [~,~,~,~,~, ExtChi, ExtTheta] = order_parameters(Q);
            tc.verifyEqual(ExtTheta, theta0, 'AbsTol', 1e-12);
            tc.verifyEqual(ExtChi,   chi0,   'AbsTol', 1e-12);
        end

        function testAnglesGaugeFix_NegativeZ(tc)
            % When true eigenvector has vz < 0, gauge-fixed result
            % should flip the sign so vz >= 0.
            theta0 = 3*pi/4;  % vz = cos(3pi/4) < 0, so should be flipped
            chi0   = pi/4;    % original chi
            n_orig = [sin(theta0)*cos(chi0); sin(theta0)*sin(chi0); cos(theta0)];
            % flip to get vz > 0 version
            n = -n_orig;
            theta_fix = acos(n(3));
            chi_fix   = atan2(n(2), n(1));
            S0 = 0.6;
            M  = S0*(n_orig*n_orig' - eye(3)/3);  % same tensor either way
            Q  = [M(1,1), M(2,2), M(3,3), M(1,2), M(2,3), M(1,3)];
            [~,~,~,~,~, ExtChi, ExtTheta] = order_parameters(Q);
            tc.verifyEqual(ExtTheta, theta_fix, 'AbsTol', 1e-12);
            tc.verifyEqual(ExtChi,   chi_fix,   'AbsTol', 1e-12);
        end

        function testBatchInputMatchesLoopOverRows(tc)
            % order_parameters should handle an (N x 6) input, giving
            % results identical to calling it row-by-row.
            rng(42);
            Qbatch = randn(8, 6) * 0.3;
            [Sy_b, Sz_b, Sx_b, S_b, R_b, Chi_b, Theta_b] = ...
                order_parameters(Qbatch);

            for i = 1:8
                [Sy_i, Sz_i, Sx_i, S_i, R_i, Chi_i, Theta_i] = ...
                    order_parameters(Qbatch(i,:));
                tc.verifyEqual(Sy_b(i),    Sy_i,    'AbsTol', 1e-12);
                tc.verifyEqual(Sz_b(i),    Sz_i,    'AbsTol', 1e-12);
                tc.verifyEqual(Sx_b(i),    Sx_i,    'AbsTol', 1e-12);
                tc.verifyEqual(S_b(i),     S_i,     'AbsTol', 1e-12);
                tc.verifyEqual(R_b(i),     R_i,     'AbsTol', 1e-12);
                tc.verifyEqual(Chi_b(i),   Chi_i,   'AbsTol', 1e-12);
                tc.verifyEqual(Theta_b(i), Theta_i, 'AbsTol', 1e-12);
            end
        end

        function testBatchSRAgainstEig(tc)
            % Batch S and R should each match eig() row by row.
            rng(7);
            Qbatch = randn(6, 6) * 0.3;
            [~,~,~, S_b, R_b] = order_parameters(Qbatch);

            for i = 1:6
                [S_ref, R_ref] = tc.SR_from_eig(Qbatch(i,:));
                tc.verifyEqual(S_b(i), S_ref, 'AbsTol', 1e-11);
                tc.verifyEqual(R_b(i), R_ref, 'AbsTol', 1e-11);
            end
        end

        function testSIsLargerThanR(tc)
            % By definition l1 >= l2 >= l3 => S = l1-l3 >= l2-l3 = R >= 0
            rng(13);
            Qbatch = randn(20, 6) * 0.3;
            [~,~,~, S, R] = order_parameters(Qbatch);
            tc.verifyGreaterThanOrEqual(S, R - 1e-12*ones(size(R)));
            tc.verifyGreaterThanOrEqual(R, -1e-12*ones(size(R)));
        end

        function testSBoundsFromProjected(tc)
            % For any symmetric tensor, the full order parameter S must
            % be >= each projected order parameter / sqrt(3):
            %   S >= Sx/sqrt(3),  S >= Sy/sqrt(3),  S >= Sz/sqrt(3)
            % (follows from Cauchy interlacing / projection inequality)
            rng(99);
            Qbatch = randn(20, 6) * 0.3;
            [Sy, Sz, Sx, S] = order_parameters(Qbatch);
            tol = 1e-12;
            tc.verifyGreaterThanOrEqual(S, Sx/sqrt(3) - tol);
            tc.verifyGreaterThanOrEqual(S, Sy/sqrt(3) - tol);
            tc.verifyGreaterThanOrEqual(S, Sz/sqrt(3) - tol);
        end

    end

    
    methods (Static, Access = private)

        function [S, R] = SR_from_eig(Q)
            % Compute S = l1-l3, R = l2-l3 from the 3x3 matrix via eig().
            % Subtracts trace/3 to match what order_parameters does.
            M = TestOrderParameters.Q_to_matrix(Q);
            M = M - trace(M)/3 * eye(3);
            lam = sort(real(eig(M)), 'descend');  % l1 >= l2 >= l3
            S = lam(1) - lam(3);
            R = lam(2) - lam(3);
        end

        function [chi, theta] = angles_from_eig(Q)
            % Compute ExtChi and ExtTheta from eig(), applying the same
            % gauge convention as order_parameters: vz >= 0.
            M = TestOrderParameters.Q_to_matrix(Q);
            M = M - trace(M)/3 * eye(3);
            [V, D] = eig(M);
            [~, k] = max(real(diag(D)));
            v = real(V(:,k));
            v = v / norm(v);
            if v(3) < 0
                v = -v;
            end
            theta = acos(min(max(v(3), -1), 1));
            chi   = atan2(v(2), v(1));
        end

        function M = Q_to_matrix(Q)
            % Unpack Q = [Qxx, Qyy, Qzz, Qxy, Qyz, Qzx] into 3x3 matrix.
            M = [Q(1), Q(4), Q(6); ...
                 Q(4), Q(2), Q(5); ...
                 Q(6), Q(5), Q(3)];
        end

    end
end