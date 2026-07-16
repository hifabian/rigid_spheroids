classdef TestBuildQuadrature < matlab.unittest.TestCase
    properties
        % Common test parameters
        kGamma = 3.5;
        thetaGamma = 75.0;
        muLogN = 5.4;
        sigmaLogN = 0.6;
        nNodes = 8;      % default node count for gamma/lognormal exactness checks
        nBounded = 24;   % node count for bounded tests (needs to resolve a wide domain well)

        % Quadrature outputs, built once per test method
        qGamma; wGamma;
        qLogN; wLogN;
        qBounded; wBounded;
        qTrapz; wTrapz;

        % trapz test fixture (kept as properties so multiple tests can reuse it)
        lvTrapz; fvTrapz;
    end

    methods (TestMethodSetup)
        function buildQuadratures(tc)
            [tc.qGamma, tc.wGamma] = build_quadrature('gamma', tc.kGamma, tc.thetaGamma, tc.nNodes);
            [tc.qLogN, tc.wLogN]   = build_quadrature('lognormal', tc.muLogN, tc.sigmaLogN, tc.nNodes);

            pdfLogN = @(x) lognPdf(x, tc.muLogN, tc.sigmaLogN);
            [tc.qBounded, tc.wBounded] = build_quadrature('bounded', pdfLogN, 1e-6, 1500, tc.nBounded);

            % A deliberately non-uniform grid with an arbitrary (not
            % perfectly normalized) density sampled on it
            tc.lvTrapz = [10; 25; 60; 80; 150; 220; 400]';
            tc.fvTrapz = lognPdf(tc.lvTrapz, tc.muLogN, tc.sigmaLogN);
            [tc.qTrapz, tc.wTrapz] = build_quadrature('trapz', tc.lvTrapz, tc.fvTrapz);
        end
    end

    methods (Test)

        % ---------------------------------------------------------- gamma
        function testGammaWeightsSumToOne(tc)
            tc.verifyEqual(sum(tc.wGamma), 1, 'AbsTol', 1e-10);
        end

        function testGammaWeightsPositive(tc)
            tc.verifyGreaterThanOrEqual(tc.wGamma, 0);
        end

        function testGammaNodesPositive(tc)
            tc.verifyGreaterThanOrEqual(tc.qGamma, 0);
        end

        function testGammaMomentsExactUpToDegree(tc)
            % Gauss-Laguerre quadrature must be EXACT (to machine precision)
            % for moments up to degree 2*n - 1
            maxDegree = 2*tc.nNodes - 1;
            for m = 0:maxDegree
                quadMoment  = sum(tc.wGamma .* tc.qGamma.^m);
                exactMoment = tc.thetaGamma^m * gamma(tc.kGamma + m) / gamma(tc.kGamma);
                tc.verifyEqual(quadMoment, exactMoment, 'RelTol', 1e-8, ...
                    sprintf('Gamma quadrature moment mismatch at degree m=%d', m));
            end
        end

        function testGammaMeanMatchesAnalytic(tc)
            meanQuad = sum(tc.wGamma .* tc.qGamma);
            tc.verifyEqual(meanQuad, tc.kGamma*tc.thetaGamma, 'RelTol', 1e-8);
        end

        % ------------------------------------------------------ lognormal
        function testLognormalWeightsSumToOne(tc)
            tc.verifyEqual(sum(tc.wLogN), 1, 'AbsTol', 1e-10);
        end

        function testLognormalWeightsPositive(tc)
            tc.verifyGreaterThanOrEqual(tc.wLogN, 0);
        end

        function testLognormalNodesPositive(tc)
            tc.verifyGreaterThanOrEqual(tc.qLogN, 0);
        end

        function testLognormalLowOrderMomentsApproxAnalytic(tc)
            % Not exact (unlike gamma/Laguerre) -- Gauss-Hermite via the
            % log-transform is only spectrally convergent for L^m, since
            % L^m is exponential (not polynomial) in the Hermite variable.
            % Low-order moments should still match closely at n=8.
            for m = 0:3
                quadMoment  = sum(tc.wLogN .* tc.qLogN.^m);
                exactMoment = exp(m*tc.muLogN + 0.5*m^2*tc.sigmaLogN^2);
                tc.verifyEqual(quadMoment, exactMoment, 'RelTol', 1e-3, ...
                    sprintf('Lognormal quadrature moment mismatch at degree m=%d', m));
            end
        end

        function testLognormalMomentsConvergeWithMoreNodes(tc)
            % Error should shrink as n grows (spectral convergence),
            % even though it is never exact for finite n.
            m = 4;
            exactMoment = exp(m*tc.muLogN + 0.5*m^2*tc.sigmaLogN^2);

            [qSmall, wSmall] = build_quadrature('lognormal', tc.muLogN, tc.sigmaLogN, 6);
            [qLarge, wLarge] = build_quadrature('lognormal', tc.muLogN, tc.sigmaLogN, 16);

            errSmall = abs(sum(wSmall .* qSmall.^m) - exactMoment);
            errLarge = abs(sum(wLarge .* qLarge.^m) - exactMoment);

            tc.verifyLessThan(errLarge, errSmall);
        end

        % --------------------------------------------------------- bounded
        function testBoundedNodesWithinBounds(tc)
            tc.verifyGreaterThanOrEqual(tc.qBounded, 1e-6);
            tc.verifyLessThanOrEqual(tc.qBounded, 1500);
        end

        function testBoundedWeightsPositive(tc)
            tc.verifyGreaterThanOrEqual(tc.wBounded, 0);
        end

        function testBoundedCapturedMassMatchesAnalyticCDF(tc)
            trueMass = lognCdf(1500, tc.muLogN, tc.sigmaLogN) - lognCdf(1e-6, tc.muLogN, tc.sigmaLogN);
            tc.verifyEqual(sum(tc.wBounded), trueMass, 'AbsTol', 1e-3);
        end

        function testBoundedNoWarningWhenMassFullyCaptured(tc)
            pdfLogN = @(x) lognPdf(x, tc.muLogN, tc.sigmaLogN);
            tc.verifyWarningFree(...
                @() build_quadrature('bounded', pdfLogN, 1e-6, 1500, tc.nBounded));
        end

        function testBoundedWarnsWhenMassIsTruncated(tc)
            % lmax = 500 only captures ~90% of the lognormal mass here
            pdfLogN = @(x) lognPdf(x, tc.muLogN, tc.sigmaLogN);
            tc.verifyWarning(...
                @() build_quadrature('bounded', pdfLogN, 1e-6, 500, tc.nBounded), ...
                'build_quadrature:mass_loss');
        end

        function testBoundedIntegratesKnownIntegralCorrectly(tc)
            % sanity check against a simple, independently-known integral:
            % uniform pdf on [0, L] has E[x] = L/2
            Lflat = 100;
            pdfFlat = @(x) ones(size(x))/Lflat;
            [qFlat, wFlat] = build_quadrature('bounded', pdfFlat, 0, Lflat, 16);
            tc.verifyEqual(sum(wFlat .* qFlat), Lflat/2, 'RelTol', 1e-8);
            tc.verifyEqual(sum(wFlat), 1, 'AbsTol', 1e-8);
        end

        % ----------------------------------------------------------- trapz
        function testTrapzNodesEqualInput(tc)
            tc.verifyEqual(tc.qTrapz, tc.lvTrapz(:));
        end

        function testTrapzWeightsSumToOne(tc)
            tc.verifyEqual(sum(tc.wTrapz), 1, 'AbsTol', 1e-10);
        end

        function testTrapzWeightsPositive(tc)
            % fv here comes from a positive pdf, so weights must stay positive
            tc.verifyGreaterThanOrEqual(tc.wTrapz, 0);
        end

        function testTrapzReproducesNormalizedTrapzIntegral(tc)
            % sum(w .* g(q)) must equal the trapz-based normalized
            % expectation of g under fv on the same grid, for ANY g,
            % since this is just algebraic re-weighting of trapz.
            testFuns = {@(x) x, @(x) x.^2, @(x) sin(x/100)};
            for i = 1:numel(testFuns)
                g = testFuns{i};
                expected = trapz(tc.lvTrapz, tc.fvTrapz .* g(tc.lvTrapz)) / ...
                           trapz(tc.lvTrapz, tc.fvTrapz);
                actual = sum(tc.wTrapz .* g(tc.qTrapz));
                tc.verifyEqual(actual, expected, 'RelTol', 1e-10);
            end
        end

        % --------------------------------------------------------- general
        function testUnknownModeThrowsError(tc)
            tc.verifyError(...
                @() build_quadrature('not_a_real_distribution', 1, 2, 3), ...
                'build_quadrature:unknown_mode');
        end

        function testAllModesReturnColumnVectorsOfMatchingLength(tc)
            tc.verifyEqual(size(tc.qGamma), size(tc.wGamma));
            tc.verifyEqual(size(tc.qLogN), size(tc.wLogN));
            tc.verifyEqual(size(tc.qBounded), size(tc.wBounded));
            tc.verifyEqual(size(tc.qTrapz), size(tc.wTrapz));

            tc.verifyEqual(numel(tc.qGamma), tc.nNodes);
            tc.verifyEqual(numel(tc.qLogN), tc.nNodes);
            tc.verifyEqual(numel(tc.qBounded), tc.nBounded);
            tc.verifyEqual(numel(tc.qTrapz), numel(tc.lvTrapz));
        end

    end
end


function p = lognPdf(x, mu, sigma)
% Self-contained lognormal pdf (avoids a Statistics Toolbox dependency
% for the test fixtures -- build_quadrature.m itself has no such
% dependency, so the tests shouldn't introduce one either).
    p = 1 ./ (x .* sigma .* sqrt(2*pi)) .* exp(-(log(x) - mu).^2 ./ (2*sigma^2));
end


function c = lognCdf(x, mu, sigma)
% Self-contained lognormal cdf, via the standard erf-based formula.
    c = 0.5 * (1 + erf((log(x) - mu) ./ (sigma*sqrt(2))));
end