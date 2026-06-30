classdef TestBuildMatrix < matlab.unittest.TestCase
    properties
        Lmax = 24;
        L2; Gxz; iLy; Gyz; iLx;
    end
    methods (TestMethodSetup)
        function buildMatrices(tc)
            [tc.L2, tc.Gxz, tc.iLy, tc.Gyz, tc.iLx] = ...
                build_matrix(tc.Lmax, 'store', false);
        end
    end

    methods (Test)
        % --- L2 ---
        function testL2Diagonal(tc)
            % L2 must be diagonal with l*(l+1) entries
            tc.verifyEqual(tc.L2, diag(diag(tc.L2)));
        end

        function testL2Values(tc)
            tc.verifyEqual(full(tc.L2(idx(0,0,tc.Lmax), idx(0,0,tc.Lmax))), 0);
            tc.verifyEqual(full(tc.L2(idx(2,0,tc.Lmax), idx(2,0,tc.Lmax))), 6);
            tc.verifyEqual(full(tc.L2(idx(4,0,tc.Lmax), idx(4,0,tc.Lmax))), 20);
        end

        function testL2PositiveSemiDefinite(tc)
            ev = eig(full(tc.L2));
            tc.verifyGreaterThanOrEqual(min(ev), -1e-12);
        end

        % --- iLy ---
        function testLyAntiHermitian(tc)
            % iLy should be anti-Hermitian: iLy + iLy' = 0
            tc.verifyEqual(full(tc.iLy + tc.iLy'), ...
                zeros(size(tc.iLy)), 'AbsTol', 1e-12);
        end

        function testLyNoSineCosineBlock(tc)
            % iLy should NOT mix positive and negative m
            % i.e. (l,m>0) <-> (l,m'<0) entries must be zero
            tc.verifyEqual(full(tc.iLy(idx(2,1,tc.Lmax), idx(2,-2,tc.Lmax))), 0);
            tc.verifyEqual(full(tc.iLy(idx(2,2,tc.Lmax), idx(2,-1,tc.Lmax))), 0);
        end

        % --- iLx ---
        function testLxAntiHermitian(tc)
            tc.verifyEqual(full(tc.iLx + tc.iLx'), ...
                zeros(size(tc.iLx)), 'AbsTol', 1e-12);
        end

        function testLxMixesSineCosine(tc)
            % iLx SHOULD mix positive and negative m
            tc.verifyNotEqual(full(tc.iLx(idx(2,1,tc.Lmax), idx(2,-2,tc.Lmax))), 0);
        end

        % --- Gxz ---
        function testGxzNoSineCosineBlock(tc)
            % Gxz couples same sign m only (like iLy)
            tc.verifyEqual(full(tc.Gxz(idx(2,1,tc.Lmax), idx(2,-2,tc.Lmax))), 0);
        end

        % --- Gyz ---
        function testGyzMixesSineCosine(tc)
            % Gyz couples opposite sign m (like iLx)
            tc.verifyNotEqual(full(tc.Gyz(idx(2,1,tc.Lmax), idx(2,-2,tc.Lmax))), 0);
        end

        % --- Rotation/Reflection consistency ---
        function testGyzIsReflectionOfGxz(tc)
            % Gyz should equal T*Gxz*T' where T is the (y=-x)-reflection
            N = size(tc.L2, 1);
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
            Gyz_rotated = T * tc.Gxz * T';
            tc.verifyEqual(full(Gyz_rotated), full(tc.Gyz), 'AbsTol', 1e-10);
        end

        function testiLxIsRotationOfiLy(tc)
            % iLx should equal T*iLy*T' where T is the z-rotation by pi/2
            N = size(tc.L2, 1);
            T = speye(N);
            for l = 0:2:tc.Lmax
                for m = 1:l
                    ii = idx(l,  m, tc.Lmax);
                    jj = idx(l, -m, tc.Lmax);
                    c = cos(m*pi/2); s = sin(m*pi/2);
                    % 2x2 rotation block for (+m, -m) pair
                    T(ii,ii) =  c; T(ii,jj) = s;
                    T(jj,ii) = -s; T(jj,jj) = c;
                end
            end
            iLx_rotated = T * tc.iLy * T';
            tc.verifyEqual(full(iLx_rotated), full(tc.iLx), 'AbsTol', 1e-10);
        end


        % --- Commutator: [L2, iLy] = 0 ---
        function testL2CommuteswithLy(tc)
            comm = tc.L2*tc.iLy - tc.iLy*tc.L2;
            tc.verifyEqual(full(comm), zeros(size(comm)), 'AbsTol', 1e-10);
        end
        
        function testL2CommuteswithLx(tc)
            comm = tc.L2*tc.iLx - tc.iLx*tc.L2;
            tc.verifyEqual(full(comm), zeros(size(comm)), 'AbsTol', 1e-10);
        end
    end
end