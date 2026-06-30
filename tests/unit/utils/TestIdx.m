classdef TestIdx < matlab.unittest.TestCase
    methods (Test)
        function testFirstBlock(tc)
            % (0,0) -> 1
            tc.verifyEqual(idx(0,0,4), 1);
        end

        function testSecondBlock(tc)
            % (2,-2) -> 2
            tc.verifyEqual(idx(2,-2,4), 2);
            tc.verifyEqual(idx(2,-1,4), 3);
            tc.verifyEqual(idx(2, 0,4), 4);
            tc.verifyEqual(idx(2, 1,4), 5);
            tc.verifyEqual(idx(2, 2,4), 6);
        end

        function testOutOfBounds(tc)
            tc.verifyEqual(idx(2, 3,4), 0);   % m > l
            tc.verifyEqual(idx(2,-3,4), 0);   % m < -l
            tc.verifyEqual(idx(6, 0,4), 0);   % l > Lmax
        end

        function testInverseWithLmdx(tc)
            % idx and lmdx should be inverses
            Lmax = 8;
            for l = 0:2:Lmax
                for m = -l:l
                    i = idx(l,m,Lmax);
                    [lr, mr] = lmdx(i);
                    tc.verifyEqual(lr, l);
                    tc.verifyEqual(mr, m);
                end
            end
        end
        
        function testMonotone(tc)
            % indices should be strictly increasing with l, then m
            Lmax = 8;
            prev = 0;
            for l = 0:2:Lmax
                for m = -l:l
                    i = idx(l,m,Lmax);
                    tc.verifyGreaterThan(i, prev);
                    prev = i;
                end
            end
        end
    end
end