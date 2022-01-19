function [dr,thr] = Resolve(M,dr_char,thr_char)
M1_1 = M(1,1); M1_2 = M(1,2);
M2_1 = M(2,1); M2_2 = M(2,2);
eval(['dr=' dr_char ';']);
eval(['thr=' thr_char ';']);
end