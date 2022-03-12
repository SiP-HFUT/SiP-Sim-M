%%%% Use matrix transformation method to calculate analytical expressions of various types of MRR transfer matrices
clc
close all
clearvars -except

%% M = tm_admrr1 (k0, k1, ...)
% syms a0 b0 a3 b3 
% M =  sym('M', [2,2]);
% x3 = [b3;a3];
% x0 = [a0;b0];
% eqs = M*x3 == x0;
% 
% % % solve the transfer matrix, M, which is defiend as '[a3; a0] = M * [b3; b0]'
% % S = solve(eqs,[a3 a0]);
% % a3_anal = simplify(S.a3);
% % a0_anal = simplify(S.a0);
% 
% % solve the transfer matrix, M, which is defiend as '[b3; b0] = M * [a3; a0]'
% S = solve(eqs,[b3 b0]);
% b3_anal = simplify(S.b3);
% b0_anal = simplify(S.b0);



%% M = tm_admrr1 (k1, k0, ...)
syms a0 b0 a3 b3 
M =  sym('M', [2,2]);
x3 = [b3;a3];
x0 = [a0;b0];
eqs = M*x0 == x3;

% solve the transfer matrix, M, which is defiend as '[a0; a3] = M * [b0; b3]'
S = solve(eqs,[a3 a0]);
a3_anal = simplify(S.a3);
a0_anal = simplify(S.a0);
