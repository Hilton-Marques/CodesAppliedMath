clc;
clear;
close all;

[A,b,Aeq,beq]=vert2lcon(eye(3));
V=lcon2vert(A,b,Aeq,beq);