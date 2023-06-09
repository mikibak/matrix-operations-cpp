clc
clear all
close all

N = [100, 500, 1000, 2000, 3000];

GaussSeidel = csvread('GaussSeidel.csv')
Jacobi = csvread('Jacobi.csv')

semilogy(N, GaussSeidel);
hold on
semilogy(N, Jacobi);
legend('Jacobi', 'GaussSeidel');
title("B:Time of execution");
ylabel("time [s]");
xlabel("Matrix size");
hold off
saveas(gcf, '/home/mikibak/matrix-operations-cpp/plots/iteracje.png');
