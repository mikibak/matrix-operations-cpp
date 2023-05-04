clc
clear all
close all

ErrorJacobi = csvread('/home/mikibak/matrix-operations-cpp/plots/error0.csv')
ErrorGaussSeidel = csvread('/home/mikibak/matrix-operations-cpp/plots/error1.csv')

semilogy(1:20, ErrorJacobi(1:20));
hold on
semilogy(1:20, ErrorGaussSeidel(1:20));
title("B: Residual error norm");
ylabel("error");
xlabel("iteration");
legend('Jacobi', 'GaussSeidel');
hold off
saveas(gcf, '/home/mikibak/matrix-operations-cpp/plots/errors.png');