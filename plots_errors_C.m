clc
clear all
close all

ErrorJacobi = csvread('/home/mikibak/matrix-operations-cpp/plots/error2.csv')
ErrorGaussSeidel = csvread('/home/mikibak/matrix-operations-cpp/plots/error3.csv')
hold on
semilogy(1:1000, ErrorJacobi);
semilogy(1:1000, ErrorGaussSeidel);
legend('Jacobi', 'GaussSeidel');
title("C: Residual error norm");
ylabel("error");
xlabel("iteration");
hold off
saveas(gcf, '/home/mikibak/matrix-operations-cpp/plots/errors_C.png');