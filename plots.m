N = [10, 50, 100, 200, 300];

GaussSeidel = csvread('GaussSeidel.csv')
Jacobi = csvread('Jacobi.csv')

hold on
plot(N, GaussSeidel);
plot(N, Jacobi);
legend('GaussSeidel', 'Jacobi');
title("Time of execution");
ylabel("time [s]");
xlabel("iterations");
hold off
saveas(gcf, '/home/mikibak/matrix-operations-cpp/plots/iteracje.png');