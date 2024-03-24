

%% 
N = 1000:1000:8000;
n = length(N);

vtime_direct = ones(1,n); 
for i = 1:n
    [A,b,x,vtime_direct(i),err_norm,index_number] = solve_direct(N(i));
end

plot_direct(N,vtime_direct);
%% 

N = 100;
[A,b,x,time_direct,err_norm,index_number] = solve_direct(N);


N = 100;
[A,b,M,bm,x,err_norm,time,iterations,index_number] = solve_Jacobi(N);

N = 100;
[A,b,M,bm,x,err_norm,time,iterations,index_number] = solve_Gauss_Seidel(N);

function [A,b,M,bm,x,err_norm,time,iterations,index_number] = solve_Gauss_Seidel(N)
% A - macierz rzadka z równania macierzowego A * x = b
% b - wektor prawej strony równania macierzowego A * x = b
% M - macierz pomocnicza opisana w instrukcji do Laboratorium 3 – sprawdź wzór (7) w instrukcji, który definiuje M jako M_{GS}
% bm - wektor pomocniczy opisany w instrukcji do Laboratorium 3 – sprawdź wzór (7) w instrukcji, który definiuje bm jako b_{mGS}
% x - rozwiązanie równania macierzowego
% err_norm - norma błędu rezydualnego rozwiązania x; err_norm = norm(A*x-b)
% time - czas wyznaczenia rozwiązania x
% iterations - liczba iteracji wykonana w procesie iteracyjnym metody Gaussa-Seidla
% index_number - Twój numer indeksu
index_number = 193410;
L1 = 0;
x = [];
M = [];
bm = [];
time = [];
iterations = [];
err_norm = 1;

[A, b] = generate_Matrix(N, L1);
D=diag(diag(A));
L=tril(A, -1);
U=triu(A, 1);
M = -(D + L)\U;
bm =(D + L)\b;

x=ones(N, 1);

tic;
for iterations=1:1000
    x=M*x + bm;
    err_norm=norm(A*x-b);
    if err_norm < 1E-12
        break
    end
end
time = toc;

err_norm =norm(A*x - b);

end

function [A,b,M,bm,x,err_norm,time,iterations,index_number] = solve_Jacobi(N)
% A - macierz z równania macierzowego A * x = b
% b - wektor prawej strony równania macierzowego A * x = b
% M - macierz pomocnicza opisana w instrukcji do Laboratorium 3 – sprawdź wzór (5) w instrukcji, który definiuje M jako M_J.
% bm - wektor pomocniczy opisany w instrukcji do Laboratorium 3 – sprawdź wzór (5) w instrukcji, który definiuje bm jako b_{mJ}.
% x - rozwiązanie równania macierzowego
% err_norm - norma błędu rezydualnego rozwiązania x; err_norm = norm(A*x-b)
% time - czas wyznaczenia rozwiązania x
% iterations - liczba iteracji wykonana w procesie iteracyjnym metody Jacobiego
% index_number - Twój numer indeksu
index_number = 193410;
L1 = 0;
x = [];
M = [];
bm = [];
time = [];
iterations = -1;
err_norm = 1;

[A, b] = generate_Matrix(N, L1);
L = tril(A, -1);
U = triu(A, 1);
D = diag(diag(A));

M = -inv(D) * (L + U);
bm = inv(D) * b;
x=ones(N, 1);

tic;
for iterations = 1:1000
    x=M*x + bm;
    err_norm = norm(A*x - b);
    if err_norm < 1E-12
        break;
    end
end
time = toc;
err_norm = norm(A*x-b);


end

%% 

function plot_direct(N,vtime_direct)
    % N - wektor zawierający rozmiary macierzy dla których zmierzono czas obliczeń metody bezpośredniej
    % vtime_direct - czas obliczeń metody bezpośredniej dla kolejnych wartości N
    plot(N, vtime_direct, 'o-b');
    title('Rozmiar macierzy vs Czas');
    xlabel('Rozmiar macierzy');
    ylabel('Czas');
    print -dpng zadanie2.png; 
end

function [A,b,x,time_direct,err_norm,index_number] = solve_direct(N)
% A - macierz z równania macierzowego A * x = b
% b - wektor prawej strony równania macierzowego A * x = b
% x - rozwiązanie równania macierzowego
% time_direct - czas wyznaczenia rozwiązania x
% err_norm - norma błędu rezydualnego rozwiązania x; err_norm = norm(A*x-b);
% index_number - Twój numer indeksu
x = [];
time_direct = [];
err_norm = 1;
index_number = 193410;
L1 = 0;
[A, b] = generate_Matrix(N, L1);
%disp(size(A));
%disp(size(b));
tic;
x=A\b;
time_direct = toc;
err_norm = norm(A*x-b);
%disp(x);
%disp(time_direct);
%disp(err_norm);
end


function [A,b] = generate_Matrix(N, convergence_factor)
    % A - macierz o rozmiarze NxN
    % b - wektor o rozmiarze Nx1
    % convergense_factor - regulacja elementów diagonalnych macierzy A, które wpływają
    %       na zbieżność algorytmów iteracyjnego rozwiązywania równania macierzowego

    if(convergence_factor<0 || convergence_factor>9)
        error('Wartość convergence_factor powinna być zawarta w przedziale [1,9]');
    end

    seed = 0; % seed - kontrola losowości elementów niezerowych macierzy A
    rng(seed); % ustawienie generatora liczb losowych

    A = rand(N, N);
    A = A - diag(diag(A)); % wyzerowanie głównej diagonalnej

    convergence_factor_2 = 1.2 + convergence_factor/10;
    diag_values = sum(abs(A),2) * convergence_factor_2;
    A = A + diag(diag_values); % nadanie nowych wartości na głównej diagonalnej

    % regulacja normy macierzy
    norm_Frobenius = norm(A,'fro');
    A = A/norm_Frobenius;

    b = rand(N,1);
end

