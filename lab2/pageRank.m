% Clear all variables and close all figures
clear all
close all

% Set the format of the output
format compact

% Call the function
[numer_indeksu, Edges, I, B, A, b, r] = page_rank();
plot_PageRank(r)
function [numer_indeksu, Edges, I, B, A, b, r] = page_rank()
    numer_indeksu = 193410;
    L1 = 1; %przedostatnia cyfra
    L2 = 4; %3 od konca cyfra
    Edges = [1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 7, 8;
             4, 6, 3, 4, 5, 5, 6, 7, 5, 6, 4, 6, 8, 4, 7, 6, 2];
    N = 8;
    I = speye(N);
    B = sparse(Edges(2,:), Edges(1,:), 1, N, N);
    A = spdiags(1./sum(B, 1)', 0, N, N);
    disp(A);
    d = 0.85;
    b = (1 - d) / N * ones(N, 1);
    M = I-(d*B*A);
    r = M \ b;
end

function plot_PageRank(r)
    figure('Name', '193410')
    bar(r);
    xlabel('Numer strony');
    ylabel('Szansa na dostanie sie do tej strony');
    title('PageRank chart');
print -dpng zadanie7.png 
end
