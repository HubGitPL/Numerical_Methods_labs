clear all
close all
format compact

n_max = 200;
a = 10;
r_max = a/2;

[circles, index_number, circle_areas, rand_counts, counts_mean] = generate_circles(a, r_max, n_max);
plot_circle_areas(circle_areas);
plot_counts_mean(counts_mean);
plot_circles(a, circles, index_number); 
print -dpng zadanie1.png 

function [circles, index_number, circle_areas, rand_counts, counts_mean] = generate_circles(a, r_max, n_max)
    index_number = 193410; % your index number
    L1 = mod(index_number, 10);
    circles = rand(n_max, 3);
    circle_areas = rand(n_max, 1);
    rand_counts = rand(n_max, 1);
    counts_mean = rand(n_max, 1);
    for i = 1:n_max
        flag = true;
        counter = 0;
        while flag
            counter = counter+1;
            R = rand * r_max;
            X = rand * (a - 2*R) + R;
            Y = rand * (a - 2*R) + R;
            flag = false; 
            for j = 1:i-1
                d = sqrt((X-circles(j, 1))^2 + (Y-circles(j, 2))^2); 
                if d < R + circles(j, 3)
                    flag = true;
                    break;
                end
            end
        end
        circles(i, 1) = X;
        circles(i, 2) = Y;
        circles(i, 3) = R;
        circle_areas(i) = pi*R^2;
        rand_counts(i) = counter;
        counts_mean(i) = sum(rand_counts(1:i))/i;
    end
    suma = 0;
    for j = 1:n_max
        suma = suma + circle_areas(j);
        circle_areas(j) = suma;
    end
end

function plot_circles(a, circles, index_number)
    figure; hold on; 
    axis equal;
    axis([0 a 0 a]);

    for i = 1:size(circles, 1)
        R = circles(i, 3); 
        X = circles(i, 1); 
        Y = circles(i, 2); 
        plot_circle1(R, X, Y); 
    end
    hold off; 
end

function plot_circle_areas(circle_areas)
    figure('Name', '193410')
    plot(circle_areas);
    xlabel('Ilosc okregow ');
    ylabel('Powierzechnia zajmowana');
    title('Zmiana powierzchni zajmowanej wraz z dodawaniem okregow');
    print -dpng zadanie3.png; 
end

function plot_counts_mean(counts_mean)
    figure('Name', '193410')
    plot(counts_mean);
    xlabel('Ilosc okregow ');
    ylabel('Średnia liczba losowań potrzebna do znalezienia i okregow');
    title('Zmiana średniej liczby losowań wraz z dodawaniem okregow');
    print -dpng zadanie5.png; 
end

function plot_circle1(R, X, Y)
    % R - promień okręgu
    % X - współrzędna x środka okręgu
    % Y - współrzędna y środka okręgu
    theta = linspace(0,2*pi);
    x = R*cos(theta) + X;
    y = R*sin(theta) + Y;
    plot(x,y)
end
