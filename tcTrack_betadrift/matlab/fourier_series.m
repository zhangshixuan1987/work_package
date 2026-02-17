function F = fourier_series(m)
%A fourier series with random phase, whose amplitude is the square root of
%the observed variance
%   
N = 15; %total number of waves retained
T = 15; % time scale corresponding to the period of the lowest frequency wave

for n = 1:N
    x_1n(n) = rand;
end


for  i = 1:150
    t(i) = 0.1* i;
    term_1 = 1.2910;
    term_2 = 0.0;
    % 
    for n = 1:N
        x_in = rand;
        term_2 = term_2 + n^(-1.5)*sin(2*3.14*(n * t(i) / T + x_1n(n)));
        %term_2 = term_2 + n^(-1.5)*6*sin(2*3.14* x_in);
    end 
    F(i) = term_1 * term_2;
end

end

