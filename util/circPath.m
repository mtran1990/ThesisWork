function loc = circPath(A, freq, t)

    loc = zeros(2,length(t));

    loc(1,:) = A*cos(2*pi*freq*t);
    loc(2,:) = A*sin(2*pi*freq*t);

end