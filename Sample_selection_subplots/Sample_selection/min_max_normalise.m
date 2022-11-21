function J = min_max_normalise(I)

I = double(I);
minimum = min(I(:));
maximum = max(I(:));

J = (I - minimum) ./ (maximum - minimum);

