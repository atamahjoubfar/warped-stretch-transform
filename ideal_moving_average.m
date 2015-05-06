function output = ideal_moving_average(input, span)
% Ideal moving average function
% By Ata Mahjoubfar
% November 2011
[size_width, size_length] = size(input);
if size_width > 1
    if size_length > 1
        error('This function works does not work on matrices.')
    end
    input = input.';
    flip_back = 1;
else
    flip_back = 0;
end
output = zeros(size(input));
extra_input = [fliplr(input(2:1+floor(span/2))), input, fliplr(input(end-floor(span/2):end-1))];
for index = 1 : length(input)
    output(index) = mean(extra_input(index:index+span-1));
end
if flip_back == 1
    output = output.';
end
end