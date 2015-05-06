% By Ata Mahjoubfar at Jalali Lab, UCLA
% value_finder finds closest elements and their indexes in a
% matrix to a particular value or in two matrixes.
% Here are forms of the funtion application:
% [closest_value_index, closest_value] = value_finder(input_matrix, search_value);
% closest_value_index = value_finder(input_matrix, search_value);
% [closest_value_index, closest_value] = value_finder(input_matrix_1, input_matrix_2);
function [OUTPUT1,OUTPUT2] = value_finder(INPUT1,INPUT2)
proximity = abs(INPUT1-INPUT2);
OUTPUT1 = find(min(proximity) == proximity);
OUTPUT2 = INPUT1(OUTPUT1);