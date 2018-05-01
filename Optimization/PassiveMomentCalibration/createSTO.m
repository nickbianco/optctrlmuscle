function createSTO(array, col_labels, storage_fpath, in_degrees)
% Saves an array to an OpenSim Storage file.
% 
%     Parameters
%     ----------
%     array : data array (size: n_rows x n_cols)
%     col_labels: labels for each column in data array
%     storage_fpath : str
%     in_degrees : bool

n_rows = size(array,1);
n_cols = size(array,2);
if n_cols ~= length(col_labels)
   error('Number of columns and number of labels must be consistent.') 
end

f = fopen(storage_fpath, 'w');
fprintf(f, '%s\n', storage_fpath);
fprintf(f, 'version=1\n');
fprintf(f, 'nRows=%i\n', n_rows);
fprintf(f, 'nColumns=%i\n', n_cols);
inDegrees = 'no';
if in_degrees
    inDegrees = 'yes';
else
fprintf(f, 'inDegrees=%s\n', inDegrees);
fprintf(f, 'endheader\n');

for i = 1:n_cols
    col_label = col_labels{i};
    if i ~= 0
        fprintf(f, '\t');
    end
    fprintf(f, '%s', col_label);
end
fprintf(f, '\n');

for i = 1:n_rows
    for j = 1:n_cols
        if j ~= 0
            fprintf(f, '\t');
        end
        fprintf(f, '%f', array(i,j));
    end
    fprintf(f, '\n');
end

fclose(f);

end