function val = read_key(file_string, key_name, fmt)

pattern_before = ['<' key_name '>'];
pattern_after = ['</' key_name '>'];
ind_start = strfind(file_string, pattern_before) + numel(pattern_before);
ind_end = strfind(file_string, pattern_after) - 1;

ind = ind_start:ind_end;
key_string = file_string(ind);

val = nan;

if isequal(fmt, 'string')
    val = key_string;
elseif isequal(fmt, 'scalar')
    val = str2double(key_string);
elseif isequal(fmt, 'array')
    ind_separator = strfind(key_string, ',');
    val = zeros(numel(ind_separator) + 1, 1);
    ind_separator = [0 ind_separator numel(key_string)+1];
    for i = 1:numel(ind_separator)-1
        val(i) = str2double(key_string(ind_separator(i)+1:ind_separator(i+1)-1));
    end
end
    
end

