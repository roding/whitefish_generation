function val = read_xml_key(file_string, key_name, fmt)

pattern_before = ['<' key_name '>'];
pattern_after = ['</' key_name '>'];
ind_start = strfind(str, pattern_before) + numel(pattern_before);
ind_end = strfind(str, pattern_after) - 1;

ind = ind_start:ind_end;
key_string = file_string(ind);

val = nan;

if fmt == 'string'
    val = key_string;
elseif fmt == 'scalar'
    val = str2double(key_string);
elseif fmt == 'array'
    ind_separator = strfind(substr, ',');
    val = zeros(numel(ind_separator) + 1, 1);
    ind_separator = [0 ind_separator numel(substr)+1];
    for i = 1:numel(ind_separator)-1
        val(i) = str2double(key_string(ind_separator(i)+1:ind_separator(i+1)-1));
    end
end
    
end

