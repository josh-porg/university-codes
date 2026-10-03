function [read_format,read_format_xyz,out_format,out_format_xyz,read_format_c2n,out_format_c2n,prec] = InOutFormat(var_per_line,var_per_line_grid,dimension_option) 

read_format = '%f'; read_format_xyz = '%f'; prec = '%20.15E'; out_format = prec;
for j = 1 : var_per_line-1
    read_format = [read_format ' %f'];
    out_format = [out_format ' ' prec];
end
read_format = [read_format '\n'];
out_format = [out_format '\n'];

for j = 1 : var_per_line_grid-1
    read_format_xyz = [read_format_xyz ' %f'];
end
read_format_xyz = [read_format_xyz '\n'];

out_format_xyz = out_format;

read_format_c2n = '%f'; out_format_c2n = '%6d';
for j = 1 : (2^dimension_option-1)
    read_format_c2n = [read_format_c2n ' %f'];
    out_format_c2n = [out_format_c2n ' %6d'];
end
read_format_c2n = [read_format_c2n '\n'];
out_format_c2n = [out_format_c2n '\n'];

end
