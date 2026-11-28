clc
clear

%global var_per_line var_per_line_grid fname_grid line_skip var_target var_cc var_node node_num element_num zone_type line_skip_grid input x_min x_max limit_value

input=100; % 1, interp from plt file 2, CI volume 3, Experiment Data 100, Tecplot ASCII (Please use this option)

var_per_line = 4; var_per_line_grid = 5; line_skip = 15; %number of header lines to skip 
line_skip_grid = 9; 
var_cc = 1:8; %total number of variables in the data
var_select = 1:8; %the # of variables you choose to read
%var_select = 1;
var_node = [];

element_num = 38523; node_num = 39065; %number of elements and nodes in the data, get it from the grid file
case_name = ''; gridName = 'grid.dat';


t_start = 150000; t_end = 169998; time_step = 1;
%t_start = 150000; t_end = 150000; time_step = 1;

t_range = t_start : time_step : t_end;

dimension_option=2; zone_num=1; %not entirely sure what dimention option is

m_folder = ['' case_name]; fname_raw = [m_folder 'test_file_ncons_']; %this (fname raw) may need to change


fname_grid = ['' gridName];

%TODO:I do not want to import more data than needed. as such i sould
%probably only select 1 or two variables
var_target = var_select;

var_name=['P    ';'U    ';'V    ';'T    ';'Y_CH4';'Y_O2 ';'Y_H2O';'Y_CO2'];

switch dimension_option
   case(2)
       zone_type = 'FEQuadrilateral';
   case(3)
       zone_type = 'FEBrick'
end

%as far as i can tell this configures vaiables for reading and writing the
%data
[read_format,read_format_xyz,out_format,out_format_xyz,read_format_c2n,out_format_c2n,prec] = InOutFormat(var_per_line,var_per_line_grid,dimension_option);

%I belive this imports the grid file such that xyz, c2n do things? idk
[xyz,c2n] = importGridFile(fname_grid,line_skip_grid,read_format_xyz,read_format_c2n,dimension_option,element_num,node_num);

%this is just assigning things i belive
ndim = dimension_option; num_c = element_num; var_num = length(var_select); It = length(t_range);

% I think this to be the part in which we create an empty snapshot style
% matrix for assigning to in the upcomeing for-loop
p = zeros(It,num_c*var_num);


parfor mm = 1:numel(t_range)
    i = t_range(mm);
    
    file_name = [fname_raw num2str(i) '.dat'];
    out = importTecASCIIdata(file_name,line_skip,read_format,var_target,var_cc,var_node,node_num,element_num,c2n,ndim);

%     for ivar = 1 : var_num
%         p(mm,(1+(ivar-1)*num_c):(ivar*num_c)) = transpose(out(:,ivar));
%     end
    v = zeros(1,num_c*var_num);
    for ivar = 1 : var_num
        v((1+(ivar-1)*num_c):(ivar*num_c)) = transpose(out(:,ivar));
    end
    p(mm,:) = v;

    filename_m = [m_folder 'flowfiled_' num2str(i) '.dat'];
    OutputTecASCIIdata(filename_m,2,dimension_option,var_name,node_num,element_num,zone_type,out_format_xyz,out_format,out_format_c2n,xyz,out,c2n);    
    disp(i);
end

save('combustorData150000_169999.mat', 'p', '-v7.3')
