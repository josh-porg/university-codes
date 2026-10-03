%{
t_start = 150000; t_end = 169999; time_step = 1;
t_range = t_start : time_step : t_end;
%}

t_start = 150; t_end = 169; time_step = .5;
t_range = t_start : time_step : t_end;

s = zeros(length(t_range),1);
p = zeros(length(t_range),1);

mm=0;
for i = t_range
    mm = mm+1;
    s(mm,1) =  i;
end

parfor mm = 1:numel(t_range)
    i=t_range(mm);
    %disp(i);
    p(mm,1) =  i;
end


