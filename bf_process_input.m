function conf = bf_process_input(input_file)
% Process the input file

conf = struct('target_dir','.','scatter','segment','order',1,'qsm_mat_file',[],...
    'qsm_cyl_table',[],'qsm_br_table',[],'segment',[],'branch',[],...
    'ssm_fun',[],'ga_init_lb',[],'ga_init_ub',[],'ga_lb',[],'ga_ub',[],'ga_int_con',[],...
    'ga_pop_size',[],'ga_gens',[],'ga_stall',[],'ga_elite',[],...
    'dt_stat1d',1,'dt_dirs',100,'dt_scale',false,'dt_w',[],...
    'ssm_tree_fun','@()read_mtg(''out.mtg'')','ga_out_dat','gaOut.dat',...
    'ga_use_par',0,'ga_rng',[],'ssm_fun_best',[],'movie',true,'qsm_merge',false,'ga_out_fun',[],...
    'ga_tol_fun',1e-6,'ga_multi',0);

fprintf('Start processing ''%s''...',input_file);

fid = fopen(input_file,'r');

line = fgetl(fid);

% Get continued lines
line = get_continued(fid,line);

while(ischar(line))
    if(~isempty(line) && line(1) ~= '#')
        C = textscan(line,'%s','delimiter','=');% process line
        left = C{1}{1}; left = strtrim(left);% left with no spaces
        right = C{1}{2}; right = strtrim(right);% right, no spaces
        %fprintf('Read: `%s'' = `%s''\n',left,right);
        conf = write_conf(left,right,conf);% save to the conf
    end
    line = fgetl(fid);
    line = get_continued(fid,line);
end

fclose(fid);
fprintf('Done.\n');
end

function line = get_continued(fid,line)
% Check and get a continued line.
k = strfind(line,'...');
while( ischar(line) && ~isempty(k) )
    line = line(1:k(1)-1);
    line = cat(2,line,strtrim(fgetl(fid)));
    k = strfind(line,'...');
end
end

function conf = write_conf(left,right,conf)
% Write the conf structure
if(strcmpi(left,'target_dir'))
    conf.target_dir = right;
elseif(strcmpi(left,'scatter'))
    C = textscan(right,'%s','delimiter',',');
    conf.scatter = C{1};
elseif(strcmpi(left,'order'))
    C = textscan(right,'%s','delimiter',',');
    conf.order = [];% to erase the default value
    for ii = 1:length(C{1})
        conf.order = cat(2,conf.order,str2double(strtrim(C{1}{ii})));
    end
elseif(strcmpi(left,'qsm_mat_file'))
    conf.qsm_mat_file = right;
elseif(strcmpi(left,'qsm_cyl_table'))
    conf.qsm_cyl_table = right;
elseif(strcmpi(left,'qsm_br_table'))
    conf.qsm_br_table = right;
elseif(strcmpi(left,'qsm_merge'))
    conf.qsm_merge = str2double(right) > 0;
elseif(strcmpi(left,'segment'))
    C = textscan(right,'%s','delimiter',',');
    for ii = 1:length(C{1})
        conf.segment = cat(2,conf.segment,str2double(strtrim(C{1}{ii})));
    end
elseif(strcmpi(left,'branch'))
    C = textscan(right,'%s','delimiter',',');
    for ii = 1:length(C{1})
        conf.branch = cat(2,conf.branch,str2double(strtrim(C{1}{ii})));
    end
elseif(strcmpi(left,'ssm_fun'))
    conf.ssm_fun = right;
elseif(strcmpi(left,'ssm_fun_best'))
    conf.ssm_fun_best = right;
elseif(strcmpi(left,'ga_init_lb'))
    C = textscan(right,'%s','delimiter',', ','MultipleDelimsAsOne',true);
    for ii = 1:length(C{1})
        conf.ga_init_lb = cat(2,conf.ga_init_lb,str2double(strtrim(C{1}{ii})));
    end
elseif(strcmpi(left,'ga_init_ub'))
    C = textscan(right,'%s','delimiter',', ','MultipleDelimsAsOne',true);
    for ii = 1:length(C{1})
        conf.ga_init_ub = cat(2,conf.ga_init_ub,str2double(strtrim(C{1}{ii})));
    end
elseif(strcmpi(left,'ga_out_dat'))
    conf.ga_out_dat = right;
elseif(strcmpi(left,'ga_lb'))
    C = textscan(right,'%s','delimiter',', ','MultipleDelimsAsOne',true);
    for ii = 1:length(C{1})
        conf.ga_lb = cat(2,conf.ga_lb,str2double(strtrim(C{1}{ii})));
    end
elseif(strcmpi(left,'ga_ub'))
    C = textscan(right,'%s','delimiter',', ','MultipleDelimsAsOne',true);
    for ii = 1:length(C{1})
        conf.ga_ub = cat(2,conf.ga_ub,str2double(strtrim(C{1}{ii})));
    end
elseif(strcmpi(left,'ga_int_con'))
    C = textscan(right,'%s','delimiter',',');
    for ii = 1:length(C{1})
        conf.ga_int_con = cat(2,conf.ga_int_con,str2double(strtrim(C{1}{ii})));
    end
elseif(strcmpi(left,'ga_out_fun'))
    C = textscan(right,'%s','delimiter',',');
    for ii = 1:length(C{1})
        conf.ga_out_fun = [conf.ga_out_fun {str2func(strtrim(C{1}{ii}))}];
    end
elseif(strcmpi(left,'ga_pop_size'))
    conf.ga_pop_size = str2double(right);
elseif(strcmpi(left,'ga_gens'))
    conf.ga_gens = str2double(right);
elseif(strcmpi(left,'ga_stall'))
    conf.ga_stall = str2double(right);
elseif(strcmpi(left,'ga_elite'))
    conf.ga_elite = str2double(right);
elseif(strcmpi(left,'ga_tol_fun'))
    conf.ga_tol_fun = str2double(right);
elseif(strcmpi(left,'ga_use_par'))
    conf.ga_use_par = str2double(right) > 0;
elseif(strcmpi(left,'ga_rng'))
    conf.ga_rng = uint64(str2double(right));
elseif(strcmpi(left,'ga_multi'))
    conf.ga_multi = str2double(right) > 0;
elseif(strcmpi(left,'dt_stat1d'))
    conf.dt_stat1d = str2double(right);
elseif(strcmpi(left,'dt_dirs'))
    conf.dt_dirs = str2double(right);
elseif(strcmpi(left,'dt_scale'))
    conf.dt_scale = str2double(right) > 0;
elseif(strcmpi(left,'dt_w'))
    C = textscan(right,'%s','delimiter',', ','MultipleDelimsAsOne',true);
    for ii = 1:length(C{1})
        conf.dt_w = cat(2,conf.dt_w,str2double(strtrim(C{1}{ii})));
    end
elseif(strcmpi(left,'movie'))
    conf.movie = str2double(right) > 0;
else
    fprintf('Unknown line: %s = %s\n',left,right);
end

end