function out = importTecASCIIdata(file_name,line_skip,read_format,var,var_cc,var_node,node_num,element_num,c2n,ndim)

    disp(['file name:' file_name]);

    fid = fopen(file_name); %this opens the file called file_name
    
    %skip over line_skip number of header lines
    for k = 1 : line_skip
        fgetl(fid);
    end
    
    var_num_cc = length(var_cc); %total number of variable in the data (var_cc is an array of vriable positions?)
    
    var_num_node = length(var_node); %this is an empty array in the sample code not sure whats going on here. in this case lenght will be zero
    
    var_num = length(var); %var is an array of the positions of the variables (e.g 1-8) that a want to exctract
    
    %this whole if statment assigns a value to N and assigns a matrix of
    %zeros to out
    if ndim < 3 %ndim is asigned in sampleO as dimention option and is assigned as 2
        %this says N is (total number of vars) * (number of elements) + (?) * (number of nodes)
        % I am unsure as to what N is but it seems to be some form of size
        % that is i think the sum total of elements and nodes in the data
        N = ceil(var_num_cc * element_num + var_num_node * node_num); %node num is number of nodes in the data
        if isempty(var_cc) == 0
            out = zeros(element_num,var_num);
        else
            out = zeros(node_num,var_num);
        end
    else
        N = ceil(var_num_cc * element_num + var_num_node * node_num);
        out = zeros(element_num,var_num);
    end
    
    %reads data from the text file, found by fid, using format specifications read_format,
    %and reads at most N elements from the file 
    import = fscanf(fid,read_format,N); 
    
%    rs = 1;
    
    for vp = 1 : var_num
        
        if ismember(vp,var_cc)
            increment = element_num;
        else
            increment = node_num;
        end
 
%        if vp == var(1)
         rs = 1 + increment*(var(vp)-1);
%        end
       
        re = rs + increment - 1;
        
        range = rs : re ;
        
        read_in = import(range);
        
        if ndim < 3
            out(:,vp) = read_in;
        else
            if ismember(var(vp),var_cc)
                out(:,vp) = read_in;
            else
                out(:,vp) = CalCellCenterValue(read_in,c2n,element_num);
            end
        end
        
%        rs = rs + increment;
    end
    
    fclose(fid);
end

function out = CalCellCenterValue(read_in,c2n,element_num)

    out = zeros(element_num,1);
    
    for ni = 1 : element_num
       
        nodes = c2n(ni,:);
        out(ni) = mean(read_in(nodes));
        
    end

end
