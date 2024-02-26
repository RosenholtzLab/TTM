function unNested=extract_save_one_pool(pool_descriptor,pool_id, savename)
    % file = './000000000139_X300_YNaN_2022-07-25_uniform_60olap_rho_1111111118/poolingDescriptor.mat'
    % descriptor = load(file);
    % texture = descriptor.textureDescriptor{1}.p;
    % savename = 'test_save_unnested.mat';
    
    texture = pool_descriptor.p;
    
    % loop to find what is struct and cell at first level
    % Fnames = fieldnames(texture);
    % for kf = 1:numel(Fnames)
    %     if isstruct([texture.(Fnames{kf})])
    %         disp(sprintf('%s is a struct',Fnames{kf}))
    %     elseif iscell([texture.(Fnames{kf})])
    %         disp(sprintf('%s is a cell',Fnames{kf}))
    %     end
    % end

    unNested = extract(texture);
    cmu = pool_descriptor.cmu;
    cmat = pool_descriptor.cmat;
    unNested.(sprintf('pool_%s_cmu',pool_id)) = cmu;
    unNested.(sprintf('pool_%s_cmat',pool_id)) = cmat;
    
    
    % Fn = fieldnames(unNested);
    % dim_array = zeros(length(Fn),1);
    %     for kF = 1:numel(Fn)
    %         % flatten row major order
    %         flat_vals = reshape(unNested.(Fn{kF}).',1,[]);
    %         size(flat_vals)
    %     end
    
    
    
    if isfile(savename)
    % add to the stats file if it exists
        save(savename, '-struct', 'unNested', '-append')
    else
     save(savename, '-struct', 'unNested')
    end
    

    %%
    function B = extract(A)
    B = struct();
    nestfun(A,num2str(0))
        function nestfun(D,C)
        % nestfun is recursive function that saves struct elements as an array in a new, un-nested struct.
        % https://www.mathworks.com/matlabcentral/answers/1627260-what-is-the-best-practice-to-recursively-extract-data-from-a-nested-structure-in-matlab
        if isstruct(D)
            F = fieldnames(D);
            for k = 1:numel(F)
                % need to put the struct values in an array
                elt = [D.(F{k})];
                if isstruct(elt)
                    % bind values to cell so the struct will not be stacked
                    elt = {D.(F{k})};
                end
                name = sprintf('pool_%s_%s_%s',pool_id,F{k},C);
                % write the values
                if ~isstruct(elt) & ~iscell(elt)
                    B.(name) = elt;
                % unpack the cell
                elseif iscell(elt)
                    nestcell(elt,F{k})
                end
            end
        end
        function nestcell(E,G)
        % nestcell is recursive function that find structs within cells
        if iscell(E)
            [x,y] = size(elt);
            for xi = 1:x
                for yi = 1:y
                    val = elt{xi,yi};
                    namec = sprintf('%s_%s_%s',G,num2str(xi),num2str(yi));
                    if isstruct(elt{xi,yi})
                        nestfun(elt{xi,yi},namec);
                    elseif iscell(elt{xi,yi})
                        nestcell(elt{xi,yi},namec);
                    else
                        disp(sprintf('not struct or cell!'))
                    end
                end
            end
        end
            end
        end 
    % disp(sprintf('done un-nesting!'))
    end
    end
    
    