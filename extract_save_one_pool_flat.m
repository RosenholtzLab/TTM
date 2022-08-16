function [texture_mat]=extract_save_one_pool_flat(pool_descriptor, pool_id,texture_mat, with_pixel_hist,save_names)
    % file = './000000000139_X300_YNaN_2022-07-25_uniform_60olap_rho_1111111118/poolingDescriptor.mat'
    % descriptor = load(file);
    % texture = descriptor.textureDescriptor{1}.p;
    % savename = 'test_save_unnested.mat';
    % h5create('stat_data_new.h5','/pool1',[3 3]); %for some reason it like the tempdir better ...
    % h5write('stat_data_new.h5','/pool1',texture_mat)
    % D = h5read('stat_data_new.h5','/pool1');
    % insert data into pool 2:
    % h5write('stat_data_new.h5','/pool1',A, [2,1,1],[1 3 2302]);
    % experiment with compression or 'chunk' in h5create
    
    texture = pool_descriptor.p;
    cmu = pool_descriptor.cmu;
    cmat = pool_descriptor.cmat;
    
    % initialize stats
%     if with_pixel_hist == 1
%         texture_mat = NaN(1,3,2306);
%     else
%         texture_mat = NaN(1,3,1789);
%     end
    
    if save_names == 1
        stat_names = strings(15,1);
    else
        stat_names = NaN;
    end
    
    % index by pooling region
    for c=1:3
        fill_id = 1;
        total_size = 0;
        current_stat = 1;
        channel_stats = texture(c);
        F = fieldnames(channel_stats);
        for k=1:numel(F)
            stat = [channel_stats.(F{k})];
            if strcmp(F{k},'pixelHist')
                hist_elts = fieldnames(stat);
                for kp=1:numel(hist_elts)
                    elt = [stat.(hist_elts{kp})];
                    elt_flat = reshape(elt,[1,numel(elt)]);
                    total_size = total_size + numel(elt);
                    texture_mat(pool_id,c,fill_id:fill_id+numel(elt)-1) = elt_flat;
                    if save_names == 1 & c == 1
                        name = sprintf('%s_%s_size_%s_pyinds_%s_%s',F{k},hist_elts{kp},num2str(numel(elt)),num2str(fill_id-1),num2str(fill_id+numel(elt)-1));
                        stat_names(current_stat,:) = name;
                        current_stat = current_stat + 1;
                    end

                    fill_id = fill_id + numel(elt);
                    
                end
            elseif ~strcmp(F{k},'pixelLPHist')
                % I believe matlab uses column major ordering
                stat_flat = reshape(stat,[1,numel(stat)]);
                total_size = total_size + numel(stat);
                texture_mat(pool_id,c,fill_id:fill_id+numel(stat)-1) = stat_flat;
                if save_names == 1 & c == 1
                    name = sprintf('%s_size_%s_pyinds_%s_%s',F{k},num2str(numel(stat)),num2str(fill_id-1),num2str(fill_id+numel(stat)-1));
                    stat_names(current_stat,:) = name;
                    current_stat = current_stat + 1;
                end
                fill_id = fill_id + numel(stat);


            end
        
        end
        texture_mat(pool_id,c,fill_id:fill_id) = cmu(c);
        texture_mat(pool_id,c,fill_id+1:fill_id+3) = cmat(:,c)';
        if save_names == 1 & c == 1
            name = sprintf('%s_size_%s_pyinds_%s_%s','cmu',num2str(1),num2str(fill_id-1),num2str(fill_id));
            stat_names(current_stat,:) = name;
            current_stat = current_stat + 1;
            fill_id = fill_id + 1;

            name = sprintf('%s_size_%s_pyinds_%s_%s','cmat',num2str(3),num2str(fill_id-1),num2str(fill_id+3-1));
            stat_names(current_stat,:) = name;
            current_stat = current_stat + 1;
        end
        fill_id = fill_id + 1;
    end
end
    
    
    