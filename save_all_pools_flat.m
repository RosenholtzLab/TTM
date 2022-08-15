function [texture_mat]=save_all_pools_new(filename,savename, image, with_pixel_hist)
    tic
    descriptor = load(filename);
    textureDescriptor = descriptor.textureDescriptor;
    num_pools = length(textureDescriptor);
    poolingRegions = descriptor.poolingRegions;
    

    
    % initialize stats
    if with_pixel_hist == 1
        texture_mat = NaN(num_pools,3,2306);
    else
        texture_mat = NaN(num_pools,3,1790);
    end
    
    % digits_needed = length(num2str(num_pools(end)));
    for i=1:1
    % i=1;
    % adding zeros to make life easier later
    % assumes no more than 10,000 pooling region
        % pool_id = num2str(i,'%05.f');
        % disp(sprintf('pool %s',pool_id))
        texture_mat = extract_save_one_pool_flat(textureDescriptor{i},i,texture_mat,with_pixel_hist, 1);
    end
    toc
    
    % specify chunk size for compression
    h5create(sprintf('%s.h5',savename),sprintf('/%s',image),size(texture_mat));
    h5write(sprintf('%s.h5',savename),sprintf('/%s',image),texture_mat);
    
    
    names = ["isConstant_size_1_pyinds_0_1"
    "constval_size_1_pyinds_1_2"
    "pixelStats_size_6_pyinds_2_8"
    "pixelLPStats_size_10_pyinds_8_18"
    "autoCorrReal_size_245_pyinds_18_263"
    "autoCorrMag_size_784_pyinds_263_1047"
    "magMeans_size_18_pyinds_1047_1065"
    "cousinMagCorr_size_80_pyinds_1065_1145"
    "parentMagCorr_size_64_pyinds_1145_1209"
    "cousinRealCorr_size_320_pyinds_1209_1529"
    "parentRealCorr_size_256_pyinds_1529_1785"
    "varianceHPR_size_1_pyinds_1785_1786"
    "pixelHist_n_size_256_pyinds_1786_2042"
    "pixelHist_x_size_256_pyinds_2042_2298"
    "pixelHist_moments_size_4_pyinds_2298_2302"
    "cmu_size_1_pyinds_2302_2303"
    "cmat_size_3_pyinds_2303_2306"];
    
    
    fileID = fopen(sprintf('%s.txt',savename),'w');
    fprintf(fileID,'%s\n',names);
    fprintf(fileID,'%d,%d\n',poolingRegions(:,1),poolingRegions(:,2));
    fclose(fileID); 
end