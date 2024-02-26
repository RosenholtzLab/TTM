function save_all_pools(filename,savename)
    tic
    descriptor = load(filename);
    textureDescriptor = descriptor.textureDescriptor;
    num_pools = length(textureDescriptor);
    digits_needed = length(num2str(num_pools(end)));
    for i=1:2
    % i=1;
    % adding zeros to make life easier later
    % assumes no more than 10,000 pooling region
        pool_id = num2str(i,'%05.f');
        disp(sprintf('pool %s',pool_id))
        extract_save_one_pool(textureDescriptor{i},pool_id,savename);
    end
    toc
end