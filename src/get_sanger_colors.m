function sc = get_sanger_colors

sc = [0.1 0.8 1;0 0 0;1 0 0;0.7 0.7 0.77;0.45 0.75 0;1 0.7 0.8];
sc = sc(floor(1:(1/16):6.95),:);
