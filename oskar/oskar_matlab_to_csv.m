freqs = 132:1:148;
ants = 1:1:256;

for freq = freqs
    fprintf("frequency: %f\n", freq)
    for iant = ants
        fprintf(".")
        load([num2str(freq) '/alphas_pola.mat']);
        alpha_te = alpha_te(iant,:,:);
        alpha_tm = alpha_tm(iant,:,:);
        x_te = squeeze(alpha_te(:,:));
        x_tm = squeeze(alpha_tm(:,:));
        suffix = [num2str(iant) '_' num2str(freq) '.txt'];
        dlmwrite(['element_pattern_spherical_wave_x_te_re_' suffix], real(x_te));
        dlmwrite(['element_pattern_spherical_wave_x_te_im_' suffix], imag(x_te));
        dlmwrite(['element_pattern_spherical_wave_x_tm_re_' suffix], real(x_tm));
        dlmwrite(['element_pattern_spherical_wave_x_tm_im_' suffix], imag(x_tm));

        load([num2str(freq) '/alphas_polb.mat']);
        alpha_te = alpha_te(iant,:,:);
        alpha_tm = alpha_tm(iant,:,:);
        y_te = squeeze(alpha_te(:,:));
        y_tm = squeeze(alpha_tm(:,:));
        dlmwrite(['element_pattern_spherical_wave_y_te_re_' suffix], real(y_te));
        dlmwrite(['element_pattern_spherical_wave_y_te_im_' suffix], imag(y_te));
        dlmwrite(['element_pattern_spherical_wave_y_tm_re_' suffix], real(y_tm));
        dlmwrite(['element_pattern_spherical_wave_y_tm_im_' suffix], imag(y_tm));
    end
    fprintf("\n")
end
