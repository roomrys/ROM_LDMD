function suffix1 = get_suffix1
global sim_number num_gauges num_times isScaled all_amr art_stab
strg = '';
if isScaled
    strg = strcat(strg, '(s)');
end
if all_amr
    strg = strcat(strg, 'A');
end
if art_stab
    strg = strcat(strg, '(AS)');
end
suffix1 = ['_s' num2str(sim_number) strg '_g' num2str(num_gauges) '_t'...
    num2str(num_times)];
