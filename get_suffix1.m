function suffix1 = get_suffix1
global sim_number num_gauges num_times isScaled
if isScaled
suffix1 = ['_s' num2str(sim_number) '(s)' '_g' num2str(num_gauges) '_t'...
    num2str(num_times)];
else
suffix1 = ['_s' num2str(sim_number) '_g' num2str(num_gauges) '_t'...
    num2str(num_times)];
end
