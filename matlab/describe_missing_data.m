function [data_index,number_of_observations,no_more_missing_observations] = describe_missing_data(data,gend,nvarobs)

[variable_index,observation_index] = find(~isnan(data));

data_index = cell(1,gend);
missing_observations_counter = NaN(gend,1);
for obs=1:gend
    idx = find(observation_index==obs);
    tmp = variable_index(idx);
    missing_observations_counter(obs,1) = nvarobs-length(tmp);
    data_index(obs) = { tmp(:) };
end
missing_observations_counter = cumsum(missing_observations_counter);

number_of_observations = length(variable_index);

if ~missing_observations_counter
    no_more_missing_observations = 0;
else
    tmp = find(missing_observations_counter>=(gend*nvarobs-number_of_observations));
    no_more_missing_observations = tmp(1);
end