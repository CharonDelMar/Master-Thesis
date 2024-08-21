
m = 3; % Dimension
num_samples = 10; % Number of sampling points
sigma = 0.5; % Standard deviation
range = [-1, 1]; % Allowed range for sampling values

samples = discrete_gaussian_sampling(m, num_samples, sigma, range);
disp(samples);

function samples = discrete_gaussian_sampling(m, num_samples, sigma, range)
    % m: Dimension of the space
    % num_samples: Number of sampling points
    % sigma: Standard deviation of the Gaussian distribution
    % range: Allowed range for sampling values (e.g., [-1, 1])

    % Initialize variables
    samples = zeros(num_samples, m);
    
    % Iterate to get the required number of samples
    for i = 1:num_samples
        while true
            % Sample m-dimensional vector
            sample = normrnd(0, sigma, 1, m);
            % Check if the sample is within the specified range
            if all(sample >= range(1) & sample <= range(2))
                samples(i, :) = sample;
                break;
            end
        end
    end
end
