input_data = rand(1000,3);
pack_1 = 1:300;
pack_2 = 301:600;
pack_3 = 601:1000;

input_data(pack_1,1) = rand(size(pack_1))*0.2;
input_data(pack_1,2) = rand(size(pack_1))*0.2;
input_data(pack_1,3) = rand(size(pack_1));

input_data(pack_2,1) = rand(size(pack_2))*0.2 + 0.3;
input_data(pack_2,2) = rand(size(pack_2));
input_data(pack_2,3) = rand(size(pack_2))*0.2 + 0.3;

input_data(pack_3,1) = rand(size(pack_3))*0.2 + 0.6;
input_data(pack_3,2) = rand(size(pack_3))*0.2 + 0.6;
input_data(pack_3,3) = rand(size(pack_3));

%% Visualize Input Data
figure(1);
hold off;
plot3(input_data(pack_1,1), input_data(pack_1,2), input_data(pack_1,3), 'ro');
hold on;
plot3(input_data(pack_2,1), input_data(pack_2,2), input_data(pack_2,3), 'go');
plot3(input_data(pack_3,1), input_data(pack_3,2), input_data(pack_3,3), 'bo');


%% Function PCA, input_data, out_dim
% Use this to switch methods
use_svd_method=1;

out_dim = 2;
[count in_dim] = size(input_data);

% SVD
if (use_svd_method)
    sub_input_data = (input_data - repmat(mean(input_data),count,1))/sqrt(count-1);
    [U,S,V] = svd(sub_input_data);
    % First out_dim columns as PCA bases
    pcaV = V(:,1:out_dim);
    output_data = input_data * pcaV;
else
    % EIG
    mean_input_data = mean(input_data);
    sub_input_data = input_data - repmat(mean_input_data, count,1);
    mean_mat = sub_input_data' * sub_input_data ./ (count - 1);
    cov_mat = mean_mat;
    [V D] = eig(cov_mat);
    % Last out_dim columns as PCA bases
    pcaV = V(:,in_dim - out_dim + 1: in_dim);
    output_data = input_data * pcaV;
end

%% End Function PCA

%% Visualize Output Data
figure(2);
hold off;
plot(output_data(pack_1,1), output_data(pack_1,2), 'ro');
hold on;
plot(output_data(pack_2,1), output_data(pack_2,2), 'go');
plot(output_data(pack_3,1), output_data(pack_3,2), 'bo');


