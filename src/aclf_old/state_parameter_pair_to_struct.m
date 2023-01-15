function [x_theta_assignment] = state_parameter_pair_to_struct(x,theta,state_names,theta_names)
%STATE_PARAMETER_PAIR_TO_STRUCT Summary of this function goes here
%   Detailed explanation goes here

    % Input Processing
    if exist("state_names")
        if length(x) ~= length(state_names)
            error(['State has dimension ' num2str(length(x)) 'but state_names has length ' num2str(length(state_names)) '.' ])
        end
    end

    if exist("theta_names")
        if length(theta) ~= length(theta_names)
            error(['theta has dimension ' num2str(length(theta)) 'but theta_names has length ' num2str(length(theta_names)) '.' ])
        end
    end

    % Constants

    % Create field names
    state_tuples = {};
    for state_dim = 1:length(x)
        if exist("state_names")
            state_tuples{2*(state_dim-1)+1} = state_names{state_dim};
        else
            state_tuples{2*(state_dim-1)+1} = ['x' num2str(state_dim)];
        end
        state_tuples{2*(state_dim-1)+2} = x(state_dim);
    end
    x_assignment = struct(state_tuples{:});

    % Create field names for theta
    theta_tuples = {};
    for theta_dim = 1:length(theta)
        if exist("theta_names")
            theta_tuples{2*(theta_dim-1)+1} = theta_names{theta_dim};
        else
            theta_tuples{2*(theta_dim-1)+1} = ['theta' num2str(theta_dim)];
        end
        theta_tuples{2*(theta_dim-1)+2} = theta(theta_dim);
    end
    theta_assignment = struct(theta_tuples{:});

    x_theta_assignment = cell2struct( ...
        [struct2cell(x_assignment);struct2cell(theta_assignment)], ...
        [fieldnames(x_assignment);fieldnames(theta_assignment)]);

end

