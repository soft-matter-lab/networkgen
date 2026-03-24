classdef assignmentmode
    properties

        auto = true;

        %%% Uniform assignment
        uniform = struct(...
        'value', 1);

        %%% Polydisperse assignment
        poly = struct();

        %%% Bimodal assignment
        bimodal = struct();

    end
end