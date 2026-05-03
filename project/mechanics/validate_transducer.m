function validate_transducer(p)
%VALIDATE_TRANSDUCER  Minimal parameter checks for the driver model.

    requiredP = {'Re','Le','Bl','Rms','Mms','Cms'};
    for k = 1:numel(requiredP)
        assert(isfield(p, requiredP{k}), 'Missing field p.%s', requiredP{k});
    end
end